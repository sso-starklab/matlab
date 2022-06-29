% segmentBehavior       into SWS/REM/RUN/MOV states by theta eeg and accelerometer/video monitoring
%
% [ sts, statedurs ] = segmentBehavior( filebase, windur, amchans, eeg, Overwrite, graphics, savetype )
%
% gets:
%   filebase          adate/fnum pair or a true filebase (required)
%   windur            {2}; sec; determine the running window size to be used
%                       for all analyses (theta, theta-non-theta, am, movement)
%   amchans           {[]}; determines which channels are AM recording
%                       :If empty, will check for am channels in the *.prm.xml 
%                       :If no am channels, will check for a whl file and use that
%                       :(-1) forces whl file source (if such exists)
%   eeg               {[]}; eeg data or a channel number to use for eeg
%                       if empty will automatically determine the channel 
%                       with the highest theta/non-theta ratio and use that
%   Overwrite           : 1: compute and overwrite
%                       : 0: compute, do not overwrite (but do write if does not exist)
%                       :-1: load/compute if a file exists/doesn't (no writing at all)
%                       :{-2}: load if existing, compute and save if not (i.e. make sure exists)
%   graphics          {1}, of local results
%   savetype          {'png'}
% 
% returns:          
%   sts: structure with 6 fields:
%   	filebase: full path
%       generator: computer, date
%       sws := immobile, no-theta
%       run := mobile, theta
%       rem := immobile, theta
%       mov := mobile, no-theta
%   statedurs: for each of the states
%
% note
%   that this assume that a good theta signal is recorded, so it may
%   determine e.g. a spindle channel as theta for neocortical recordings, use
%   wisely.. (at some future point, should include a check for a theta peak
%   in the spectrum, and extend to cortical phenomena - spindles/HVS, up/down
%   states, etc)
%
% does
%   segments data according to theta signal and movement:
%   theta is taken from the channel with the highest theta/delta ratio (3/6 SD)
%   movement is taken either from am recordings (0.03/0.05 m/sec2) or from whl file (0.25/2 cm/sec)
%   segments are defined at the EEG sampling rate
%
% calls:
%   ParseArgPairs, verb, replacetok                     (general)
%   LoadXml, LoadStims, LoadVals, msave
%       mergevals, makesrslen, get_stimchans            (format)
%   bestTheta, thetaRatio,                              (lfp)
%   am2states, whl2states                               (movement)
%   parse, setdiffranges, uniteranges, 
%       intersectranges, resampleranges
%       getdatainranges, plotranges, dilutesegments     (sets)
%   fig_out                                             (graph)
%
% files:
%   required:             filebase.eeg, filebase.xml
%   optional:             filebase.val*, filebase.whl
%   intermediate output:  filebase.phs, filebase.am
%   output:               filebase.sts.sws, *.rem, *.run, *.mov
%
% 18-nov-12 ES

% revisions
% 11-dec-12 (1) windur externalized
%           (2) par structure call support
%           (3) val/evt file merging
% 27-dec-12 (1) minor bug fixes (mergeeg call, intersectranges)
%           (2) save intermediate states as well
% 21-jan-13 (1) fallback to no-AM fixed
% 07-apr-13 (1) thetaRatio blockalized and externalized
% 08-dec-14 (1) get amchans using get_stimchans.m
% 21-jul-15 (1) added check for corrupted *whl file
% 24-feb-18 (1) cleaned up
% 31-aug-18 (1) modified default amTH to be -1
% 10-sep-18 (1) modified default amTH to be [ 0.02 -1 ]
% 18-aug-19 cleaned up

% TODO:
% (4) exclude stim times during bestTheta computation                                   /not important if short stim/
% (6) presentation: scale the distribution and raw data to same ylim                    /better resolution if ylim differs/..

% DETAILED ALGORITHM:
% (1) Determine the channel with the "best" theta signal from all recorded
%   channels; define "best" by the overall theta-nontheta ratio. 
%   Save this channel number in a filebase.phs file (see phAnalysis for the format)
%   Use an external module "bestTheta" for this.
% (2) Use the selected channel for determining theta/delta ratio
%   and partition the data into theta/non-theta periods
%   Use an external module for this, "thetaRatio.m"
% (3) Determine the movement/immobility periods
%   -If amchans are defined by the amchans argument/get_date_parameters, 
%       Segment the data into movement/immobility based on the RMS of the AM 
%       (based on externally-defined absolute number; if using a data-based, 
%       there will be a threshold even in files in which the animal was 
%       immobile/mobile throughout..). Use an external routine am2states.
%   -If there are no am recordings but there is a whl file, use that - 
%       segment using an absolute threshold on the speed. Use an external
%       routine whl2states.m for this (use the am because during e.g. grooming the
%       animal is not moving but the head is; if there is really no theta then 
%       this would eventually be classified as "mov" rather than "run")
% (4) Combine the data from the two steps:
%       -immobility (no head acceleration/movement) and theta => REM
%       -immobility and no-theta => SWS
%       -mobility (head acceleration or horizontal movement) and theta => RUN
%       -mobility and non-theta => MOV (could be grooming, turning etc)
%       -others - undefined (could be a substantial portion)..
%   If there are neither whl nor am data: assume that in home cage and classify as
%       -theta => REM (true REM + in-cage running)
%       -non-theta => SWS (true SWS + e.g. awake grooming etc that does not generate theta)
% (5) Temporally smooth the epochs (in the literature its usually 1-2 sec; I'm
%       using 2 sec min duration and 0.5 sec min ISI). this is done with
%       an external module "dilutesegments" that first combines adjacent segments
%       then removes short ones
% (6) Save the outcome in 4 separate files (sws, rem, run, *.sts.mov)
%       and plot the summary in a *png (or other) figure

function [ sts, statedurs ] = segmentBehavior( filebase, varargin )

%----------------------------------------------------------------------%
% Initialize output
%----------------------------------------------------------------------%
sts                         = [];
statedurs                   = [];

%----------------------------------------------------------------------%
% Constants:
%----------------------------------------------------------------------%
% path to bash scripts:
process_merge               = '/usr/local/bin/mergeeeg';

% eeg theta selection + segmentation:
thetaBP                     = [ 5 11 ];                                                             % [Hz] theta band
nonthetaBP                  = [ 2 4 ];                                                              % [Hz] non-theta (delta) band

% eeg theta selection:
BT_Overwrite                = -2;                                                                   % use the existing theta selection
BT_specMode                 = 'welch';                                                              % 'mt' and 'timedomain' are just slower
BT_graphics                 = 1;                                                                    % plot
neurochans                  = [];                                                                   % which channels to focus on

% eeg segmentation:
npoles                      = 2;                                                                    % relevant for the phase estimation/timedomain
tdTH                        = [ 3 6 ];                                                              % non-theta/theta
tdMode                      = 'rms';                                                                % 'rms' or 'hilb'
minisiSECtd                 = 0.5;                                                                  % [s] minimal inter-event-interval
filtMode                    = 'iir';                                                                % 'iir' or 'fir';

% movement (am/whl) segmentation:
minisiSECmov                = 0.5;                                                                  % minimal inter-event-interval (sec)

% am segmentation:
%amTH                = [ 0.03 0.05 ];    % [m/s^2]; RMS of the vector sum; 0.03-0.05 should be used
amTH                        = [ 0.02 -1 ];                                                          % [m/s^2]; use data-driven segmentation; apparently 0.295/0.326 are nicely splitting the data (31-aug-18)
amOverwrite                 = 1;                                                                    % eventually switch to -2: overwrite any existing am file (am computed before 19-nov-12 were Z-scores, without vector summation)
amGraphics                  = [ 1 0 ];                                                              % do not plot
amRMS                       = 0.05;                                                                 % [s]; *am file resolution (segmentation is by windur though)

% whl segmentation:
whlTH                       = [ 0.25 2 ];      	                                                    % [cm/s];  speed TH 
whlOverwrite                = -2;                                                                   % do not overwrite if existing
whlGraphics                 = 0;                                                                    % logical; do not plot these instances
whlScale                    = 0.5;                                                                  % [cm/pixel]; modify this as needed (in Rutgers, camera was higher, 0.44)
whlRotate                   = [];                                                                   % [rad]; whl rotation
whlSpeedSD                  = [ 1 6 ];                                                              % [SD]; whl extreme speed removal
whlInterpolate              = 1;                                                                    % logical; whl interpolation

% general;
verbose                     = 1;
statenames                  = { 'sws', 'rem', 'run', 'mov', 'nul', 'imb', 'mob', 'not', 'the' };
sts2write                   = [ 1 : 4 6 : 9 ];
nstates                     = length( sts2write );
sts2plot                    = 1 : 5;

%----------------------------------------------------------------------%
% arguments
%----------------------------------------------------------------------%
nargs                       = nargin;
mfname                      = upper( mfilename );
if nargs < 1 || isempty( filebase )
    error( 'filebase required' )
end
[ windur, neurochans, amchans, eeg...
    , Overwrite, graphics, savetype ] = ParseArgPairs(...
    { 'windur', 'neurochans', 'amchans', 'eeg'...
    , 'Overwrite', 'graphics', 'savetype' }...
    , { 2, [], [], [], -2, 1, 'png' } ...
    , varargin{ : } );
tdWindow                    = windur;
mindurSECtd                 = windur;
mindurSECmov                = windur;
amWindow                    = windur;

% filebase
if isa( filebase, 'struct' ) && isfield( filebase, 'FileName' )
    par                     = filebase;
    filebase                = par.FileName;
elseif isa( filebase, 'char' ) || isa( filebase, 'cell' ) && length( filebase ) == 2
    par                     = LoadXml( filebase );
else
    fprintf( 1, '%s: incompatible format for 1st argument (filebase)\n', mfname )
    return    
end
[ pathname, filename, extname ] = fileparts( filebase );
if isempty( extname )
    ismerged                = 1;
else
    ismerged                = 0;
end
eegfname                    = [ filebase '.eeg' ];
whlfname                    = [ filebase '.whl' ];
if ~exist( pathname, 'dir' )
    fprintf( 1, 'missing source directory %s\n', pathname )
    return
end

% parameters
Fs                          = par.lfpSampleRate;
nchans                      = par.nChannels;
if isempty( amchans )
    amchans                 = get_stimchans( par, [], 'am' );
end
if ~isequal( amchans, -1 )
    amchans( amchans > nchans | amchans < 1 ) = [];
end
if isempty( eeg )
    eegchan                 = [];
elseif length( eeg ) == 1                                                                           % force using this eeg channel
    eegchan                 = eeg;
    eeg                     = [];
else                                                                                                % force using these eeg data
    eegchan                 = NaN;
end

% check whether to overwrite
stsfname                    = cell( 1, nstates );
stsexist                    = true( 1 );
for i                       = 1 : nstates
    stsfname{ i }           = [ filebase '.sts.' statenames{ sts2write( i ) } ];
    stsexist                = stsexist && exist( stsfname{ i }, 'file' );
end
if Overwrite < 0 && stsexist
    verb( sprintf( '%s: Loading existing state data for %s...'...
        , mfname, [ filename extname ] ), -verbose )
    sts.filebase            = filebase;
    info                    = dir( stsfname{ 1 } ); 
    sts.generator           = { 'unknown', datestr( info.date, 'ddmmmyy' ) };
    for i                   = 1 : nstates
        statename           = statenames{ sts2write( i ) };
        cmd                 = sprintf( 'sts.%s = load( ''%s'' );', statename, stsfname{ i } );
        eval( cmd );
        cmd                 = sprintf( 'statedurs( i ) = sum( diff( sts.%s, [], 2 ) + 1 ) / Fs;', statename );
        eval( cmd );
    end
    verb( 'Done!', verbose );
    return
else
    verb( sprintf( '%s: Segmenting %s into brain states'...
        , mfname, [ filename extname ] ), verbose )
end

% where to save summary figure
if length( graphics ) == 1
    graphics( 2 )           = BT_graphics;
end
if any( graphics ~= 0 )
    homedir                 = strfind( pathname, 'dat' ); 
    if isempty( homedir )
        homedir             = [ fileparts( pathname ) '/' ];
    else
        homedir             = pathname( 1 : homedir - 1 ); 
    end
    figdir                  = sprintf( '%sfigs', homedir );
    if ~exist( figdir, 'dir' )
        mkdir( homedir, 'figs' );
    end
    figname                 = sprintf( '%s/%s%s.sts', figdir, filename, extname );
end

% make sure eeg file exists
if ~exist( eegfname, 'file' ) && ismerged
    fprintf( 1, 'attempting to merge eeg files for %s...\n', filebase )
    if exist( process_merge, 'file' )
        pwd0                = pwd;
        eval( sprintf( 'cd( ''%s'' )', fileparts( pathname ) ) )
        eval( sprintf( '!%s %s', process_merge, filename ) )
        eval( sprintf( 'cd( ''%s'' )', pwd0 ) )
    end
end

% merge the val/evt files
if ~exist( eegfname, 'file' )
    fprintf( 1, 'missing eeg file %s\n', eegfname )
    return
elseif ismerged
    Vals                    = LoadStims( filebase );
    if isempty( Vals )
        Vals                = LoadVals( filebase );
    end
%     if isempty( Vals )
%         mergevals( par );
%     end
end
    
%----------------------------------------------------------------------%
% Computation I
% determine the channel with the highest theta/delta ratio 
% (should exclude stimulation times by val files, not critical though)
%----------------------------------------------------------------------%
if isempty( eegchan )
    verb( sprintf( '%s: Determining "best" theta channel...', mfname ), verbose )
    eegchan                 = bestTheta( filebase, par, neurochans, [], thetaBP, nonthetaBP, BT_Overwrite, graphics( 2 ), BT_specMode );
end

%----------------------------------------------------------------------%
% Computation II
% Segment file: theta/non-theta ratio by time-domain eeg filtering
%----------------------------------------------------------------------%
verb( sprintf( '%s: segmenting file by theta/delta ratio...', mfname ), verbose )

% compute the ratio
if ~isempty( eeg )
    tdRatioPlot             = thetaRatio( eeg,      eegchan, nchans, thetaBP, nonthetaBP, Fs, filtMode, npoles, tdWindow, tdMode );
else
    tdRatioPlot             = thetaRatio( eegfname, eegchan, nchans, thetaBP, nonthetaBP, Fs, filtMode, npoles, tdWindow, tdMode );
end
a                           = memmapfile( eegfname, 'Format', 'int16' );
neeg                        = length( a.Data ) / nchans;
clear a

% partition and apply post-hoc logic
verb( 'post-processing...', -verbose )
tmat0                       = parse( find( tdRatioPlot <= tdTH( 1 ) ) );
tmat1                       = parse( find( tdRatioPlot >= tdTH( 2 ) ) );
tmat0                       = dilutesegments( tmat0, mindurSECtd * Fs, minisiSECtd * Fs );          % no-theta
tmat1                       = dilutesegments( tmat1, mindurSECtd * Fs, minisiSECtd * Fs );          % theta
tmatNaN                     = setdiffranges( [ 1 neeg ], uniteranges( tmat0, tmat1 ) );             % neither
verb( 'done TD-ratio.', verbose )

%----------------------------------------------------------------------%
% Computation III
% Segment file by movement/no-movement (by whl file/am data) - optional
%----------------------------------------------------------------------%
if ~isempty( amchans ) && amchans( 1 ) ~= -1
    % movement state by AM data:
    verb( sprintf( '%s: segmenting file by AM...\n', mfname ), -verbose )
    [ mmat0, mmat1, am, amTH ] = am2states( filebase, 'channels', amchans, 'nchans', nchans ...
        , 'Fs', Fs, 'TH', amTH, 'graphics', amGraphics, 'Overwrite', amOverwrite ...
        , 'windur', amWindow, 'mindurSEC', mindurSECmov, 'minisiSEC', minisiSECmov, 'amT', amRMS );
elseif exist( whlfname, 'file' )
    try
        whl                 = load( whlfname );
        whlOK               = 1;
    catch
        whlOK               = 0;
    end
    if whlOK
        % movement state by WHL data:
        verb( sprintf( '%s: segmenting file by WHL...\n', mfname ), -verbose )
        [ mmat0, mmat1, mv ] = whl2states( filebase, whlOverwrite, whlGraphics...
            , whlTH, mindurSECmov, minisiSECmov...
            , whlScale, whlRotate, whlSpeedSD, whlInterpolate );
        % convert to EEG sampling rate
        mmat0               = resampleranges( mmat0, round( Fs / mv.Fs ), 1 );
        mmat1               = resampleranges( mmat1, round( Fs / mv.Fs ), 1 );
    end
else
    mmat0                   = [];
    mmat1                   = [];
end
mmatNaN                     = setdiffranges( [ 1 neeg ], uniteranges( mmat0, mmat1 ) );             % neither

%----------------------------------------------------------------------%
% Computation IV
% Combine and define brain states
%----------------------------------------------------------------------%
imb                         = mmat0;                                                                % immobile
mob                         = mmat1;                                                                % mobile
not                         = tmat0;                                                                % no-theta
the                         = tmat1;                                                                % theta
if ~isempty( mmat0 ) || ~isempty( mmat1 ) || ~isempty( amchans ) && amchans( 1 ) ~= -1
    sws                     = intersectranges( mmat0, tmat0 );                                      % SWS := immobile, no-theta
    rem                     = intersectranges( mmat0, tmat1 );                                      % REM := immobile, theta
    run                     = intersectranges( mmat1, tmat1 );                                      % RUN := mobile, theta
    mov                     = intersectranges( mmat1, tmat0 );                                      % MOV := mobile, no-theta
else                                                                                                % assume home cage, immobile
    sws                     = tmat0;                                                                % no theta (may include grooming, sniffing..)
    rem                     = tmat1;                                                                % theta (may include walking..)
    run                     = [];
    mov                     = [];
end
nul                         = setdiffranges( [ 1 neeg ], uniteranges( [ sws; rem ], [ run; mov ] ) );

%----------------------------------------------------------------------%
% Summarize
%----------------------------------------------------------------------%
t0                          = sum( diff( tmat0, [], 2 ) + 1 ) / Fs;
t1                          = sum( diff( tmat1, [], 2 ) + 1 ) / Fs;
tNaN                        = sum( diff( tmatNaN, [], 2 ) + 1 ) / Fs;

m0                          = sum( diff( mmat0, [], 2 ) + 1 ) / Fs;
m1                          = sum( diff( mmat1, [], 2 ) + 1 ) / Fs;
mNaN                        = sum( diff( mmatNaN, [], 2 ) + 1 ) / Fs;

statedurs                   = [ sum( diff( sws, [], 2 ) + 1 ) / Fs
    sum( diff( rem, [], 2 ) + 1 ) / Fs
    sum( diff( run, [], 2 ) + 1 ) / Fs
    sum( diff( mov, [], 2 ) + 1 ) / Fs 
    sum( diff( nul, [], 2 ) + 1 ) / Fs 
    sum( diff( imb, [], 2 ) + 1 ) / Fs 
    sum( diff( mob, [], 2 ) + 1 ) / Fs 
    sum( diff( not, [], 2 ) + 1 ) / Fs 
    sum( diff( the, [], 2 ) + 1 ) / Fs ];

verb( sprintf( '%s: TD durations: low/high/neither: %d/%d/%d sec'...
    , mfname, round( t0 ), round( t1 ), round( tNaN ) ), verbose )
verb( sprintf( '%s: MOV durations: immobile/mobile/neither: %d/%d/%d sec'...
    , mfname, round( m0 ), round( m1 ), round( mNaN ) ), verbose )
verb( sprintf( '%s: State durations: SWS/REM/RUN/MOV: %d/%d/%d/%d sec; null: %d sec'...
    , mfname, round( statedurs( 1 ) ), round( statedurs( 2 ) )...
    , round( statedurs( 3 ) ), round( statedurs( 4 ) ), round( statedurs( 5 ) ) ), verbose )

%----------------------------------------------------------------------%
% Graphics
%----------------------------------------------------------------------%
if abs( graphics( 1 ) )
    
    verb( sprintf( '%s: plotting summary...\n', mfname ), -verbose )

    tim                     = ( 1 : neeg )' / Fs;
    DSF                     = 32;
    didx                    = 1 : DSF : length( tim );
    COL0                    = [ 1 0 0 ];
    COL1                    = [ 0 0 1 ];
    colors                  = [ 1 0 0; 1 0.5 0; 0 0 1; 0 0.5 1; 0.3 0.3 0.3 ];
    LW                      = 3;
    
    fig = figure;

    %--------------
    % TD plotting - ratio, distribution
    subplot( 3, 1, 1 )
    plot( tim( didx ), tdRatioPlot( didx ) );
    xlim( tim( [ 1 end ] ) )
    hold on
    ypos                    = max( ylim );
    if ~isempty( tmat0 )
        ph                  = plotranges( ceil( tmat0 / Fs ), ypos * 0.9, [], 'linewidth', LW, 'color', COL0 );
    end
    if ~isempty( tmat1 )
        ph                  = plotranges( ceil( tmat1 / Fs ), ypos * 1, [], 'linewidth', LW, 'color', COL1 );
    end
    ylim( [ min( ylim ) ypos * 1.1 ] )
    hold off
    ylabel( 'T/D ratio' )
    xlabel( 'Time [s]' )
    title( replacetok( [ filename extname ], '\_', '_' ) )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    pos                     = get( gca, 'position' ); 
    set( gca, 'position', [ pos( 1 : 2 ) 0.5 pos( 4 ) ] )
    
    subplot( 3, 3, 3 )
    TH                      = tdTH;
    xx                      = tdRatioPlot;
    str = 'tdRatio';
    [ h, bins ]             = hist( xx, 1000 );
    hold on 
    line( h, bins, 'color', [ 1 1 1 ] * 0 ), 
    line( xlim, TH( 1 ) * [ 1 1 ], 'color', COL0 );
    if length( TH ) == 2
        line( xlim, TH( 2 ) * [ 1 1 ], 'color', COL1 );
    end
    title( sprintf( '%s; TH=%0.3g/%0.3g', str, TH( 1 ), TH( end ) ) );
    ylim( bins( find( ( cumsum( h ) / sum( h ) ) > 0.99, 1, 'first' ) ) * [ -0.01 1 ] + [ 0 eps ] )
    h0                      = hist( getdatainranges( xx, tmat0 ), bins );
    h1                      = hist( getdatainranges( xx, tmat1 ), bins );
    b0                      = barh( bins( : ), h0( : ) ); set( b0, 'edgecolor', COL0, 'facecolor', COL0 );
    b1                      = barh( bins( : ), h1( : ) ); set( b1, 'edgecolor', COL1, 'facecolor', COL1 );
    set( gca, 'xtick', [], 'xticklabel', [] )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    set( gca, 'position', [ 0.7 pos( 2 ) 0.2 pos( 4 ) ] )
    
    %--------------
    % AM/WHL plotting: RMS, distribution
    subplot( 3, 1, 2 )
    if exist( 'am', 'var' )
        movFs               = Fs;
        plot( tim( didx ), am( didx ) );
        xlim( tim( [ 1 end ] ) )
        ylabel( 'Accel [m/s^2]' )
    elseif exist( 'mv', 'var' )
        spd                 = mv.spd;
        movFs               = mv.Fs;
        timW                = ( 1 : length( spd ) ) / movFs;
        plot( timW, spd );
        xlim( timW( [ 1 end ] ) )
        ylabel( 'Speed [cm/s]' )
    end
    xlim0                   = xlim;
    hold on
    ypos                    = max( ylim );
    if ~isempty( mmat0 )
        ph                  = plotranges( ceil( mmat0 / Fs ), ypos * 0.9, [], 'linewidth', LW, 'color', COL0 );
    end
    if ~isempty( mmat1 )
        ph                  = plotranges( ceil( mmat1 / Fs ), ypos * 1, [], 'linewidth', LW, 'color', COL1 );
    end
    xlim( xlim0 );
    ylim( [ min( ylim ) ypos * 1.1 ] )
    hold off
    xlabel( 'Time [s]' )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    pos                     = get( gca, 'position' );
    set( gca, 'position', [ pos( 1 : 2 ) 0.5 pos( 4 ) ] )
    
    subplot( 3, 3, 6 )
    if exist( 'am', 'var' )
        TH                  = amTH;
        xx                  = am;
        str                 = 'Accel [m/s^2]';
    elseif exist( 'mv', 'var' )
        TH                  = whlTH;
        xx                  = spd;
        str                 = 'Speed [cm/s])';
    else 
        TH = [];
    end
    if ~isempty( TH )
        [ h, bins ]         = hist( xx( ~isnan( xx ) ), 1000 );
        hold on
        line( h, bins, 'color', [ 1 1 1 ] * 0 ),
        line( xlim, TH( 1 ) * [ 1 1 ], 'color', COL0 );
        if length( TH ) == 2
            line( xlim, TH( 2 ) * [ 1 1 ], 'color', COL1 );
        end
        title( sprintf( '%s; TH=%0.3g/%0.3g', str, TH( 1 ), TH( end ) ) );
        ylim( bins( find( ( cumsum( h ) / sum( h ) ) > 0.99, 1, 'first' ) ) * [ -0.01 1 ] )
        
        h0                  = hist( getdatainranges( xx, mmat0 ), bins );
        h1                  = hist( getdatainranges( xx, mmat1 ), bins );
        b1                  = barh( bins( : ), h1( : ) ); set( b1, 'edgecolor', COL1, 'facecolor', COL1 );
        b0                  = barh( bins( : ), h0( : ) ); set( b0, 'edgecolor', COL0, 'facecolor', COL0 );
        set( gca, 'xtick', [], 'xticklabel', [] )
        set( gca, 'box', 'off', 'tickdir', 'out' )
        set( gca, 'position', [ 0.7 pos( 2 ) 0.2 pos( 4 ) ] )
    end
    
    %--------------
    % summary - ranges, distribution
    subplot( 3, 1, 3 )
    ylim( [ 0.6 1 ] );
    hold on
    ypos                    = max( ylim );
    if ~isempty( sws )
        ph                  = plotranges( ceil( sws / Fs ), ypos * 0.7, [], 'linewidth', LW, 'color', colors( 1, : ) );
    end
    if ~isempty( rem )
        ph                  = plotranges( ceil( rem / Fs ), ypos * 0.8, [], 'linewidth', LW, 'color', colors( 2, : ) );
    end
    if ~isempty( run )
        ph                  = plotranges( ceil( run / Fs ), ypos * 0.9, [], 'linewidth', LW, 'color', colors( 3, : ) );
    end
    if ~isempty( mov )
        ph                  = plotranges( ceil( mov / Fs ), ypos * 1.0, [], 'linewidth', LW, 'color', colors( 4, : ) );
    end
    ylim( [ min( ylim ) ypos * 1.1 ] )
    xlim( tim( [ 1 end ] ) )
    hold off
    set( gca, 'ytick', [], 'yticklabel', [] )
    xlabel( 'Time [s])' )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    pos = get( gca, 'position' ); 
    set( gca, 'position', [ pos( 1 : 2 ) 0.5 pos( 4 ) ] )
    % add separator lines
    if ismerged 
        [ nsamplesfiles, fnums ] = makesrslen( filebase, 'eeg', -2 );
        filedurs            = nsamplesfiles / Fs;
        line( [ 1 1 NaN ]' * cumsum( filedurs( : ) ).', [ ylim NaN ]' * ones( 1, length( filedurs ) ), 'color', [ 0 0 0 ] );
        xx                  = cumsum( [ 0; filedurs ] );
        xpos                = ( xx( 1 : end - 1 ) + xx( 2 : end ) ) / 2;
        for i               = 1 : length( fnums )
            text( xpos( i ), ypos * 1.05, sprintf( '%d', fnums( i ) ), 'HorizontalAlignment', 'center' ); 
        end
    end
    
    subplot( 3, 3, 9 )
    hold on
    for i                   = sts2plot
        bh                  = bar( i, statedurs( i ) );
        set( bh, 'edgecolor', colors( i, : ), 'facecolor', colors( i, : ) );
    end
    set( gca, 'xtick', 1 : length( sts2plot ), 'xticklabel', statenames( sts2plot ) );
    xlim( [ 1 length( sts2plot ) ] + [ -0.5 0.5 ] )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    set( gca, 'position', [ 0.7 pos( 2 ) 0.2 pos( 4 ) ] )
    
end

% organize output structure
sts.filebase                = filebase;
sts.generator               = { computer, datestr( now, 'ddmmmyy' ) };
for i                       = 1 : nstates
    cmd                     = sprintf( 'sts.%s = %s;', statenames{ i }, statenames{ i } );
    eval( cmd );
end
    
%----------------------------------------------------------------------%
% Save
%----------------------------------------------------------------------%
if graphics( 1 ) > 0
    fig_out( fig, 1, [ figname '.' savetype ], savetype ),
    print( fig, [ '-d' savetype ], [ figname '.' savetype ] )
    verb( sprintf( '%s: Saved figure %s', mfname, [ figname '.' savetype ] ), verbose )
end

if Overwrite == 1 || ( Overwrite ~= -1 && ~stsexist )
    verb( sprintf( '%s: saving state segmented files...', mfname ), -verbose )
    for i                   = 1 : nstates
        cmd                 = sprintf( 'msave( ''%s'', %s );', stsfname{ i }, statenames{ sts2write( i ) } );
        eval( cmd );
    end
    verb( 'done!', verbose )
end

return

% EOF
