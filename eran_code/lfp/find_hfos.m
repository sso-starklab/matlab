% FIND_HFOS             detect HFO epochs in continuous data
%
% rips = find_hfos( filebase, chans, bperiods )
% find_hfos( ..., OverWrite, ripBP, ripTH, mindur, minisi, diffOrd, behavDur )
%
% filebase:     can be: -full path
%                       -a par structure
%                       -overloaded with a { date, fnum } cell;
%               in the latter case, a negative fnum indicates a merged directory
% chans:        can be a list of channels or
%               a negative number, indicating, an electrode group number
% bperiods:     -can be a list of time stamps (2-column matrix, sampling rate of
%               source file) from which to compute baseline (e.g. SWS)
%               -overloaded as an [-M -SD] pair for external specification of mean and SD for detection,
%               -or left empty, in which case the brain states (run/sws/rem/mov) will be
%               determined using movement data and theta/delta ratio
%
% other arguments are:
% Overwrite:    {1} to recompute and overwrite
%               0 to just compute (with writing but not overwriting)
%               -1 to just load/compute if a file exists/doesn't (no writing)
%               -2 load if exists, compute and save if doesn't
% ripBP:        {[80 250]}, [Hz] bandpass
% ripTH:        {[ 2 5 inf ]}, [SDs]; threshold, for event edges/peak/outliers
% mindur/minisi {0.015,0.015}, [sec]; minimal event duration/inter-event interval
% difford:      {0}; differntiation order; relevant for multi-channel input
%               only; 2 indicates CSD; 1 indicates difference; -x indicates
%               to use channel x as a reference for all others
% behavDur      {2}, [sec]; passed to segmentBehavior
%
% rips: structure with fields
%
%      filebase: full path
%         chans: channel/s numbers (1-based)
%        nsamps: [samples] at Fs
%     generator: cell array: { computer date }
%            Fs: [Hz]
%            bp: [Hz]; [ lowcutoff highcutoff ]
%            TH: [SD]; [ edges threshold outliers ] (SD)
%           RMS: [sec]
%        mindur: [sec]
%        minisi: [sec]
%             M: [mV] or arbitrary units; used for detection
%            SD: [mV] or arbitrary units; used for detection
%           seg: [samples]; above TH( 2 )
%             t: [samples]; time of peak power
%           pow: same units as M; mean power during seg
%            sd: [SD], mean over entire seg
%             f: [Hz], by number of troughs
%         edges: [samples]; above TH( 1 )
%         trigs: [samples]; t, aligned to closest trough
%         peaks: [index samples]; actually troughs
%
% calls         ParseArgPairs, LoadXml, LoadStims, LoadVals             (format)
%               MakeEvtFile                                             (evt format)
%               verb, num3str, clipmat, minmax                          (general)
%               segmentBehavior                                         (lfp)
%               makefir, firfilt, ma_rms, local_max, makegausslpfir     (ssp)
%               dilutesegments, getdatainranges, intersectranges, 
%                   isoverlap parse, resampleranges, setdiffranges
%                   uniteranges                                         (sets)
%               struct_select                                           (structs)
%
% see also      pt_avg (plots LFP/CSD/spectrum triggered on the HFOs)
%
% called by     find_all_hfos, get_hfo_power, get_hfo_times, plotHFOs

% 17-jul-12 ES

% revisions
% 11-dec-12 (1) merge vals when processing merged eeg files
%           (2) behavDur externalized
% 15-mar-13 if no SWS and no bperiods, issues a warning but uses the ENTIRE
%               file for baseline computations
% 30-mar-13 reference channel case fixed
% 19=may-13 bug in trigs clipping fixed
% 01-oct-13 algorithm changed. the following changes were made:
%           1. filtering. its now done (optionally) by difference-of-gaussian FIRs
%           (rather than a bandpass FIR or a sequence of IIR)
%           2. power computation. its now done (optionally) by low-pass filtering of the
%           rectified signal (same as RMS, but using a smoother filter).
%           3. baseline variance computation. its now computed from power
%           based on a clipped version of the filtered trace. Thus two
%           running power signals are computed (this reduces the bias of
%           ripple rate on ripple detection threshold)
%           -also, now all arguments are given as pairs
% 13-aug-19 cleaned a bit and renamed as find_hfos.m
% 18-aug-19 cleaned up

% DETAILED ALGORITHM:
% 1. partition data into states. use SWS for baseline computation,
%   excluding any periods during which there was stimulation.
%   alternatively use externally-defined mean and SD
% 2. pre-process the input signal. can use a single eeg channel, the 1st or 2nd
%   order spatial derivative (for removing a reference channel or computing
%   the CSD), and/or use any other source data (e.g. csd file, dat file...)
% 3. filter the data in the band-pass and compute a running RMS; average
%   the power over all channels and compute mean and SD for the SWS epochs
% 4. detect candidate HFOs (regions of high power)
% 5. apply post-hoc logic (combine adjacent regions, reject short ones;
%   expand edges; reject outliers)
% 6. align to min trough, compute statistics (amp, duration, freq, peak times)
%
% details of step 3:
% -the straight-forward method - band-pass filter, then compute power (RMS)
% by a MA, then take the stats from the entire SWS - has a couple of
% issues: (1) false positives - due to the filtering, which injects new
% oscillations (around any delta/step function-like artifact). This can be
% ameliorated by MA/median/gaussian filtering; gaussian has the best
% frequecny response and is used as the default; (2) false negatives - due
% to including the high-SD events (e.g. the real ripples) in the baseline
% stats, the number of detected events will be biased strongly by the
% amplitude and number of high-amp events. This can be minimized by not
% including high amp events in the baseline, and is applied by clipping (2
% SD) the band-pass before power computation. (3) noisy threshold crossing
% with a MA-based power computation - this is minimized by low-pass
% filtering with a Gaussian (instead of a box-car). 
%
% of course wavelets can be used instead. However they would have the same
% problem with false positives (delta function has flat spectrum), so they
% are not even included in this procedure.

function rips = find_hfos( filebase, chans, varargin )

% flow control
mfname                      = upper( mfilename );

% initialize output
rips                        = [];

%----------------------------------------------------------------------%
% PARAMETERS
%----------------------------------------------------------------------%

% constants
suffix                      = 'eeg';

% arguments
nargs                       = nargin;
if nargs < 2 || isempty( filebase ) || isempty( chans )
    rips                    = [];
    return
end
[ bperiods, Overwrite, ripBP, ripTH, mindurSECrip, minisiSECrip...
    , diffOrd, behaviorDurSEC, padBuffer, filtMode, powerMode, clipBase ...
    , savebase, verbose ]   = ParseArgPairs(...
    { 'bperiods', 'Overwrite', 'ripBP', 'ripTH', 'mindurSECrip', 'minisiSECrip'...
    , 'diffOrd', 'behaviorDurSEC', 'padBuffer', 'filtMode', 'powerMode', 'clipBase'...
    , 'savebase', 'verbose' }...
    , { [], -2, [ 80 250 ], [ 2 5 inf ], 0.015, 0.015...
    , 0, 2, [ -0.01 0.01 ], 'dog', 'LP', 1 ...
    , 'spw', 0 }...
    , varargin{ : } );

% file parameters
if isa( filebase, 'struct' ) && isfield( filebase, 'FileName' )
    par                     = filebase;
    filebase                = par.FileName;
elseif isa( filebase, 'char' ) || isa( filebase, 'cell' ) && length( filebase ) == 2
    par                     = LoadXml( filebase );
else
    fprintf( 1, '%s: incompatible format for 1st argument (filebase)\n', mfname )
    return
end
pathname                    = fileparts( filebase );
if ~exist( pathname, 'dir' )
    fprintf( 1, '%s: missing directory %s\n', mfname, pathname )
    return
end

% par file parameters
nBits                       = par.nBits;
nchans                      = par.nChannels;
P2P                         = par.VoltageRange;
switch suffix
    case { 'eeg', 'csd' }
        Fs                  = par.lfpSampleRate;
    case 'dat'
        Fs                  = par.SampleRate;
end
if chans( 1 ) < 0
    egroup                  = abs( chans( 1 ) );
    if egroup <= length( par.ElecGp )
        chans               = par.SpkGrps( egroup ).Channels + 1; % Consider only the channels used for spike sorting
    else
        fprintf( 1, '%s: Missing group %d\n', mfname, egroup )
        return
    end
else
    egroup                  = 0;
end

% ripple detection parameters:
ripWindow                   = pi / mean( ripBP );      % sec (0.01904 for default)
if length( ripTH ) == 1
    ripTH                   = [ 2 ripTH inf ];
elseif length( ripTH ) == 2
    ripTH                   = [ ripTH inf ];
end
ripTH                       = sort( ripTH );                  % SDs
EXPAND                      = [ -Fs Fs ];                    % post-hoc expansion; samples
filtMode                    = lower( filtMode );
switch filtMode
    case 'fir'
        hRip                = makefir( ripBP, Fs, [], 'bandpass' );
    case 'median'
        medWin1             = floor( Fs ./ ripBP( 1 ) / 2 ) * 2 + 1;
        hRip2               = makegausslpfir( ripBP( 2 ), Fs, 6 );
    case 'dog'
        hRip1               = makegausslpfir( ripBP( 1 ), Fs, 6 );
        hRip2               = makegausslpfir( ripBP( 2 ), Fs, 6 );
    otherwise
        error( 'unsupported filtMode' )
end
powerMode = upper( powerMode );
switch powerMode
    case 'RMS'
    case 'LP'
        %sd = Fs / ( 2 * pi * ( 1 / ripWindow ) ); x = -ceil( 6 * sd ) : ceil( 6 * sd ); gwin = 1/( 2 * pi * sd ) * exp( -(x.^2/2/sd.^2 ) ); powerWin = gwin / sum( gwin );
        powerWin            = makegausslpfir( 1 / ripWindow, Fs, 6 );
    otherwise
        error( 'unsupported powerMode' )
end

% referencing
if nchans == 1 || isempty( diffOrd ) || ~isa( diffOrd, 'double' )
    diffOrd                 = 0;
end
if length( diffOrd ) > 1
    diffOrd                 = diffOrd( 1 );
end
if diffOrd < 0 && ismember( abs( diffOrd ), chans )
    diffOrd                 = 0;
end

% output files
if egroup > 0
    schan                   = egroup;
else
    schan                   = chans( 1 );
end
savename1                   = sprintf( '%s.%s.%s', filebase, savebase, num3str( chans( 1 ), 3 ) );
if schan < 100 % assume max 999 elecs/elec groups
    savename2               = sprintf( '%s.evt.r%s', filebase, num3str( chans( 1 ), 2 ) );
else
    savename2               = sprintf( '%s.evt.%s', filebase, num3str( chans( 1 ), 3 ) );
end
if Overwrite < 0 && exist( savename1, 'file' )
    verb( sprintf( '%s: loading %s...', mfname, savename1 ), verbose )
    load( savename1, '-mat' )
    return
end

%----------------------------------------------------------------------%
% PREPARATIONS (eeg, val, segmentation)
%----------------------------------------------------------------------%

% make sure eeg and val files exist
if ~strcmp( suffix, 'eeg' ) && ~exist( [ filebase '.' suffix ], 'file' )
    fprintf( 1, '%s: %s file missing', mfname, suffix )
    return
end
eegfname                    = [ filebase '.' suffix ];
Vals                        = LoadStims( filebase );
a                           = memmapfile( eegfname, 'Format', 'int16' );
neeg                        = length( a.data ) / nchans;

% get baseline periods (SWS)
if isempty( bperiods ) % no baseline periods supplied
    
    % load sts.sws periods
    if ~exist( [ filebase '.sts.sws' ], 'file' )
        segmentBehavior( filebase, behaviorDurSEC );
    end
    verb( sprintf( '%s: loading sts.sws file...', mfname ), verbose )
    bperiods                = load( [ filebase '.sts.sws' ] );
    
    % now reject epochs with stimulation of any kind
    if ~exist( 'Vals', 'var' )
        Vals                = LoadVals( filebase );
    end
    if ~isempty( Vals )
        verb( sprintf( '%s: removing stimulus times...', mfname ), -verbose )
        vals                = resampleranges( Vals( :, 1 : 2 ), Fs, par.SampleRate );      % stim epochs
        pad                 = [ floor( padBuffer( 1 ) * Fs ) ceil( padBuffer( 2 ) * Fs ) ];
        vals                = [ vals( :, 1 ) + pad( 1 ) vals( :, 2 ) + pad( 2 ) ];
        bperiods            = setdiffranges( bperiods, vals );                         % sws w/o stim
        bperiods            = dilutesegments( bperiods, behaviorDurSEC * Fs, 0 );      % no merging, just remove short epochs
    end
    
    if isempty( bperiods )
        warning( '%s: NOTE: no SWS epochs!! using the entire duration for baseline computation!!!', upper( mfilename ) )
        bperiods            = [ 1 neeg ];
    end
    verb( sprintf( '%d/%d sec to be used for baseline.'...
        , round( sum( diff( bperiods, [], 2 ) + 1 ) / Fs ), round( neeg / Fs ) ), verbose )
    
elseif sum( bperiods( : ) < 0 ) > 1 && length( bperiods ) == 2
    
    % assume that bperiods is overloaded to be [ M SD ] computed externally
    % and do not recompute based on the EEG
    M                       = abs( bperiods( 1 ) );
    SD                      = abs( bperiods( 2 ) );
    
end
bperiods                    = intersectranges( bperiods, [ 1 neeg ] );

% load eeg
verb( sprintf( '%s: loading %s data...', mfname, suffix ), verbose )
if diffOrd < 0
    chans                   = [ chans( : ); abs( diffOrd ) ];
end
if diffOrd > 0 && length( chans ) < ( diffOrd + 1 )
    chans0                  = chans;
    % add extra channels
    % the way its written, it assumes that channels are consecutive, and
    % 1st derivative is channel minus deeper channel.
    % should change this to something like what I did in computeBaseline
    if diffOrd == 1
        echans              = chans - 1;
    else % ==2
        echans              = chans + [ -1 1 ];
    end
    probe                   = get_probe( par );                                     % make sure they are all neurochannels
    neurochans              = probe( ~isnan( probe ) );
    ridx                    = ~ismember( echans, neurochans );
    if sum( ridx )
        echans( ridx )      = [];
    end
    if diffOrd == 1 && isempty( echans )
        chans               = [ chans chans ];
    elseif diffOrd == 2 && length( echans ) == 1
        chans               = sort( [ chans chans echans ] );
    else
        chans               = sort( [ chans echans ] );
    end
end

if length( chans ) == 1
    eeg                     = a.data( chans : nchans : end );
    eeg                     = single( eeg );
    eeg                     = eeg / 2 .^ nBits * P2P;
else
    n                       = length( a.data );
    if diffOrd == 0
        eeg                 = zeros( n / nchans, 1, 'single' );
        for i               = 1 : length( chans )
            tmp             = a.data( chans( i ) : nchans : end );
            eeg             = eeg + single( tmp ) / 2 .^ nBits * P2P;
        end
        eeg                 = eeg / length( chans );
    else
        % an alternative if only the mean is required and we do not use blocks
        % this is NOT the same as first filtering (if median) and first power (in any case)
        idx                 = ( chans( : ) * ones( 1, n / nchans ) + ones( length( chans ), 1 ) * [ 0 : nchans : ( n - nchans ) ] )';
        eeg                 = single( a.data( idx ) ) / 2 .^ nBits * P2P;
    end
    
end
clear a

%----------------------------------------------------------------------%
% DETECT RIPPLE EVENTS by time-domain filtering and TH:
%----------------------------------------------------------------------%
verb( sprintf( '%s: detecting HFOs (%d channels)...', mfname, length( chans ) ), -verbose )
neeg                        = size( eeg, 1 );

% compute difference if requested:
if diffOrd == 1                                                                     % compute spatial derivative
    eeg                     = diff( eeg, 1, 2 );
    chans                   = chans0;
elseif diffOrd == 2                                                                 % compute csd
    eeg                     = -diff( eeg, 2, 2 );
    chans                   = chans0;
elseif diffOrd < 0                                                                  % remove a reference channel
    ref                     = chans == abs( diffOrd );
    eeg                     = eeg( :, ~ref ) - eeg( :, ref ) * ones( 1, length( chans ) - 1 );
    chans                   = setdiff( chans, abs( diffOrd ) );
end

% algorithm:
% 1. band pass filter using DOG (lowpass at the higher cut-off, then at the lower cut-off, then subtract)
switch filtMode
    case 'fir'
        rip                 = firfilt( eeg, hRip );
    case 'median'
        eegNotHigh          = firfilt( double( eeg ), hRip2 );
        eegLo               = medfilt1( eegNotHigh, medWin1 );
        rip                 = single( eegNotHigh - eegLo );
    case 'dog'
        eegNotHigh          = firfilt( eeg, hRip2 );
        eegLo               = firfilt( eegNotHigh, hRip1 );
        rip                 = eegNotHigh - eegLo;
end
if clipBase
    clipTH                  = ripTH( 2 );
    mmx                     = mean( rip, 1 )' * [ 1 1 ] + clipTH * std( rip, [], 1 )' * [ -1 1 ];
    ripC                    = rip;
    for i                   = 1 : size( rip, 2 )
        ripC( :, i )        = clipmat( rip( :, i ), mmx( i, 1 ), mmx( i, 2 ) );
    end
end
switch powerMode
    case 'RMS'
        ripPower0           = ma_rms( rip, round( Fs * ripWindow ) );
        if clipBase
            ripPower1       = ma_rms( ripC, round( Fs * ripWindow ) );
        end
    case 'LP'
        ripPower0           = firfilt( abs( rip ), powerWin );
        if clipBase
            ripPower1       = firfilt( abs( ripC ), powerWin );                     % LP(|clip(BP(x))|)
        end
end
chan                        = 1;
if size( eeg, 2 ) > 1
    rip                     = mean( rip, 2 );                                       % average over channels
    ripPower0               = mean( ripPower0, 2 );
end
if ~exist( 'M', 'var' ) && ~exist( 'SD', 'var' )
    if clipBase
        ripPowerSWS         = getdatainranges( ripPower1, bperiods );               % exclude the high-amp events from the mean
    else
        ripPowerSWS         = getdatainranges( ripPower0, bperiods );
    end
    if isempty( ripPowerSWS )                                                       % if there is no SWS in the file, cannot detect
        M                   = inf;
        SD                  = inf;
    else
        M                   = mean( ripPowerSWS, 1 );
        SD                  = std( ripPowerSWS, 1 );
    end
end

% detect candidate events
TH                          = M( chan ) + ripTH( 1 : 2 ) * SD( chan );
rmat                        = parse( find( ripPower0 > TH( 2 ) ) );

% apply post-hoc temporal constraints to ripple epochs:
rmat                        = dilutesegments( rmat, mindurSECrip * Fs, minisiSECrip * Fs );
ridx                        = sum( ( ( rmat + EXPAND( 1 ) ) < 1 ) + ( ( rmat + EXPAND( 2 ) ) > neeg ), 2 ) > 0;
rmat( ridx, : )             = [];
nrips                       = size( rmat, 1 );

%----------------------------------------------------------------------%
% compute statistics, peaks, and expand edges
%----------------------------------------------------------------------%
tpeak                       = zeros( nrips, 1 );
amp                         = zeros( nrips, 1 );
f                           = zeros( nrips, 1 );
rmate                       = zeros( nrips, 2 );
trigs                       = zeros( nrips, 1 );
allpeaks                    = [];
for i                       = 1 : nrips
    widx                    = rmat( i, : ) + EXPAND;
    % determine time of peak and mean power in that range
    pow                     = ripPower0( rmat( i, 1 ) : rmat( i, 2 ), chan );
    peakidx                 = round( sum( ( 1 : length( pow ) )' .* pow ) / sum( pow ) );
    tpeak( i, : )           = rmat( i, 1 ) + peakidx - 1;
    amp( i, : )             = mean( pow );
    % expand edges to lower TH (this approx. doubles the duration for 2/5 SDs)
    mat                     = parse( find( ripPower0( widx( 1 ) : widx( 2 ), chan ) >= TH( 1 ) ) );
    vec                     = Fs + 1 + [ 0 diff( rmat( i, : ) ) ];
    sidx                    = isoverlap( mat, vec );
    if sum( sidx ) == 1
        rmate( i, : )       = mat( sidx, : ) + widx( 1 ) - 1;
    else
        rmate( i, : )       = uniteranges( vec, mat( sidx, : ) ) + widx( 1 ) - 1;
    end
    % detect all peaks (actually troughs..) in the range and determine event frequency
    idx                     = rmat( i, 1 ) : rmat( i, 2 );
    rr                      = rip( rmate( i, 1 ) : rmate( i, 2 ), chan );
    peaks                   = local_max( rr, 'min' );
    f( i, : )               = Fs / mean( diff( peaks ) );
    % edit the occurence time to the closest trough
    tpeaks                  = peaks + rmate( i, 1 ) - 1;
    [ ~, minidx ]           = min( abs( tpeaks - tpeak( i ) ) );
    if isempty( minidx )% || length( peaks ) == 1
        fprintf( '%s: PROBLEM!!\n', upper( mfilename ) )
        trigs( i )          = NaN;
    else
        trigs( i )          = tpeaks( minidx );
        allpeaks            = [ allpeaks; [ i * ones( length( peaks ), 1 ) tpeaks ] ];
    end
end

% reject outliers by highest (3rd) TH (should be a very high number)
sd                          = ( amp - M( chan ) ) / SD( chan );
rmv                         = isnan( f ) | sd > ripTH( 3 ) | isnan( trigs );
if sum( rmv )
    tpeak( rmv )            = [];
    amp( rmv )              = [];
    sd( rmv )               = [];
    f( rmv )                = [];
    rmat( rmv, : )          = [];
    rmate( rmv, : )         = [];
    trigs( rmv )            = [];
    nrips                   = nrips - sum( rmv );
    allpeaks( ismember( allpeaks( :, 1 ), find( rmv ) ), : ) = [];
end
tmat                        = [ find( ~rmv ) ( 1 : sum( ~rmv ) )' ];
for i                       = 1 : size( tmat, 1 )
    allpeaks( ismember( allpeaks( :, 1 ), tmat( i, 1 ) ), 1 ) = tmat( i, 2 );
end
verb( sprintf( '%s: %d events detected (%d sec).', mfname, nrips, ceil( neeg / Fs ) ), verbose )

%----------------------------------------------------------------------%
% summarize and save
%----------------------------------------------------------------------%

% summarize
rips.filebase               = filebase;
rips.chans                  = chans;
rips.nsamps                 = neeg;
rips.generator              = { computer, datestr( now, 'ddmmmyy' ) };
rips.Fs                     = Fs;
rips.bp                     = ripBP;
rips.TH                     = ripTH;
rips.RMS                    = ripWindow;
rips.mindur                 = mindurSECrip;
rips.minisi                 = minisiSECrip;
rips.filtMode               = filtMode;
rips.powerMode              = powerMode;
rips.clipBase               = clipBase;
rips.M                      = M( chan );
rips.SD                     = SD( chan );
rips.seg                    = rmat;                     % samples
rips.t                      = tpeak;                    % samples
rips.pow                    = amp;                      % A2D units
rips.sd                     = sd;                       % SDs
rips.f                      = f;                        % Hz
rips.edges                  = rmate;                    % samples
rips.trigs                  = trigs;                    % samples
rips.peaks                  = allpeaks;                 % actually, troughs

% detect overlapping events
idx                         = find( ( rips.edges( 2 : end, 1 ) - rips.edges( 1 : end - 1, 2 ) ) < 0 );
idx                         = unique( [ idx; idx + 1 ] );
mat                         = parse( idx );
% now go over overlapping segments. for each, determine the global
% properties. keep those in say the first index and mark the others for
% deletion. afterwards delete all but the first in each segment
% note that this may be due to repetitive stimulation, and then the
% expansion is not really adequate. 
irmv                        = [];
newpeaks                    = [];
for i                       = 1 : size( mat, 1 )
    i0                      = mat( i, 1 );
    cidx                    = ( mat( i, 1 ) : mat( i, 2 ) )';
    irmv                    = [ irmv; cidx( 2 : end ) ];
    % determine time of peak and mean power in that range
    segs                    = rips.seg( cidx, : );
    [ ~, maxidx ]           = max( diff( segs, 1, 2 ) );
    seg                     = rips.seg( cidx( maxidx ), : );
    edges                   = minmax( rips.edges( cidx, : ) );
    pow                     = ripPower0( edges( 1 ) : edges( 2 ), chan );
    peakidx                 = round( sum( ( 1 : length( pow ) )' .* pow ) / sum( pow ) ); % center of mass of entire (expanded) region
    tpeak                   = edges( 1 ) + peakidx - 1;
    amp                     = mean( pow );
    sd                      = ( amp - M( chan ) ) / SD( chan );
    % detect all peaks (actually troughs..) in the range and determine event frequency
    rr                      = rip( edges( 1 ) : edges( 2 ), chan );
    peaks                   = local_max( rr, 'min' ); % may require that each trough surpasses the lower TH (may add stability)
    f                       = Fs / mean( diff( peaks ) );
    % edit the occurence time to the closest trough
    tpeaks                  = peaks + edges( 1 ) - 1;
    [ ~, minidx ]           = min( abs( tpeaks - tpeak ) );
    trigs                   = tpeaks( minidx );
    newpeaks                = [ newpeaks; [ i0 * ones( length( peaks ), 1 ) tpeaks ] ];
    % update:
    rips.edges( i0, : )     = edges;
    rips.seg( i0, : )       = seg;
    rips.pow( i0 )          = amp;
    rips.sd( i0 )           = sd;
    rips.t( i0 )            = tpeak;
    rips.f( i0 )            = f;
    rips.trigs( i0 )        = trigs;
end

if size( mat, 1 ) > 0
    % remove the unnecessary events (the ones not updated)
    kidx                    = true( size( rips.edges( :, 1 ) ) );
    kidx( irmv )            = 0;
    rips                    = struct_select( rips, kidx );
    % update the peaks: remove the old peaks
    rmv                     = ismember( rips.peaks( :, 1 ), idx );
    rips.peaks( rmv, : )    = [];
    % add the new ones
    allpeaks                = rips.peaks;
    allpeaks                = sortrows( [ allpeaks; newpeaks ], 1 );
    tmat                    = [ unique( allpeaks( :, 1 ) ) ( 1 : length( rips.t ) )' ];
    for i                   = 1 : size( tmat, 1 )
        allpeaks( ismember( allpeaks( :, 1 ), tmat( i, 1 ) ), 1 ) = tmat( i, 2 );
    end
    rips.peaks              = allpeaks;
end
nrips                       = length( rips.trigs );

% save evt and sws files
if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( savename1, 'file' ) )
    verb( sprintf( '%s: saving %s...', mfname, savename1 ), verbose )
    save( savename1, 'rips', '-v6' )
    if strcmp( savebase, 'spw' ) % otherwise do not make an event file here!
        MakeEvtFile( rips.t(:), savename2, 'hfo', Fs, 1 );
    end
end
verb( sprintf( '%s: done.', mfname ), verbose )

return

% EOF
