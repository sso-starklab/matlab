% postprocess_spikes            following manual sorting
%
% CALL                          rc = postprocess_spikes( filebase )
%
% GETS                          filebase        full path including file prefix
%                                               e.g. '/Volumes/slab1/mP23/dat/mP23_18/mP23_18'
%                                               the files required are mP23_18.eeg, mP23_18.clu.1, ...
%
% OPTIONAL ARGUMENTS (given as parameter/value pairs)
%
%                               doCheckClu      {1}, logical flag
%                               OverwriteParse  {-2}, Overwrite indicator for stimulus parsing
%                                                1: compute and overwrite
%                                                0: compute, do not overwrite (but do write if does not exist)
%                                               -1: load/compute if a file exists/doesn't (no writing at all)
%                                               -2: load if existing, compute and save if not (i.e. make sure exists)
%                               graphicsParse   {1}, plot a summary of all stimuli
%                                               
% RETURNS                       rc              vector of logical flags for each of the 7 processing steps (see below)
%
% ASSUMPTIONS:
%                               (1) manual clustering has been done
%                               (2) get_channel_setup has been updated properly
%                               (3) filebase is correct (see filebaseLookup; not necessary but useful)
% DOES
% First, preparations:
%           (1) make sure clu is backed up in man
%           (2) make sure an eeg file exists
%           (3) create an *prm.xml file using the information in get_channel_setup
%
% Then, analyses:
%           (1) reasligns res, prunes clu, and re-extracts spk and fet
%           (2) parse stimuli and determine levels for each type (parseNchannels)
%           (3) extract digital events from analog channel (parseDigital)
%           (4) compute all auto and cross correlations + significance (spikes2spikes)
%               and compute cluster quality metric (spikes_stats)
%           (5) segment into brain states (segmentBehavior)
%           (6) detect spontateous HFOs (find_all_hfos)
%           (7) compute spike-field statistics (spikePhase)
%           (8) generate stimulus PSTH (multipeth_make)
%
% CALLS                                            
%           (preparations) ParseArgPairs, LoadXml, create_prm_files
%           (1) reaslignres, uniqueres, dat2spkfet
%           (2) parseNchannels, get_stimlevels
%           (3) get_stimchans, parseDigital, uhist
%           (4) run_spikes_stats, spikes_stats_plot
%           (5) eeg2whl, get_emap, segmentBehavior
%           (6) find_all_hfos, plotProbeHFOs
%           (7) hfoAnalysisSpikingNew, hfoAnalysisSpikingCountNew
%           (8) spikePhase
%           (9) load_spikes, multipeth_make

% Note on minAmp selection
% If minAmp is not specified (left []), the routine parseOneChannel will
% receive an argument minAmpRelative of -0.0005. This means that anything
% below 0.1% of the max possible voltage (5V), i.e. 5 mV (equivalent to
% 0.5mA), will initially be treated as noise. Then the noise variance will
% be computed, and 3 SDs from the mean will be set as the threshold. This
% will yield, in most cases, a threshold of about 2 mV. 
%
% If one wishes to change this value, there are two options:
% (1) modify the definition of noise. this can be done by:
% minAmp                  = 0.002;    % 0.2% of the max, i.e. 10 mV (for 5V max)
%
% (2) circumvent the entire process and force a threshold value. this is done using a negative sign:
% minAmp                  = -0.01;    % threshold will be exactly 10 mV

% 
% to be added:
% (1) cell type determination by GWN (gwnAnalysis)
% (2) chirp analysis (wnAnalysis)
% (3) place field and phase precession (ppAnalysis)

% 10-jan-12 ES

% revisions
% 14-aug-19 ES+AL created function postprocess_spikes, based on postprocess_spk
% 04-sep-19 redundant if doSST removed
% 08-sep-19 uhist for parseDigital 'visualization' called with only rising edges (1)
% 12-sep-19 *lfp renaming in PC fixed
% 17-sep-19 toFlip added and defaulted
% 14-oct-19 modified check of mov
% 16-jan-20 added check for channel number ordering
% 14-may-20 added get_stimlevels also for RAMP 
% 15-may-20 (1) added support for segmentBehavior using *whl file (even if there are AM channels)
%           (2) added support for running find_all_hfos using a non-SWS state as a baseline 
% 25-jun-20 (1) swapped order of toflip: now, flipLFP is toflip(1), and flipSPK is toflip(2)
%           (2) changed call to plotProbeHFOs: now called with flipLFP and not with toflip
% 02-jul-20 (1) added selectionModeHFOs {'lfpPower'} and OverwriteSPWsearch {0}
%               can use 'csdPower' or 'csdSource'; if a re-run, recommended
%               to also modify OverwriteSPWsearch to be -2
% 05-jul-20 (1) added byParSST, defaults to 1 (see spikes_stats)
% 27-aug-20 (1) added doRealignRes, defaults to 1
%           (2) added doLinearTrack, defaults to 0
%           (3) added doCellType, defaults to 0

function rc = postprocess_spikes( filebase, varargin )

%------------------------------------------------------------------------
% constants

% conventions
shiftbits               = 2^15;                 % required for conversion from int16 to uint16 (in some animals, e.g. mC41, bugs in the system can be compensated by setting this to 0)
%toflip                  = [ 1 0 ];             % flip channels to conserve geometric ordering from channel 1 (bottom) to channel N (top) of each shank given that *xml is ordered in the same manner
spkNotDetrended         = 0;                    % typically, spikes are detrended during extraction, and the cell type classifiers assume this. If no detrending was done, this flag should be 1

% parsing parameters
parseSuffix                 = 'eeg';                                              % 'dat' is better temporally but slower and noisier
%minAmps                     = -0.0005;                                           % for most cases, can use -0.0005 (0.5 mA hard threshold); otherwise, can use 0.001 (1 mA, soft threshold)
%minAmps                     = [ -0.016 -0.0005 -0.0005 -0.019 -0.0005 ];                            % values for mP23_18
minDur                      = []; 
minDC                       = []; 
minCC                       = 0.45;

% HFO plotting parameters
ripFreq                     = [ 120 180 ];
ripPowSD                    = 6;
ripDurMS                    = [ 40 60 ];

% spike phase parameters
spectralParameters          = { 2 100 50 'log' };                               % old analyses used 20 bins between 2-40 Hz, linearly spaced

% PETH parameters
PSTH_stimVal                 = [ 0.001 5 ];                                                     % 1 mA to 5 A
PSTH_stimDurs                = [ 0.01 0.04; 0.04 0.08; 0.08 0.125; 0.125 0.25; 0.25 1 ];        % [s]
PSTH_uflag                   = 1;
PSTH_cmp                     = 'eq';
PSTH_stimTypes               = { 'PULSE', 'PSINE' };
PSTH_ilevels                 = { 'B', 'C' };
PSTH_sources                 = { 'LED', 'LD' };

%------------------------------------------------------------------------
% flow control
mfname                      = upper( mfilename );

% initialize output
rc                          = zeros( 1, 10 );
%[ 1-realsign; 2-Stim; 3-digital; 4-segmentBehavior; 5-LinearTrack; 6-SSTsource; 7-doHFO; 8-spikePhase; 9-Stim PSTH; 10-cellType ]

% argument handling
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ prmOW, doCheckClu ...
    , doParseStim, OverwriteParse, graphicsParse, minAmps ...
    , OverwriteGetLevels, graphicsGetLevels ...
    , doRealignRes, realignResMethod ...
    , doLinearTrack ...
    , doCellType ...
    , doDigital, OverwriteDIG ...
    , doSST, OverwriteSST, ignoredChannels, byParSST ...
    , doSB, OverwriteSB, graphicsSB, useWhlSB ...
    , doMOV, OverwriteMOV, graphicsMOV ...
    , doHFO, OverwriteHFOs, graphicsHFOs, bstateHFOs ...
    , selectionModeHFOs, OverwriteSPWsearch ...
    , doHFOspk, OverwriteHFOspk, OverwriteHFOphs ...
    , doPHASE, OverwritePHASE, graphicsPHASE ...
    , doMMP ...
    , toFlip, ilevel, savetype ...
    ]                       = ParseArgPairs(...
    { 'prmOW', 'doCheckClu' ...
    , 'doParseStim', 'OverwriteParse', 'graphicsParse', 'minAmps' ...
    , 'OverwriteGetLevels', 'graphicsGetLevels' ...
    , 'doRealignRes', 'realignResMethod' ...
    , 'doLinearTrack' ...
    , 'doCellType' ...
    , 'doDigital', 'OverwriteDIG' ...
    , 'doSST', 'OverwriteSST', 'ignoredChannels', 'byParSST' ...
    , 'doSB', 'OverwriteSB', 'graphicsSB', 'useWhlSB'...
    , 'doMOV', 'OverwriteMOV', 'graphicsMOV' ...
    , 'doHFO', 'OverwriteHFOs', 'graphicsHFOs', 'bstateHFOs' ...
    , 'selectionModeHFOs', 'OverwriteSPWsearch' ...
    , 'doHFOspk', 'OverwriteHFOspk', 'OverwriteHFOphs'...
    , 'doPHASE', 'OverwritePHASE', 'graphicsPHASE'...
    , 'doMMP' ...
    , 'toFlip', 'ilevel', 'savetype' ...
    }...
    , { -2, 1 ...
    , 1, -2, 1, -0.0005 ...
    , -2, 1 ...
    , 1, 'extremum' ...
    , 0 ...
    , 0 ...
    , 1, -2 ...
    , 1, -2, [], 1 ...
    , 1, -2, [ 1 1 ], 0 ...
    , 1, -2, 1 ...
    , 1, -2, -1, 'sws' ...
    , 'lfpPower', 0 ...
    , 1, -2, -2 ...
    , 1, -2, -1 ...
    , 1 ...
    , [ 1 0 ], 'B', 'png' ...
    }...
    , varargin{ : } );

% get the xml files
par                         = LoadXml( [ filebase '.xml' ] );
[ pathname, filename, extname ] = fileparts( filebase );
filename                   = [ filename extname ];
delim                      = strfind( filename, '_' );
if isempty( delim )
    delim                  = strfind( filename, '-' );
end
if ~isempty( delim )
    setup                  = filename( 1 : ( delim - 1 ) );
else
    setup                  = filename;
end
delim                      = strfind( setup, 'es' );
if ~isempty( delim )
    setup                  = setup( delim + 2 : end );
end
if ~exist( pathname, 'dir' )
    fprintf( '%s: missing filebase %s\n', upper( mfilename ), filebase )
    return
end
close all

flipLFP                     = toFlip( 1 );
flipSPK                     = toFlip( 2 );

%------------------------------------------------------------------------%
% verify that clu files exist and are backed up in man
%------------------------------------------------------------------------%
if doCheckClu
    
    ok                      = 0;
    clufiles                = dir( [ filebase '*clu*' ] );
    nspkGrps                = length( par.SpkGrps );
    if length( clufiles )   ~= nspkGrps
        fprintf( '%s: Missing/extra clu files - check that all *temp.clu* have been deleted!\n', upper( mfilename ) );
    else
        clufilebkp          = dir( [ pathname '/man/*clu*' ] );
        if isequal( [ clufiles.name ], [ clufilebkp.name ] )
            ok = 1;
        else
            fprintf( '%s: Apparently some clu files were not backed up or there are extra *clu* files in the man directory!\n', mfname )
        end
    end
    if ~ok
        fprintf( '\t\tCHECK and RERUN!!\n' )
        return
    end
  
end

%------------------------------------------------------------------------%
% make sure eeg file has an *eeg suffix, if not, rename the *lfp
%------------------------------------------------------------------------%
eegfname                    = [ filebase '.eeg' ];
lfpfname                    = [ filebase '.lfp' ];
if ~exist( eegfname, 'file' ) && exist( lfpfname, 'file' )                                          % rename lfp to eeg
    if isunix
        cmd                 = sprintf( '!mv %s %s', lfpfname, eegfname );
    elseif ispc
        cmd                 = sprintf( '!move %s %s', lfpfname, eegfname );
    else
        cmd                 = '';
    end
    eval( cmd );
end
if ~exist( eegfname, 'file' )
    fprintf( '%s: missing file %s\n', upper( mfilename ), eegfname )
    return
end

%------------------------------------------------------------------------%
% load/make prm/xml file
%------------------------------------------------------------------------%
par                         = create_prm_files( filebase, setup, prmOW );                           % requires an updated get_channel_setup.m
if prmOW > 1
    return
end
nGroups                     = length( par.SpkGrps );
for i                       = 1 : nGroups
    njumps                  = length( unique( diff( par.SpkGrps( i ).Channels ) ) ) - 1;
    if njumps > 0
        fprintf( '%s: channels of spike group %d not ordered properly\n', upper( mfilename ), i )
        return
    end
end

%------------------------------------------------------------------------%
% realign res by clusters
%------------------------------------------------------------------------%
if doRealignRes
    byParSST =0;
    % realign and re-extract
    rc1                     = realignres( filebase, 'method', realignResMethod, 'byPar', byParSST );
    rc2                     = resunique( filebase );
    nsamples                = par.SpkGrps( 1 ).nSamples;                                    % assume same for all spike groups
    peaksample              = par.SpkGrps( 1 ).PeakSample;                                  % assume same for all spike groups
    rc3                     = dat2spkfet( filebase, 'nsamples', nsamples, 'peaksample', peaksample );

    % summarize
    rc1                     = all( rc1( : ) );
    rc2                     = all( rc2( : ) );
    rc3                     = all( rc3( : ) );
    rc( 1 )                 = all( [ rc1 rc2 rc3 ] );
    
end

%------------------------------------------------------------------------%
% Stimulus statistics
%------------------------------------------------------------------------%
if doParseStim
       % for uLED
   uLED = 0;
   if uLED
            minAmps = [-0.0005/500 -0.0005/500 -0.0005/500 -0.0005/500 -0.0005/500 -0.0005/500 -0.0005/500 -0.0005/500 -0.0005/700 -0.0005/500 -0.0005/500 -0.0005/500];
            vcenters            = 1e-5 : 1e-5 : 6e-4;
            dcenters            = [ 1 2 5 10 20 60 100 200 400 800 1600 9600 ] / 1000;             % s
            fcenters            = [ 1 : 12 14 : 2 : 100 110 : 10 : 200 ];                          % Hz
            fmat                = [ 0 10; 0 40; 0 100; 0 200; 10 0; 40 0; 100 0; 200 0 ];          % [ Hz Hz ]
            ccenters            = sort( diff( fmat, [], 2 ) )';
            bins = cell( 1, 4 );
        
            bins{ 1 } = vcenters;
            bins{ 2 } = dcenters;
            bins{ 3 } = fcenters;
            bins{ 4 } = ccenters;
            OverwriteCategorize         = 1;
   end
    % parse
    stimchans                   = sort( get_stimchans( par ) );
    stims                       = parseNchannels( filebase, stimchans ...
        , 'Overwrite', OverwriteParse, 'suffix', parseSuffix, 'graphics', graphicsParse ...
        , 'minAmpRelative', minAmps, 'minDurationSEC', minDur, 'minDutyCycle', minDC ...
        , 'minCC', minCC );
    stim_plot( stims );
    
    % partition into discrete levels
    stimTypes                   = { 'PULSE', 'WN', 'ZAP', 'PSINE', 'SINE', 'RAMP' };
    for i                       = 1 : length( stimTypes )
        stimType                = stimTypes{ i };
        [ levels, stypes, nvals ]   = get_stimlevels( filebase ...
            , stimType, 'graphics', graphicsGetLevels, 'Overwrite', OverwriteGetLevels ,'durRange', [ 0.005 inf ] );
        fprintf( 1, '%s:\n', stimType )
        disp( [ stypes num2cell( [ nvals levels ] ) ] )
    end
    
    % summarize
    rc( 2 )                     = 1;
end

%------------------------------------------------------------------------%
% digital channels
%------------------------------------------------------------------------%
digChan                     = get_stimchans( par, [], 'digitalin' );
if doDigital && ~isempty( digChan )
    % for uLED
   %  shiftbits = 0;
    digfname                = sprintf( '%s.dig', filebase );
    dexists                 = exist( digfname, 'file' );
    if OverwriteDIG < 0 && dexists
        fprintf( 1, 'Loading existing parsed digital data for %s\n', filebase )
        load( digfname, '-mat', 'mat' );
    else
        mat                 = parseDigital( filebase, 'chan', digChan, 'shiftbits', shiftbits, 'Overwrite', OverwriteDIG );
    end
    
    [ aa, bb ]              = uhist( mat( mat( :, 3 ) == 1, 2 ) );
    fprintf( 1, 'Summary of digital events (event number, number of events):\n' )
    disp( [ bb' aa' ] )
    
    rc( 3 )                 = 1;
end

%------------------------------------------------------------------------%
% Segmentation into states based on movement and eeg
%------------------------------------------------------------------------%

% movment data 
if doMOV
    
    % compute
    mov                         = eeg2whl( filebase, 'Overwrite', OverwriteMOV, 'graphics', graphicsMOV );
    
    % post-hoc check that contents make sense
    if isempty( mov ) 
        fprintf( 1, 'Note - no Spotter channels defined\n' )
    elseif all( std( mov.pos ) < 1 )
        fprintf( 1, 'Note - sub-pixel variance in estimated position trace - check Spotter\n' )
    end
    
        
end

% segmentation into states based on am/movement and eeg
if doSB || doHFO || doPHASE
    
    try
        map                     = get_emap( filebase );
        neurochans              = map( :, 3 );
        neurochans              = neurochans( ismember( neurochans, setdiff( neurochans, ignoredChannels ) ) );
    catch
        neurochans              = sort( get_stimchans( par, [], 'neuronal' ) );
        fprintf( 1, 'Note - neuronal channels determined based on Anatomical groups in *xml file\n' )
    end
    if ~exist( [ filebase '.phs' ], 'file' )
        OverwriteSB             = 1;
    end
    if useWhlSB
        amchans                 = -1;
    else
        amchans                 = [];
    end
    segmentBehavior( filebase, 'neurochans', neurochans ...
            , 'amchans', amchans, 'Overwrite', OverwriteSB, 'graphics', graphicsSB );
    rc( 4 )                     = 1;
    
end

%------------------------------------------------------------------------%
% linear track and update stim
%------------------------------------------------------------------------%
if doLinearTrack
    
    % summarize
    rc( 5 )                     = 1;
    
end


%------------------------------------------------------------------------%
% Spiking statistics 
%------------------------------------------------------------------------%
% compute CCH, ACH, waveform stats, cell types
if doSST
    [ ~, sst ]          = run_spikes_stats( filebase, 'Overwrite', OverwriteSST, 'byPar', byParSST );
    rc( 6 )                 = 1;
    spikes_stats_plot( sst, filebase, 'savetype', savetype, 'graphics', 1 ...
        , 'flipLFP', flipLFP, 'flipSPK', flipSPK, 'hpfUsed', spkNotDetrended );                         % plot the units and save the figures
    
end


%------------------------------------------------------------------------%
% LFP oscillations and spike-field statistics
%------------------------------------------------------------------------%
if doHFO
    
    % detect and plot the HFOs
    try
        hstats              = find_all_hfos( filebase ...
            , 'Overwrite', OverwriteHFOs, 'graphics', graphicsHFOs ...
            , 'bstate', bstateHFOs, 'selectionMode', selectionModeHFOs, 'OverwriteSPWsearch', OverwriteSPWsearch ...
            , 'filtMode', 'dog', 'ignoredChannels', ignoredChannels );
        rc( 7 )             = 1;
    catch
        rc( 7 )             = 0;
    end 
    if graphicsHFOs && exist( [ filebase '.sps' ], 'file' )
        plotProbeHFOs( filebase, flipLFP, ignoredChannels ); % flipLFP is used here to flip the probe
    end

    % if ripples repeat the sst, now with the HFOs
    if ~exist( 'hstats', 'var' )
        load( [ filebase '.sps' ], '-mat', 'stats' )
        hstats = stats;
    end
    if any( inrange( hstats( :, 5 ), ripFreq ) ...
            | hstats( :, 7 ) >= ripPowSD ...
            | inrange( hstats( :, 9 ), ripDurMS ) )
        try
            load( [ filebase '.sst' ], '-mat', 'sst' )
            spikes_stats_plot( sst, filebase, 'savetype', savetype, 'graphics', 1 ...
                , 'flipLFP', flipLFP, 'flipSPK', flipSPK, 'hpfUsed', spkNotDetrended );
        catch
            fprintf( 1, 'Non-terminal error: could not load/plot the *sst\n' )
        end
    end
    
end

% spike-ripple analyses
if doHFO && doHFOspk
    
    % phase histograms
    hfoAnalysisSpikingNew( filebase, 'spontaneous' ...
        , 'Overwrite', OverwriteHFOspk, 'OverwritePHS', OverwriteHFOphs ...
        , 'figdir', 1, 'graphics', [ 0 1 0 0 ], 'ilevel', ilevel );
    
    % rate histograms
    hfoAnalysisSpikingCountNew( filebase, 'graphics', 1 ...
        , 'ilevel', ilevel, 'savetype', 'png' );
    
end

% spike phases relative to the multi-level LFP of the "best theta" channel:
if doPHASE
    try
        egroups             = 1 : length( par.SpkGrps );                                            % here consider ALL egroups, also those discarded by probe/neurochans
        load( [ filebase '.phs' ], 'eegchan', '-mat' );                                             % the 'best theta' channel
        for bstate          = { 'all', 'sws', 'the' }
            spikePhase( filebase, egroups( : ), 'ilevel', ilevel, 'phaseSource', eegchan...
                , 'periods', bstate, 'spectralMethod', spectralParameters ...
                , 'graphics', graphicsPHASE, 'Overwrite', OverwritePHASE, 'savetype', savetype );
        end
        rc( 8 )             = 1;
    catch
    end
end

%------------------------------------------------------------------------%
% PSTHs
%------------------------------------------------------------------------%

if doMMP
    %for uLED
    % PSTH_stimVal = [ 1e-5 6e-4 ];
    shanknums                       = [];                   % each shank stimulated by its own local source
    [ stimchans, ~, ~, source ]     = get_stimchans( par );
    idx                             = ismember( source, PSTH_sources );
    stimchans                       = stimchans( idx );
    
    for j                   = 1 : length( PSTH_ilevels )
        ilev                = PSTH_ilevels{ j };
        s                   = load_spikes( filebase, shanknums, ilev );
        for k               = 1 : length( PSTH_stimTypes )
            stimType        = PSTH_stimTypes{ k };
            for i = 1 : size( PSTH_stimDurs, 1 )
                stimDur = PSTH_stimDurs( i, : );
                for si   = 1 : length( stimchans )
                    astimchan = stimchans( si );
                    multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums, 'stimTypes', stimType, 'valRange', PSTH_stimVal, 'durRange', stimDur, 'channels', astimchan, 'uflag', PSTH_uflag, 'multi', PSTH_cmp );
                end
            end
        end
    end
    
    rc( 9 ) = 1;
    
end

%------------------------------------------------------------------------%
% cell type classification
%------------------------------------------------------------------------%
if doCellType

    % DCanalysis
    % mF108, F105
    % 4 blue LED, camk::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.046 0.048;0.048 0.052;0.046 0.05;0.048 0.05;0.048 0.05;0.048 0.05];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED','stimType', 'PSINE', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');

    % mP20
    % 4 blue LED, pv::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 0; % '0' effect on PV
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.046 0.048;0.048 0.052;0.046 0.05;0.048 0.05;0.048 0.05;0.048 0.05];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.01 0.04 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2');
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED','stimType', 'PSINE', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2');

    uflag = 0;
    multi = 'ge';
    simOnly =0;
    s = celltypeClassification( filebase,'slevel', slevel,'simOnly',simOnly,'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2', 'stimType', 'SINE' );
    
     
    % mP101, 
    % 4 blue LED, pv::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 0; % '0' effect on PV
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.04 0.05;0.05 0.052;0.02 0.03;0.048 0.05];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2');
    uflag = 0;
    multi = 'ge';
    simOnly =0;
    s = celltypeClassification( filebase,'slevel', slevel,'simOnly',simOnly,'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2', 'stimType', 'SINE' );

    % mA234,mS234 
    % 3 blue LED, CamKII::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.03 0.032;0.03 0.032;0.03 0.032];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');

    % mC400 
    % 12 blue uLED, CamKII::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    ilevel = 'B';
    slevel = [6.4e-4 7.5e-5; 6.4e-4 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-4;6.4e-5 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-5;6.4e-4 7.5e-5]; % ranges of pulse power
    slevel = [];
    s = celltypeClassification( filebase,'slevel', slevel,'graphics',1, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');
    
    % mF84, 
    % 1 blue LED, CamKII::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.020 0.022];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.125 0.250 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');
    
    % mC41, mF79, mF93
    % 5 blue LED, CamKII::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.009 0.0095;0.009 0.0095;0.009 0.0095;0.009 0.0095;0.009 0.0095];% mC41, mF93
    slevel = [0.050 0.052;0.050 0.052];% mF79
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');    

    % mDL5, mP23, mB142
    % 3 blue  (33,34,35), 2 red (37,38), PV::ChR2, PV::jaws
    supFlag = 0; % tag by increasing firing rate
    celltype= 0; % '0' effect on INT
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.025 0.03;0.038 0.04;0.018 0.02; 0.035 0.038];
    s = celltypeClassification( filebase,'slevel', slevel,'graphics',0, 'sourceType', {'LED','LD'}, 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2');
    s = celltypeClassification( filebase,'slevel', slevel,'graphics',1, 'sourceType', {'LED','LD'}, 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.08 0.125 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::chr2');

    supFlag = 1; % tag by decreasing firing rate
    slevel = [0.033 0.041;0.033 0.041];
    s = celltypeClassification( filebase,'slevel', slevel,'graphics',1, 'sourceType', 'LD', 'Overwrite',1, 'wavRange', [ 550 700 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::jaws');
    s = celltypeClassification( filebase,'slevel', slevel,'graphics',1, 'sourceType', 'LD', 'Overwrite',1, 'wavRange', [ 550 700 ], 'durRange', [ 0.250 1 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'pv::jaws');

    % mV99, 
    % 1 blue LED, VIP::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    celltype= 0; % '0' effect on PV
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.020 0.022];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'vip::chr2');
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'graphics', 0, 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.01 0.04 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'vip::chr2');

    % mK01, 
    % 1 blue LED, CCK::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    celltype= 0; % '0' effect on PV
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.020 0.022];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'cck::chr2');
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',0, 'wavRange', [ 400 500 ], 'durRange', [ 0.125 0.250 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'cck::chr2');

    % mDS1 mDS2
    % 2 blue LED, CamKII::ChR2
    supFlag = 0; % tag by increasing firing rate
    celltype= 1; % '1' effect on PYR
    ilevel = 'B';
    slevel = []; % ranges of pulse power
    slevel = [0.015 0.020;0.015 0.020];
    slevel = [0.006 0.007; 0.006 0.007];
    slevel = [0.015 0.025; 0.015 0.025];
    s = celltypeClassification( filebase,'slevel', slevel, 'sourceType', 'LED', 'Overwrite',0, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype',celltype, 'supFlag', supFlag, 'opsinType', 'camkii::chr2');
    
    rc( 10 ) = 1;

end


return

% EOF


%------------------------------------------------------------------------%
% August 2019
%------------------------------------------------------------------------%

% parsing modified extensively - simplified, supports non-linear mapping;
% does not support arbitrary templates yet.
% user control is collapsed into the argument minAmps, to be described in
% detail by AL
%
% for instance, for mP23_18, we used:
minAmps                     = [ -0.016 -0.0005 -0.0005 -0.019 -0.0005 ];
postprocess_spikes( filebase, 'minAmps', minAmps )

% to do - in this function:
% (1) parse: explain minAmps
% (2) write detailed documentation for this routine
% (6) spikes2spikes - consider adding check_mono

% optional - in this function:
% (1) digital:  should add a table

% to do - other routines:
% (1) get_channel_setup
% (2) filebaseLookup
