% hfoAnalysisSpiking            compute LFP phase and spiking during HFOs
%
% [ uphs tims s saux ] = hfoAnalysisSpiking( filebase, trigMode )
% 
% returns
% uphs: unwrapped (multi-cycle) phases
% tims: sample for each phase (@spkFs)
% s:    information about spikes during events. 
%           fields include: clu/res/map (as in load_spikes.m), plut:
%           phs - multi-cycle phase of each spike
%           bins - multi-cycle phase bins 
%           phsbin - multi-cycle, multi-event binned spikes
%       To recover the event numbers, write
%           floor( s.phsbin / length( s.bins ) )
%       To recover the within-event bin, write:
%           rem( s.phsbin - 1, length( s.bins ) ) + 1
%       this can be used to compute phase histograms, preferred phases,
%           phase-resolved cross-correlation and so on
% saux:  some quick results done internally. include multi-cycle phase,
%           single-cycle phases, preferred phases + p-values, and
%           single-unit-all-event/single-event-all-unit rasters. 
% Does:
% 1. determines which channel to use for the LFP phase
% 2. determines which events to compute the phase for
% 3. computes the phase and unwrap (see eegPhase.m, calcPhase.m)
% 4. saves the phase information to disk
% 5. determines which units to assign phases to
% 6. assigns phases to units
% 7. computes phase-resolved raster plots and phase-histograms. Also collapse
%       histograms to a single cycle. 
% 8. summarize - statistics and graphics for the LFP phases (1 plot)
%    statistics and graphics for the spiking (3 plots)
% 9. save to disk
%
% Defaults:
% 1. the channel on which the 'best' spontaneous HFOs occur
% 2. all spontaneous events that do not overlap with any external stimulation
% 3. use a hilbert transform on the piece-wise upsampled eeg
% 4. save to the dat directory
% 5. all well-isolated units (on all shanks)
% 6. use 20 phase bins (adequate resolution, about 0.3 ms for ripples)
% 7. Does not compute cross-correlations,
%       but does return the spike phases for subsequent
%       routines (e.g. PRCC)
% 8. plot the gain, use a uni-modal statistical test
% 9. save to the mat/figs directories
%
% Notes:
% (1) To prevent loading the EEG/DAT data multiple times, its best to
% modularize the use of this routine. A first run should be done with
% trigMode 'all', but without any spike processing. This will include
% segments from (1) same-channel/no-stim HFOs (2) same-shank stim (not
% necessarily unique). Alignment is by ~peak trough for HFOs and by stim
% onset for induced events. The routine then saves the phases and can
% exist.
% Subsequent calls will attempt to use the pre-saved phases. If a subset of
% periods is not available in an already-saved *phs* file, it will compute
% the phases for those but not save (unless an OverwritePHS of -3 is
% defined, in which case the missing segments will be added to the existing
% file)
% (2) This routine is general enough to be used for e.g. theta phase in the
% place field as well. Should be called with (1) proper events - e.g. the
% run part of each lap (2) the proper band-pass - e.g. theta frequency. The
% only issue is that then the data may become very large (when spikes are at
% spkFs), so it might in fact be better to keep the original eeg sampling
% rate and down-sample the spikes.
%
% files:
% required - *eeg, *xml, *clu, *res
% intermediate - *phs.###_###.###, where the # stand for lowpass_highpass
%       [Hz] and channel number (3 digits). Other parameters are only specified
%       within the file.
% 
% calls:
% (blab)          LoadXml
% (circs)         circ_mean, ray_test
% (formats)       get_egroup, get_emap, get_stimchans, get_triggers, LoadVals, LoadStims
% (general)       clipmat, minmax, num3str, ParseArgPairs, pol2car, remrnd, replacetok, scaleto, uhist, verb
% (graph)         alines, barwerror, circ, fig_out, myjet, textf, tilefig
% (hfoAnalysis)   hfoAnalysisSpikingPlot1
% (lfp)           eegPhase, plotHFOs, pt_avg
% (spikes)        determine_units, get_spikes, plot_raster
% (ssp)           fft_upsample, makegaussfir, makegausslpfir, makefir
% (stats)         bounds, calc_sem, calc_spearman
% (sets)          intersectranges, inranges, isoverlap, parse, resampleranges, setdiffranges, sortranges
%
% see also          hfoAnalysisSpikingCount

% 26-jul-13 ES

% revisions
% 28-jul-13 (1) compute phase, unwrapped phase, and cycle number by
%               eegPhase. use subroutines calcPhase, unwrapPhase, and
%               monotonic (see eegPhase.m)
%           (2) summmarize the statistics for the cycle number and phase
%                   distribution (which is not uniform)
% 29-jul-13 (1) save phase file
%           (2) spiking analysis
% 02-oct-13 (1) modified a bit to default to DOG filtering for phase
% 30-oct-13 (1) for spontaneous, remove periods during stim. 
%               for all, compute BL rates based on non-stim times
%           (2) expanded output to enable accumulation
%           (3) added eprops, for spontaneous this is 
%               [ SD f dur pow ]
% 17-nov-13 (1) major modifications..
% 12-jan-14 updates for "external" mode
% 17-aug-19 cleaned a bit and renamed hfoAnalysisSpikingNew
% 09-mar-21 cleaned properly and renamed back as hfoAnalysisSpiking
% 11-mar-21 additional cleaning and documentation
%           removed dependency on makesrslen
% 24-mar-21 modified graphical callback routine

function [ stats, s, uphs, tims, strs ] = hfoAnalysisSpiking( filebase, trigMode, varargin )

uphs                                = [];
tims                                = [];
s                                   = [];
stats                               = [];
fig                                 = [];

%------------------------------------------------------------------------
% constants
%------------------------------------------------------------------------
% flow control
verbose                             = 1;

% histograms
maxDurSEC                           = 0.5;                                  % remove super-long spontaneous events 
cycSupport                          = 0.99;                                 % compute histograms for 99% of the cycles
cycRange                            = [ -10 10 ];                           % but never more than 10 to each side

% graphics
colors                              = [ 0 0 0.7; 1 0 0 ];
blackColor                          = [ 0 0 0 ]; 
sepStyle                            = '--';
ng                                  = 4;

%------------------------------------------------------------------------
% arguments
%------------------------------------------------------------------------
nargs                               = nargin;
if nargs < 1 || isempty( filebase )
    error( 'missing argument filebase' ); 
end
if nargs < 2 || isempty( trigMode )
    trigMode                        = 'induced';
end
trigMode = lower( trigMode );
if ~ismember( trigMode, { 'induced', 'spontaneous', 'external', 'all' } )
    error( 'erroneous specification of trigMode' )
end
[ suffix, freqBP, filtmode, phsmode, OverwritePHS, nbins, iperiods ...
    , stimType, wavRange, sourceType, durRange, valRange, uStim, cmp, simOnly, minEvents ...
    , channel, periods, trigs, eprops, Clu, Res, Map ...
    , padBuffer ...
    , doLFP, staWin, tempSD, spatBin, whitenFlag, normalizeSWS, recomputeBL, scalef, refchan, nSD, plotPTA, minFract, vflagPTA...
    , doSPK, ilevel, clustr, shanknums, smoothphase, toplot ...
    , Overwrite, graphics, savetype...
    ] = ParseArgPairs(...
    { 'suffix', 'ripBP', 'filtmode', 'phsmode', 'OverwritePHS', 'nbins', 'iperiods' ...
    , 'stimType', 'wavRange', 'sourceType', 'durRange', 'valRange', 'uStim', 'cmp', 'simOnly', 'minEvents' ...
    , 'channel', 'periods', 'trigs', 'eprops', 'Clu', 'Res', 'Map' ...
    , 'padBuffer'...
    , 'doLFP', 'staWin', 'tempSD', 'spatBin', 'whitenFlag', 'normalizeSWS', 'recomputeBL', 'scalef', 'refchan', 'nSD', 'plotPTA', 'minFract', 'vflagPTA'...
    , 'doSPK', 'ilevel', 'clustr', 'shanknums', 'smoothphase', 'toplot'...
    , 'Overwrite', 'graphics', 'savetype'...
    }...
    , { 'eeg', [ 80 250 ], 'gaussian', 'hilbert', -3, 20, [] ...
    , { 'PULSE', 'PSINE' }, [ 400 500 ], { 'LED', 'LD' }, [ 0.05 0.07 ], [ 0 inf ], 1, 'eq', 0, 5 ...
    , [], [], [], [], [], [], [] ...
    , [ -0.01 0.01 ] ...
    , 0, [], 0.0001, 1, 0, 1e3, 0, [], [], inf, 0, 0, 0 ...
    , [], 'B', 'clu', [], 1, 'gain'...
    , -2, [ 0 0 0 0 ], 'png'...
    }...
    , varargin{ : } );

if isempty( doSPK )
    if isequal( trigMode, 'all' )
        doSPK                       = 0;
    else
        doSPK                       = 1;
    end
end
if ~ismember( lower( suffix ), { 'eeg', 'dat' } )
    error( 'erroneous specification of suffix' )
end
if ~ismember( lower( filtmode ), { 'lowpass', 'highpass', 'bandpass', 'median', 'gaussian' } )
    error( 'erroneous specification of filtmode' )
end
if ~ismember( lower( phsmode ), { 'hilbert', 'troughs', 'peaks', 'wavelet' } )
    error( 'erroneous specification of phsmode' )
end
if ~ismember( lower( toplot ), { 'rate', 'gain', 'occurrence' } )
    error( 'erroneous specification of toplot' )
end
if ~ismember( OverwritePHS, -3 : 1 )
    error( 'erroneous specification of OverwritePHS' )
end

%------------------------------------------------------------------------
% preparations
%------------------------------------------------------------------------
% filebase
[ pathname, filename, extname ]     = fileparts( filebase ); 
filename                            = [ filename extname ];
if ~exist( pathname, 'dir' )
    msg = sprintf( '%s: Missing %s directory file!!', upper( mfilename ), pathname );
    throw( MException( 'hfoAnalysisSpiking:MissingPRMfile', msg ) );
end

% files
if ~exist( [ filebase '.prm.xml' ], 'file' )
    msg = sprintf( '%s: Missing %s.prm.xml file!!', upper( mfilename ), filebase );
    throw( MException( 'hfoAnalysisSpiking:MissingPRMfile', msg ) );
end
if ~exist( [ filebase '.eeg' ], 'file' )
    msg = sprintf( '%s: Missing %s.eeg file!!', upper( mfilename ), filebase );
    throw( MException( 'hfoAnalysisSpiking:MissingEEGfile', msg ) );
end
if ~exist( [ filebase '.sst' ], 'file' )
    msg = sprintf( '%s: Missing %s.sst file!!', upper( mfilenamed ), filebase );
    throw( MException( 'hfoAnalysisSpiking:MissingSSTfile', msg ) );
end

% par file parameters
par                                 = LoadXml( filebase );
spkFs                               = par.SampleRate;
nchans                              = par.nChannels;
switch suffix
    case 'eeg'
        Fs                          = par.lfpSampleRate; 
    case 'dat'
        Fs                          = par.SampleRate; 
end

% paths
delim                               = strfind( filebase, '/dat/' );
if isempty( delim )
    fprintf( '%s: Cannot save fig and/or data\n', upper( mfilename ) )
end
if isa( graphics, 'char' ) && exist( graphics, 'dir' )
    figdir                          = graphics;
    graphics                        = ones( 1, ng );
else
    figdir                          = [ filebase( 1 : delim ) 'figs/hfo' ];
    if ~exist( fileparts( figdir ), 'dir' )
        mkdir( fileparts( fileparts( figdir ) ), 'figs' )
    end
    if ~exist( figdir, 'dir' )
        mkdir( fileparts( figdir ), 'hfo' )
    end
end
matdir                              = [ filebase( 1 : delim ) 'mat/hfo' ];
if ~exist( fileparts( matdir ), 'dir' )
    mkdir( fileparts( fileparts( matdir ) ), 'mat' )
end
if ~exist( matdir, 'dir' )
    mkdir( fileparts( matdir ), 'hfo' )
end
figname                             = sprintf( '%s/%s.hfo.spiking.%s', figdir, filename, trigMode );
savename                            = sprintf( '%s/%s.hfo.spiking.%s', matdir, filename, trigMode );
savename2                           = sprintf( '%s/%s.hfo.spiking.%s.details', matdir, filename, trigMode );
strs                                = { figname, savename, savename2 };

lg                                  = length( graphics );
graphics                            = [ graphics zeros( 1, ng - lg ) ];

% set up the fir/s
switch lower( filtmode )
    case 'bandpass'
        % bandpass filtering
        firs{ 1 }                   = makefir( freqBP, spkFs, [], 'bandpass' );
    case 'highpass'
        % highpass then lowpass filtering
        firs{ 1 }                   = makefir( [ freqBP( 1 ) NaN ], spkFs, [], 'hipass' );
        firs{ 2 }                   = makefir( [ NaN freqBP( 2 ) ], spkFs, [], 'lopass' );
    case 'median'
        % high-pass by subtracting the output of a median (nonlinear) filter:
        firs{ 1 }                   = complex( 0, ceil( spkFs/freqBP( 1 ) ) );
        firs{ 2 }                   = makefir( [ NaN freqBP( 2 ) ], spkFs, [], 'lopass' );
    case 'lowpass'
        % high-pass by subtracting the output of a low-pass fir:
        firs{ 1 }                   = complex( 0, makefir( [ NaN freqBP( 1 ) ], spkFs, [], 'lopass' ) );
        firs{ 2 }                   = makefir( [ NaN freqBP( 2 ) ], spkFs, [], 'lopass' );
    case 'gaussian' 
        % high-pass by subtracting the output of a low-pass fir:
        firs{ 1 }                   = complex( 0, makegausslpfir( freqBP( 1 ), spkFs ) );
        firs{ 2 }                   = makegausslpfir( freqBP( 2 ), spkFs );
end

% check if need to recompute at all
if Overwrite < 0 && exist( savename, 'file' ) && trigMode( 1 ) ~= 'i'       % cannot determine file name yet if induced..
    fprintf( '%s: loading %s...\n', upper( mfilename ), savename )
    load( savename, '-mat', 's' );
    fprintf( '%s: loading %s...\n', upper( mfilename ), savename2 )
    load( savename2, '-mat', 'stats' );
    return
end

%------------------------------------------------------------------------
% channel and periods for phase determination
%------------------------------------------------------------------------
verb( sprintf( '%s: Working on %s: determining events and periods...', upper( mfilename ), filename ), verbose )

% select channel
if isempty( channel )
    switch trigMode
%         case { 'induced', 'external' }
%             if isempty( channel )
%                 msg                 = sprintf( '%s: Must specify a channel for externally-defined events!', upper( mfilename ) );
%                 throw( MException( 'hfoAnalysisSpiking:MissingExternalChannel', msg ) );
%             end
        case { 'induced','spontaneous', 'all' }
            % select channel based on highest SD over all selected channels
            spsfile                 = [ filebase '.sps' ];
            if ~exist( spsfile, 'file' )
                msg                 = sprintf( '%s: Must have an *sps file for automatic channel selection (spontaneous events)!', upper( mfilename ) );
                throw( MException( 'hfoAnalysisSpiking:MissingSpsFile', msg ) );
            end
            sps                     = load( spsfile, '-mat', 'stats' );
            % new method - max SD, but only consider 110-210 Hz range!
            vidx                    = find( inrange( sps.stats( :, 5 ), [ 110 210 ] ) );
            [ ~, maxidx ]           = max( sps.stats( vidx, 7 ) );
            channel                 = sps.stats( vidx( maxidx ), 2 );
    end
end
if isempty( channel )
    fprintf( '%s: no channel matching the requirements!\n', upper( mfilename ) );
    return
end
channel                             = channel( 1 );                         % neuronal channel by which to determine the phase
fprintf( '%s: channel %d selected!\n', upper( mfilename ), channel )

% select events
switch trigMode
    case { 'induced', 'all' }
        if isempty( periods )
            % get the shank of the chosen channel 
            emap                    = get_emap( filebase );
            shank                   = emap( ismember( emap( :, 3 ), channel ), 1 );
            % get the valid target shanks
            [ allstimchans, alltargets, ~, allsources, allwavelengths ] = get_stimchans( filebase );
            idx                     = inrange( allwavelengths, wavRange ) & ismember( allsources, sourceType );
            allstimchans            = allstimchans( idx );
            alltargets              = alltargets( idx );
            [ ualltargets, uidx ]   = unique( alltargets, 'last' );
            if ~isequal( ualltargets, alltargets )
                alltargets          = ualltargets;
                allstimchans        = allstimchans( uidx );
            end
            if ( all( alltargets > 100 ) || isequal( alltargets, 0 ) )
                if length( allstimchans ) == 1
                    trigchan        = allstimchans;
                else
                    error( 'Targets above 100 not supported in the present framework' )
                end
            else
                trigchan            = allstimchans( alltargets == shank ); % stimulation channel
            end
            if isempty( trigchan )
                fprintf( '%s: no trigger channel matches the requirements!\n', upper( mfilename ) );
                return
            end
            
            % get the external events
            params                  = { 'types', stimType, 'durs', durRange, 'vals', valRange };
            [ ~, tims, durs, ignvals, ~ ]   = get_triggers( filebase, trigchan, uStim, cmp, simOnly, params );
            periods                 = [ tims tims + durs * spkFs - 1 ] / spkFs; % [sec]
            if isempty( periods ) || size( periods, 1 ) < minEvents
                fprintf( '%s: no periods matching the requirements!\n', upper( mfilename ) );
                return
            end
            trigs                   = periods( :, 1 ) + 1 / mean( freqBP ) / 2;
            nn                      = length( trigs );
            eprops                  = [ ignvals NaN * ones( nn, 1 ) durs NaN * ones( nn, 1 ) ];% SD freq duration power
            if strcmp( trigMode, 'all' )
                rips                = plotHFOs( filebase, channel );
                periods2            = rips.edges / rips.Fs;
                trigs2              = rips.trigs / rips.Fs;
                if size( uniteranges( periods, periods2 ), 1 ) ~= ( size( periods, 1 ) + size( periods2, 1 ) )
                    error( '%s: internal error - overlapping spontaneous and induced events!!', upper( mfilename ) )
                end
                periods             = [ periods; periods2 ];
                trigs               = [ trigs; trigs2 ];
                [ trigs, sidx ]     = sort( trigs );
                periods             = periods( sidx, : );
            end
        end
    case 'spontaneous'
        if isempty( periods )
            spwfname                = [ filebase '.spw.', num3str( channel ) ];
            if exist( spwfname, 'file' )
                load( spwfname, '-mat', 'rips' )
            else
                rips                = plotHFOs( filebase, channel );
            end
            periods                 = rips.edges / rips.Fs;
        end
        trigs                       = rips.trigs / rips.Fs;
        eprops                      = [ rips.sd rips.f diff( rips.edges, 1, 2 ) + 1 rips.pow ]; % SD freq duration power
        ridx                        = diff( periods, 1, 2 ) > maxDurSEC;
        trigs( ridx )               = [];
        periods( ridx, : )          = [];
        eprops( ridx, : )           = [];
        
    case 'external'
        if isempty( periods ) || isempty( trigs )
            msg                     = sprintf( '%s: Must specify periods and trigs for externally-defined events!', upper( mfilename ) );
            throw( MException( 'hfoAnalysisSpiking:MissingExternalEvents', msg ) );
        end
        if isempty( eprops ) || size( eprops, 1 ) ~= size( trigs, 1 ) || size( eprops, 2 ) ~= 4
            nn                      = size( periods, 1 );
            nans                    = NaN * ones( nn, 1 );
            durs                    = diff( periods, 1, 2 );
            eprops                  = [ nans nans durs nans ];% SD freq duration power
        end
end
% clean periods
if ~isempty( iperiods )
    [ ~, idx1, ~ ]                  = intersectranges( periods, iperiods );
    periods                         = periods( idx1, : );
    trigs                           = trigs( idx1, : );
end
% should also make sure that ranges do not overlap etc..
if ~isequal( size( sortranges( periods ) ), size( periods ) )
    [ speriods, idx1 ]              = sortranges( periods );
    fprintf( '%s: %d overlapping periods in %s!! merging\n'...
        , upper( mfilename ), size( periods, 1 ) - length( idx1 ), filename )
    periods                         = speriods;
    trigs                           = trigs( idx1, : );
    eprops                          = eprops( idx1, : );
    eprops( :, 3 )                  = diff( speriods, 1, 2 ) + 1 / spkFs;
end

if isempty( periods ) || size( periods, 1 ) < minEvents
    fprintf( '%s: no periods matching the requirements!\n', upper( mfilename ) );
    return
end

nevents                             = size( periods, 1 );
if nevents ~= numel( trigs )
    msg                             = sprintf( '%s: Periods/trigs size mismatch!', upper( mfilename ) );
    throw( MException( 'hfoAnalysisSpiking:EventSpecificationMismatch', msg ) );
end
if ~exist( 'eprops', 'var' )
    eprops                          = []; 
end

% determine units (needed for the saved filename)
shankclu                            = determine_units( filebase, shanknums, ilevel );
shanknums                           = unique( shankclu( :, 1 ) );

% spiking shank (string for induced *mat name/fig titles)
if ~isempty( shanknums )
    shstr                           = ( [ repmat( 'S', size( shanknums( : ) ) ) num2str( shanknums( : ) ) ]' );
    shstr                           = [ replacetok( shstr( : ).', '0', ' ' ) '_' ];
else
    shstr                           = '';
end

% target file names
if trigMode( 1 ) == 'i'
    savebase                        = sprintf( '%s.hfo_spiking_%s', filename, trigMode );
    % trigger channels
    tstr                            = ( [ repmat( 'T', size( trigchan( : ) ) ) num2str( trigchan( : ) ) ]' );
    tstr                            = [ tstr( : ).' '_' ];
    chstr                           = sprintf( 'C%d_', channel );
    % trigger properties
    if isa( stimType, 'cell' )
        stmtype                     = stimType{ 1 };
        for i                       = 2 : length( stimType )
            stmtype                 = sprintf( '%s_%s', stmtype, stimType{ i } );
        end
    else
        stmtype                     = stimType;
    end
    if any( isinf( valRange ) )
        valRangeStr                 = minmax( eprops( :, 1 ) );
    else
        valRangeStr                 = valRange;
    end
    stmstr                          = sprintf( '%s_w%dw%d_v%0.3gv%0.3g_d%dd%d', stmtype...
        , round( wavRange( 1 ) ), round( wavRange( 2 ) )...
        , remrnd( valRangeStr( 1 ), 0.0001, 'round' ), remrnd( valRangeStr( 2 ), 0.0001, 'round' )...
        , ceil( durRange( 1 ) * 1000 ), ceil( durRange( 2 ) * 1000 ) ); % val in 0.1mV, time in 1 ms
    stmstr                          = replacetok( stmstr, '#', '.' );
    dstr                            = sprintf( '%s%s%s%s', tstr, chstr, shstr, stmstr );
    corename                        = sprintf( '%s.%s', savebase, dstr );
    figname                         = [ figdir '/' corename ];
    savename                        = [ matdir '/' corename ];
    savename2                       = [ matdir '/' corename '.details' ];
    strs                            = { figname, savename, savename2 };
end

if Overwrite < 0 && exist( savename, 'file' ) % cannot determine file name yet if induced..
    fprintf( '%s: loading %s...\n', upper( mfilename ), savename )
    load( savename, '-mat', 's' );
    fprintf( '%s: loading %s...\n', upper( mfilename ), savename2 )
    load( savename2, '-mat', 'stats' );
    return
end

%------------------------------------------------------------------------
% get the phases for those periods
%------------------------------------------------------------------------
% detailed design of this section:
% -check if there's a saved file, load it
% -check if the desired data segments are there
%   if not: -determine which segments are missing and compute phases for those
%           -combine with the loaded segments (if any) 
%           -save the combined data (if desired)
% -extract the desired segments for further internal processing
tims                                = [];
uphs                                = [];

sourcefile                          = [ filebase '.' suffix ];
bpstr                               = replacetok( sprintf( '%3d_%3d', freqBP( 1 ), freqBP( 2 ) ), '0', ' ' );
phsfname                            = sprintf( '%s.phs.%s.%s', filebase, bpstr, num3str( channel ) );
periodsUps                          = resampleranges( periods * Fs, spkFs, Fs );% []->spkFs
% to go back: 
periodsMissing                      = periodsUps;                           % [ samples @ spkFs ]
trigsMissing                        = trigs;                                % []

% load
if OverwritePHS < 0 && exist( phsfname, 'file' )
    verb( sprintf( '%s: loading %s...', upper( mfilename ), phsfname ), verbose )
    PHS                             = load( phsfname, '-mat' );
    uphs                            = PHS.uphs;
    tims                            = PHS.tims;
    % check that overlapping exactly
    mat                             = parse( tims );
    if isequal( mat, periodsUps )
        periodsMissing              = [];
    else
        [ tmp, sidx ]               = setdiffranges( periodsMissing, mat );
        if isequal( periodsMissing( sidx, : ), tmp )
            periodsMissing          = periodsMissing( sidx, : );
            trigsMissing            = trigsMissing( sidx );
        end
    end
end

% recompute
periodsMissing                      = resampleranges( periodsMissing, Fs, spkFs ) / Fs; % spksFs->[sec]
if isempty( periodsMissing )
    timsNew                         = [];
    uphsNew                         = [];
else
    verb( sprintf( '%s: Computing %s phases (%d-%d) from %s file: %0.3g sec of data...\n'...
        , upper( mfilename ), phsmode, freqBP( 1 ), freqBP( 2 ), suffix, sum( diff( periodsMissing, 1, 2 ) ) ), -verbose )
    [ ~, timsNew, uphsNew ]         = eegPhase( sourcefile, channel, periodsMissing...
        , firs, [], Fs, spkFs, nchans, phsmode, 1 / Fs * 2, trigsMissing );
    verb( 'Done', verbose )
    if sum( diff( resampleranges( periodsMissing * Fs, spkFs, Fs ), 1, 2 ) + 1 ) ~= length( timsNew )
        error( '%s: internal error - check eegPhase!!', upper( mfilename ) )
    end
end

% now combine all and sort
if exist( 'PHS', 'var' )
    tims                            = [ tims; timsNew ];
    uphs                            = [ uphs; uphsNew ];
    if ~issorted( tims )
        [ ~, sidx ]                 = sort( tims );
        tims                        = tims( sidx );
        uphs                        = uphs( sidx );
    end
else
    tims                            = timsNew;
    uphs                            = uphsNew;
end

% save the phases and times:
if OverwritePHS == 1 || ( OverwritePHS ~= -1 && ~exist( phsfname, 'file' ) ) || ( OverwritePHS == -3 && ~isempty( timsNew ) )
    Fs0                             = Fs;
    Fs                              = spkFs;
    generator                       = { computer, datestr( now, 'ddmmmyy' ) };
    save( phsfname, 'filebase', 'generator', 'Fs', 'firs', 'phsmode', 'suffix', 'uphs', 'tims' );
    verb( sprintf( '%s: !!!! Saved %s... !!!!! ', upper( mfilename ), phsfname ), verbose )
    Fs                              = Fs0;
end

% determine exclusion epochs
vals                                = LoadStims( filebase );                % @spkFs
if isempty( vals )
    vals                            = LoadVals( filebase );
end
pad                                 = [ floor( padBuffer( 1 ) *spkFs ) ceil( padBuffer( 2 ) * spkFs ) ];
xperiods                            = [ vals( :, 1 ) + pad( 1 ) vals( :, 2 ) + pad( 2 ) ];
xperiods                            = sortranges( xperiods );
if strcmp( trigMode, 'spontaneous' )
    
    % first, determine the periods to remove - anything that overlaps with
    % the xperiods should be removed
    % then, determine which samples to keep - anything that is inside the
    % remaining periods
    
    % to do all this - upsample the periods first:
    periodsUps                      = resampleranges( periods * Fs, spkFs, Fs );
    
    % then determine overlap with xperiods:
    rmvperiods                      = isoverlap( periodsUps, xperiods );
    periodsUps( rmvperiods, : )     = [];
    periods                         = resampleranges( periodsUps, Fs, spkFs ) / Fs;
    
    % now determine which samples of tims to keep:
    ridx                            = true( size( tims ) );
    ridx( inranges( tims, periodsUps ) ) = 0;
    nevents0                        = nevents;
    nevents                         = size( periods, 1 );
    fprintf( '%s: excluding %d/%d sec: keeping %d/%d events\n'...
        , upper( mfilename ), round( sum( ridx ) / spkFs )...
        , round( length( tims ) / spkFs ), nevents, nevents0 )
    tims( ridx )                    = [];
    uphs( ridx )                    = [];
    kidx                            = inranges( trigs, periods );
    trigs                           = trigs( kidx );
    eprops                          = eprops( kidx, : ) ;
    
    % check again at end:
    if ~isequal( resampleranges( periods * Fs, spkFs, Fs ), parse( tims ) )
        error( 'internal resmapling error..' )
    end
end

% extract the relevant segments only
fprintf( '%s: extracting params and computing cycle statistics..\n', upper( mfilename ) )
[ mat, midx ]                       = parse( tims );
[ ~, sidx ]                         = intersectranges( mat, periodsUps );
if length( sidx ) ~= size( mat, 1 )
    midx                            = midx( sidx, : );
    idx                             = zeros( sum( diff( midx, 1, 2 ) + 1 ), 1 );
    i0                              = 0;
    for i                           = 1 : size( midx, 1 )
        ni                          = diff( midx( i, : ) ) + 1;
        tidx                        = ( midx( i, 1 ) : midx( i, 2 ) );
        idx( i0 + ( 1 : ni ) )      = tidx;
        i0                          = i0 + ni;
    end
    tims                            = tims( idx );
    uphs                            = uphs( idx );
end

% could still be edge effect
mat                                 = parse( tims );
dmat                                = setdiffranges( mat, periodsUps );
ridx                                = inranges( tims, dmat );
tims( ridx )                        = [];
uphs( ridx )                        = [];
ok                                  = isequal( setdiffranges( parse( tims ), dmat ), periodsUps );
fprintf( '%s: %d! After trimming: %d phase samples (%d desired) in %d epochs\n'...
    , upper( mfilename ), ok, sum( diff( setdiffranges( parse( tims ), dmat ), 1, 2 ) + 1 )...
    , sum( diff( periodsUps, 1, 2 ) + 1 ), size( periodsUps, 1 ) )

% convert tims + phs -> cycs, and tims -> evnt
[ mat, midx ]                       = parse( tims );
cycs                                = zeros( size( uphs ), 'single' );
evnt                                = zeros( size( uphs ), 'single' );
for i                               = 1 : size( mat, 1 )
    idx                             = midx( i, 1 ) : midx( i, 2 );
    cycnum                          = floor( uphs( idx ) / ( 2 * pi ) );
    cycnum                          = [ -flipud( monotonic( flipud( -cycnum( cycnum < 0 ) ) ) ); cycnum( cycnum == 0 ); monotonic( cycnum( cycnum > 0 ) ) ];
    cycs( idx )                     = cycnum;
    evnt( idx )                     = i;
end
% to convert uphs -> phs: phs = mod( uphs, 2 * pi );

% some elementary stats on the cycs paramters:
ucyc                                = unique( cycs );
cedges                              = unique( [ ucyc - 0.5; ucyc + 0.5 ] ); 
ccount                              = histc( cycs, cedges ); 
ccount( end )                       = [];
if diff( minmax( cycs ) ) + 1 == length( ccount )
    aa                              = ccount;
else
    aa                              = uhist( cycs );
end
nucycs                              = length( ucyc );
nval                                = zeros( nucycs, 1 );
mdur                                = nval;
sdur                                = nval;
for i                               = 1 : nucycs
    fidx                            = find( cycs == ucyc( i ) );
    mat                             = parse( fidx );
    nval( i )                       = size( mat, 1 );
    dur                             = diff( mat, 1, 2 );
    mdur( i )                       = mean( dur ) / spkFs;
    sdur( i )                       = calc_sem( dur ) / spkFs;
end

if ~doSPK
    return
end

%------------------------------------------------------------------------\
% get the spikes
%------------------------------------------------------------------------\
verb( sprintf( '%s: Loading unit properties...', upper( mfilename ) ), verbose )
nclu                                = size( shankclu, 1 );
if nclu == 0
    fprintf( 1, '%s: NO VALID units at %s for %s\n', upper( mfilename ), ilevel, filebase )
    return
end
verb( sprintf( '%s: %d %s-class units ', upper( mfilename ), size( shankclu, 1 ), ilevel ), verbose )
if isempty( Clu ) || isempty( Res ) || isempty( Map )
    [ Clu, Res, Map ]               = get_spikes( filebase, shankclu, clustr );
end
if isempty( Clu )
    verb( sprintf( 'NO SPIKES in %s, shanks %s\n', filebase, num2str( shanknums ) ), verbose )
    return
end
verb( sprintf( '%s: %d spikes', upper( mfilename ), length( Clu ) ), verbose )

%------------------------------------------------------------------------\
% allocate a phase to each spike
%------------------------------------------------------------------------\
verb( sprintf( '%s: Computing phase histograms...', upper( mfilename ) ), verbose )

% get the total samples
info                                = dir( sourcefile );
lens                                = info.bytes / ( par.nBits / 8 ) / par.nChannels;
totsamps                            = round( sum( lens ) * Fs / par.lfpSampleRate ); % total samples @ eegFs

% set up some parameters (for the histograms)
mmcycs                              = bounds( cycs, 1 - cycSupport  );
mmcycs                              = [ max( [ mmcycs( 1 ) cycRange( 1 ) ] ) min( [ mmcycs( 2 ) cycRange( 2 ) ] ) ];
ucyc0                               = unique( cycs );
tmp                                 = ucyc0( inrange( ucyc0, mmcycs ) );
if sum( diff( tmp ) > 1 )
    tmp( 1 : find( diff( tmp ) > 1 ) ) = [];
    mmcycs                          = intersectranges( minmax( tmp ), mmcycs );
end
minPhase                            = mmcycs( 1 ) * 2 * pi;
maxPhase                            = ( mmcycs( 2 ) + 1 ) * 2  * pi;
binsize                             = 2 * pi / nbins;
phsEdges                            = ( minPhase - binsize / 2 : binsize : maxPhase - binsize / 2 )'; % ignore the very last bin to have integer multiples
phsBins                             = ( phsEdges( 1 : end - 1 ) + phsEdges( 2 : end ) ) / 2;
phsBins                             = round( phsBins * nbins / ( 2 * pi ) ) / nbins * ( 2 * pi );
totbins                             = length( phsBins );
totcycs                             = totbins / nbins;
mcycdur                             = sum( nval .* mdur ) / sum( nval );    % [s] - averaged over all cycles
phsBinSize                          = 1 / nbins * mcycdur;

% get the baseline rates
xsamples                            = sum( diff( xperiods, 1, 2 ) + 1 );    % @spkFs
ridx                                = inranges( Res, xperiods );
CluBase                             = Clu;
CluBase( ridx )                     = [];
[ count, clunum ]                   = uhist( CluBase );                     % this is not precise: should exclude here (1) stim (2) hfo events
basesamples                         = totsamps * spkFs / Fs - xsamples;     % spkFs
ratesBL                             = zeros( 1, nclu );
ratesBL( ismember( Map( :, 1 ), clunum ) ) = count / basesamples * spkFs;

% also get the mean rates
insamples                           = sum( diff( periodsUps, 1, 2 ) + 1 );  % @spkFs
kidx                                = inranges( Res, periodsUps );
CluIn                               = Clu( kidx );
[ count, clunum ]                   = uhist( CluIn );
ratesMean                           = zeros( 1, nclu );
ratesMean( ismember( Map( :, 1 ), clunum ) ) = count / insamples * spkFs;

% get the multi-cycle, multi-event spike phases (rasters)
ptims                               = parse( tims );
idx0                                = inranges( Res, ptims );               % spikes in events
clu                                 = Clu( idx0 );
res                                 = Res( idx0 );
[ ~, idx ]                          = ismember( res, tims );                % phases/events for those
spkPhs                              = uphs( idx );
spkEvt                              = double( evnt( idx ) );
[ spkPhs, sidx ]                    = sort( spkPhs );
spkEvt                              = spkEvt( sidx );
clu                                 = clu( sidx );
res                                 = res( sidx );
binMatrix                           = phsEdges( [ ( 1 : length( phsEdges ) - 1 )' ( 2 : length( phsEdges ) )' ] );
[ kidx, spkBins ]                   = inranges( spkPhs, binMatrix );        % binned phases
ignPhs                              = sort( spkPhs( setdiff( 1 : length( spkEvt ), kidx ) ) )';
if ~isempty( ignPhs )
    verb( sprintf( '%s: %d/%d spikes out of phase range'...
        , upper( mfilename ), length( ignPhs ), length( spkEvt ) ), verbose )
end
spkEvt                              = spkEvt( kidx );                       % event identity
spkTim                              = totbins * spkEvt + spkBins;           % event-resolved phase bin
spkPhs                              = spkPhs( kidx );                       % unwrapped phase [rad]
clu                                 = clu( kidx );                          % clu
res                                 = res( kidx );                          % time [samples]@spkFs
[ spkTim, sidx ]                    = sort( spkTim );                       % sort again by time
spkBins                             = spkBins( sidx );
spkEvt                              = spkEvt( sidx );
spkPhs                              = spkPhs( sidx );
clu                                 = clu( sidx );
res                                 = res( sidx );
s                                   = struct( 'filebase', filebase, 'Fs' ...
    , spkFs, 'shankclu', shankclu, 'map', Map...
    , 'clu', clu, 'res', res...
    , 'phsbin', spkTim, 'phs', spkPhs, 'bins', phsBins );

% make multi-cycle histograms
phsHists                            = zeros( totbins, nclu );
for j                               = 1 : size( shankclu, 1 )
    uidx                            = Map( ismember( Map( :, 2 : 3 ), shankclu( j, 1 : 2 ), 'rows' ), 1 );
    phsbin                          = rem( s.phsbin( s.clu == uidx ) - 1, totbins ) + 1;
    if isempty( phsbin )
        continue
    end
    h                               = histc( phsbin, ( 1 : totbins + 1 ) - 0.5 ); 
    phsHists( :, j )                = h( 1 : totbins );
end

% derive wrapped histograms
phshists                            = permute( sum( reshape( phsHists, [ nbins totcycs nclu ] ), 2 ), [ 1 3 2 ] );
phsbins                             = circ_mean( reshape( phsBins, [ nbins totcycs ] )' )';
phsbins( phsbins >  2 * pi - binsize / 2 ) = phsbins( phsbins >  2 * pi - binsize / 2 ) - 2 * pi;
phsbins                             = clipmat( phsbins, [ 0 2 * pi ] );

% compute preferred phases 
[ mPhase, rPhase, sPhase ]          = circ_mean( phsBins( : ), phsHists );
sPhase                              = sPhase ./ sqrt( sum( phsHists ) );
pval                                = ray_test( phsBins( : ), phsHists );


%------------------------------------------------------------------------\
% adjust, summarize and save
%------------------------------------------------------------------------\
phsbins                             = round( phsbins * nbins / ( 2 * pi ) ) / nbins * ( 2 * pi );
phsBins                             = round( phsBins / ( 2 * pi ) * nbins ) * ( 2 * pi ) / nbins;

% expand histograms to the range cycRange to enable accumulation:
allbins                             = cycRange( 1 ) * nbins : 1 : ( cycRange( 2 ) + 1 ) * nbins;
usebins                             = round( phsBins / ( 2 * pi / nbins ) );
binidx                              = ismember( allbins, usebins );
amat                                = NaN * ones( nclu, length( allbins ) );
amat( :, binidx )                   = phsHists';

% direct way to evaluate the time/bin, accounting for cycle assymetry etc:
bcenters                            = allbins * ( 2 * pi / nbins );
alledges                            = [ bcenters - pi / nbins bcenters( end ) + pi / nbins ];
hh                                  = histc( uphs, alledges ); 
hh( end )                           = [];                                   % number of samples / phase bin
avec                                = ceil( hh / ( spkFs * phsBinSize ) );

% compute gain
gain                                = bsxfun( @rdivide, bsxfun( @rdivide, amat', avec ) / phsBinSize, ratesBL );

% make a raster for each unit (expanded to fullwidth):
rbins                               = find( binidx ); % mapping from spkBins/usebins -> allbins
nfullbins                           = length( allbins );
rast                                = cell( nclu, 1 );
for i                               = 1 : nclu
    cidx                            = clu == Map( i, 1 );
    eidx                            = spkEvt( cidx );
    bidx                            = spkBins( cidx );
    rast{ i }                       = sparse( rbins( bidx ), eidx, 1, nfullbins, nevents );
end

% make a multi-unit raster for each event:
r                                   = cell( nevents, 1 );
nspks                               = zeros( nevents, 1 );
for j                               = 1 : nevents
    r{ j }                          = sparse( nfullbins, nclu );
    for i                           = 1 : nclu
        r{ j }( :, i )              = rast{ i }( :, j );
    end
    nspks( j )                      = full( sum( r{ j }( : ) ) );
end

% summarize in a struct:
stats                               = struct( 'filename', filename );
stats.phsBins                       = allbins * ( 2 * pi / nbins );         % all possible bins
stats.phsbins                       = phsbins';                             % 1 cycles
stats.periods                       = ptims;                                % actual periods @ spkFs
stats.trigs                         = round( trigs * spkFs );
stats.eprops                        = eprops;
stats.trigfname                     = repmat( { filename }, [ nevents 1 ] );
stats.trigchan                      = ones( nevents, 1 ) * channel;
stats.Fs                            = spkFs * ones( nclu, 1 );
stats.channel                       = channel * ones( nclu, 1 );
stats.filename                      = repmat( { filename }, [ nclu 1 ] );
stats.phsBinSize                    = phsBinSize * ones( nclu, 1 );
stats.denom                         = ones( nclu, 1 ) * avec';
stats.shankclu                      = shankclu;
stats.nevents                       = nevents * ones( nclu, 1 );
stats.phsHists                      = amat;
stats.gain                          = gain';
stats.phshists                      = phshists';                            % 20 bins
stats.ratesMean                     = ratesMean';
stats.ratesBL                       = ratesBL';
stats.rast                          = rast;
stats.mPhase                        = mPhase';
stats.rPhase                        = rPhase';
stats.pval                          = pval';

% add the depth relative to layer.. (relevant actually only for CA1..)
sst                                 = spikes_stats_depth( filebase, 'graphics', 0, 'Overwrite', -2 );%, 'flipLFP', toflipRelative );
kidx                                = ismember( sst.shankclu( :, 1 : 2 ), stats.shankclu( :, 1 : 2 ), 'rows' );
stats.depth                         = sst.depth( kidx );

fprintf( 1, '%s\nMean PYR phase: %0.3g (%d units)\nMean INT phase: %0.3g (%d units)\n%s\n'...
    , repmat( '*', [ 1 80 ] )...
    , circ_mean( mPhase( shankclu( :, 3 ) == 1 )' ) * 180 / pi, sum( shankclu( :, 3 ) == 1 )...
    , circ_mean( mPhase( shankclu( :, 3 ) == 0 )' ) * 180 / pi, sum( shankclu( :, 3 ) == 0 )...
    , repmat( '*', [ 1 80 ] ) );

if ( Overwrite == 1 || ~isempty( savename ) && ~exist( savename, 'file' ) && Overwrite ~= -1 )
    verb( sprintf( '%s: saving %s', upper( mfilename ), savename ), verbose )
    save( savename, 's', '-v6' );
    verb( sprintf( '%s: saving %s', upper( mfilename ), savename2 ), verbose )
    save( savename2, 'stats', '-v6' );
end


%------------------------------------------------------------------------\
% plot the results:
%------------------------------------------------------------------------\

% REPLACE THIS PART BY AN EXTERNAL ROUTINE!!

pidx                                = shankclu( :, 3 );
if smoothphase % smooth the histograms
    win                             = makegaussfir( 1, 1 );
end

cyclims                             = minmax( ucyc0( nval > 0.2 * size( periods, 1 ) ) ); % focus on the cycles with >=20% occurrence
cyclims                             = max( abs( cyclims ) ) * [ -1 1 ] + [ 0 1 ]; % make symmetric around cycle0
cyclims                             = double( cyclims );
cyclims                             = intersectranges( mmcycs, cyclims );
if ~strcmp( trigMode, 'spontaneous' )
    cyclims                         = intersectranges( [ -1 max( cyclims ) ], cyclims );
end
    
% get the lfp for those 
if doLFP
    if isempty( staWin )
        verb( sprintf( '%s: Loading LFP data...', upper( mfilename ) ), verbose )
        staWin                      = remrnd( 1 / mean( freqBP ) * max( abs( cyclims ) ) * 1.5, 1 / Fs, 'ceil' ) * [ -1 1 ];
    end
     % approx; should warp the LFP.. as is, this is VERY confusing and best avoided.
    specMode                        = 'eeg';
    specChan                        = [];
    
    [ ~, gchans ]                   = get_egroup( par, channel );
    if find( gchans == channel ) == 1
        csdchans                    = gchans( [ 1 1 2 ] );
    elseif find( gchans == channel ) == length( gchans )
        csdchans                    = gchans( length( gchans ) + [ -1 -1 0 ] );
    else
        csdchans                    = gchans( find( gchans == channel ) + ( -1 : 1 ) );
    end
    [ avgcsd, ~, pttim, xcsdhat ]   = pt_avg( filebase, csdchans, trigs...
        , 'stawinSEC', staWin, 'tempSD', tempSD, 'spatBin', spatBin, 'suffix', specMode...
        , 'graphics', plotPTA, 'scalefactor', scalef, 'specChans', specChan...
        , 'whitenFlag', whitenFlag, 'normalizeSWS', normalizeSWS...
        , 'savetype', savetype, 'refchan', refchan...
        , 'recomputeBL', recomputeBL, 'nSD', nSD, 'minFract', minFract, 'vflag', vflagPTA );
    avgcsd                          = avgcsd( :, 2 );
    xcsdhat                         = xcsdhat( 2, :, : );
    
    xlfphat                         = permute( xcsdhat, [ 2 3 1 ] );
    xlfp                            = fft_upsample( avgcsd, spkFs / Fs );
    ptphs                           = pttim / 1000 / phsBinSize / nbins + 0.5;
    xphs                            = fft_upsample( ptphs, spkFs / Fs );
    mat                             = parse( find( diff( xphs ) > 0 ) );
    [ ~, midx ]                     = max( diff( mat, 1, 2 ) );
    kidx                            = mat( midx, 1 ) : mat( midx, 2 );
    xphs                            = xphs( kidx );
    xlfp                            = xlfp( kidx );
end

if sum( graphics ) > 0
    titstr                          = sprintf( '%s, channel %d (%s); %d events, %0.3g sec'...
        , replacetok( filename, '\_', '_' ), channel, shstr( 1 : end - 1 )...
        , nevents, sum( diff( periods, 1, 2 ) ) );
end

% plot the phase/cycle statistics:
if graphics( 1 )
    
    phsedges                        = 0 : binsize : 2 * pi;
    h                               = NaN * ones( nbins, nucycs );
    for i                           = 1 : nucycs
        fidx                        = cycs == ucyc( i );
        h0                          = hist( mod( uphs( fidx ), 2 * pi ), phsedges );
        h0( end )                   = [];
        h( :, i )                   = h0;
    end
    
    fig( 1 )                        = figure;
    subplot( 2, 2, 1 )
    bh                              = bar( ucyc, aa / spkFs );
    xlim( ucyc( [ 1 end ] ) ),
    ylabel( 'Time in cycle [s]' ),
    title( sprintf( '%d events (%0.3g h)', size( periods, 1 ), totsamps / Fs / 3600 ) );
    set( bh, 'EdgeColor', [ 0 0 1 ], 'FaceColor', [ 0 0 1 ] );
    xlim( cyclims * 2 )
    
    subplot( 2, 2, 3 )
    bh                              = bar( ucyc, nval ); 
    xlim( ucyc( [ 1 end ] ) )
    ylabel( 'Number of cycles' )
    set( bh, 'EdgeColor', [ 0 0 1 ], 'FaceColor', [ 0 0 1 ] );
    xlabel( 'Cycle number' )
    alines( size( periods, 1 ) * ( 0 : 0.2 : 1 ), 'y', 'color', [ 0 0 0], 'linestyle', '--' );
    xlim( cyclims * 2 )
    
    subplot( 2, 2, 2 )
    barwerror( ucyc, mdur * 1000, sdur * 1000 );
    xlim( cyclims + [ -0.5 0.5 ] )
    xlabel( 'Cycle number' )
    ylabel( 'Cycle duration [ms]' )
    if doLFP
        line( xphs, scaleto( xlfp, ylim ), 'color', blackColor, 'linewidth', 2 )
    end
    
    subplot( 2, 2, 4 )
    phsDist                         = bsxfun( @rdivide, h, sum( h ) ) * nbins;
    bidx                            = inrange( ucyc, cyclims );
    imagesc( ucyc( bidx ), phsbins( 2 : end ), phsDist( 2 : end, bidx ) )
    axis xy
    ylabel( 'Phase [rad]' )
    ch                              = colorbar( 'horiz' );
    set( ch, 'tickdir', 'out', 'box', 'off' )
    
    for i = 1 : 4
        subplot( 2, 2, i )
        set( gca, 'tickdir', 'out', 'box', 'off' )
    end
    colormap( myjet )
    textf( 0.5, 0.975, titstr );
    
end

if graphics( 2 )

    fig                             = hfoAnalysisSpikingPlot1( stats, [] ...
        , 'minEvents', 10, 'trigMode', trigMode...
        , 'exttitle', { '' } );
    textf( 0.5, 0.975, titstr );
    
end

if graphics( 3 )
    
    [ ah, fig( 3 ) ]                = tilefig( ceil( sqrt( nclu ) ), ceil( nclu / ceil( sqrt( nclu ) ) ), 1, 0.9 );
    xticks                          = phsBins / ( 2 * pi );
    xidx                            = inrange( xticks, cyclims );
    for i                           = 1 : nclu
        subplot( ah( i ) )
        plot_raster( rast{ i }, phsBins/ (2*pi), [], [], 1 ); 
        peth                        = full( sum( rast{ i }, 2 ) ./ stats.denom( i, : )' ) / phsBinSize;
        if smoothphase
            peth                    = firfilt( peth, win );
        end
        if isnan( pidx( i ) )
            acolor                  = [ 1 0 0.7 ];
        else
            acolor                  = colors( pidx( i ) + 1, : );
        end
        line( xticks( xidx ), scaleto( peth( xidx ), ylim ), 'color', acolor )
        xlim( cyclims )
        alines( cyclims( 1 ) : cyclims( 2 ), 'x', 'color', blackColor, 'linestyle', sepStyle );
        axis off
        title( sprintf( '%d.%d', shankclu( i, 1 ), shankclu( i, 2 ) ) )
    end
    textf( 0.5, 0.975, titstr );
    for i                           = ( nclu + 1 ) : length( ah )
        subplot( ah( i ) )
        axis off
    end
    
end

if graphics( 4 ) 
    
    % sort by population synchrony:
    phs                             = stats.phsBins( : ) * ones( 1, nclu );
    pyridx                          = stats.shankclu( :, 3 ) == 1;
    phspyr                          = stats.phsBins( : ) * ones( 1, sum( pyridx ) );
    intidx                          = stats.shankclu( :, 3 ) == 0;
    phsint                          = stats.phsBins( : ) * ones( 1, sum( intidx ) );
    
    mphs                            = NaN * ones( nevents, 1 );
    mR                              = mphs;
    rPyr                            = cell( nevents, 1 );
    rInt                            = cell( nevents, 1 );
    nspksPyr                        = zeros( nevents, 1 );
    mphsPyr                      	= zeros( nevents, 1 );
    mRpyr                           = zeros( nevents, 1 );
    nspksInt                        = zeros( nevents, 1 );
    mphsInt                         = zeros( nevents, 1 );
    mRint                           = zeros( nevents, 1 );
    for i                           = 1 : nevents
        ii                          = find( r{ i } );
        [ mphs( i ), mR( i ) ]      = circ_mean( phs( ii ) ); 
        rPyr{ i }                   = r{ i }( :, pyridx );
        nspksPyr( i )               = numel( find( rPyr{ i } ) );
        ii                          = find( rPyr{ i } );
        [ mphsPyr( i ), mRpyr( i ) ]    = circ_mean( phspyr( ii  ) ); 
        rInt{ i }                   = r{ i }( :, intidx );
        nspksInt( i )               = numel( find( rInt{ i } ) );
        ii                          = find( rInt{ i } );
        [ mphsInt( i ), mRint( i ) ]    = circ_mean( phsint( ii ) ); 
    end
    
    % note this is interesting - U-shaped correlation
    % i.e. synchrony is high at low and high spike counts but not in the
    % middle. Should correlate with ripple amplitude!
    % i.e. correlate population synchrony w/ ripple amp!!
    % of course the way to do it is by PYR synchrony (becasue INT are at a
    % very different pahse..)
    
    if sum( nspksPyr == median( nspksPyr ) ) >= 3
        % select the 3 best PYR synchrony, median PYR spike count events:
        fidx                        = find( nspksPyr == median( nspksPyr ) );
        [ ~, idx ]                  = sort( mRpyr( fidx ), 'descend' ); 
        maxidx                      = sort( fidx( idx( 1 : 3 ) ) );
    else
        % otherwise just take the 3 max PYR spike count ones:
        [ ~, maxidx ]               = sort( nspksPyr, 'descend' );
    end
    
    fig( 4 )                        = figure;
    for i                           = 1 : 3
        j                           = maxidx( i );
        subplot( 3, 3, i )                                                  % example raster (single ripple event)
        r0                          = r{ j };
        r0( :, pidx == 1 )          = 0;
        r1                          = r{ j };
        r1( :, pidx == 0 )          = 0;
        plot_raster( r0, phsBins / ( 2 * pi ), [], [], 2, colors( 1, : ) );
        plot_raster( r1, phsBins / ( 2 * pi ), [], [], 2, colors( 2, : ) );
        alines( cyclims( 1 ) : cyclims( 2 ), 'x', 'color', blackColor, 'linestyle', sepStyle );
        alines( find( diff( shankclu( :, 1 ) ) ) + 0.5, 'y',  'color', blackColor, 'linestyle', sepStyle );
        xlim( cyclims )
        if doLFP
            line( ptphs, scaleto( mydetrend( xlfphat( :, j ) ), ylim ), 'color', blackColor );
        end
        xlabel( 'Cycle' ), ylabel( 'Unit' )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        title( sprintf( '%s, channel %d; event %d (%d spikes)', replacetok( filename, '\_', '_' ), channel, j, nspks( j ) ) )
        axis square
        set( gca, 'tickdir', 'out', 'box', 'off' )
    end
    
    % polar plot
    pTH                             = 0.05;
    
    subplot( 3, 3, 4 )
    phi                             = mPhase; 
    R                               = rPhase;
    [ x, y ]                        = pol2car( R, phi );
    mR                              = max( R );
    sR                              = sort( R );
    m                               = ceil( mR / 0.1 );
    lims                            = [ -1 1 ] * 0.1 * m;
    lh                              = zeros( m + 2, 1 );
    for i                           = 1 : m
        lh( i )                     = circ( 0, 0, i * 0.1, [ 0 0 0 ] ); 
    end
    xlim( lims )
    ylim( lims )
    axis square
    lh( i + 1 )                     = line( lims, [ 0 0 ] );
    lh( i + 2 )                     = line( [ 0 0 ], lims );
    set( lh, 'linestyle' , '--', 'color', [ 0 0 0 ] )
    hold on
    title( sprintf( 'R: %0.3g', mR ) )
    mR                              = zeros( 1, 2 );
    sidx                            = pval( : ) <= pTH;
    for ct                          = [ 1 0 ]
        cidx                        = ~sidx & shankclu( :, 3 ) == ct;
        if sum( cidx )
            ph                      = plot( x( cidx ), y( cidx ), 'o' );
            set( ph, 'color', colors( ct + 1, : ) )
        end
        cidx                        = sidx & shankclu( :, 3 ) == ct; % plot only sig.
        if ~sum( cidx )
            continue
        end
        ph                          = plot( x( cidx ), y( cidx ), '.' );
        set( ph, 'color', colors( ct + 1, : ) )
        [ mm, mR( ct + 1 ) ]        = circ_mean( colvec( phi( cidx ) ), colvec( R( cidx ) ) ); % take the modulation depths into account
        [ xx, yy ]                  = pol2car( max( lims ), mm );
        line( [ 0 xx ], [ 0 yy ], 'color', colors( ct + 1, : ) );
    end
    axis off
    
    
    MS1                             = 12;
    MS2                             = 12;
    
    subplot( 3, 3, 5 )
    ph                              = plot( eprops( :, 1 ), nspksInt, '.b', eprops( :, 1 ), nspksPyr, '.r' ); 
    set( ph( 1 ), 'markersize', MS1, 'color', colors( 1, : ) ), 
    set( ph( 2 ), 'markersize', MS2, 'color', colors( 2, : ) )
    axis tight
    axis square 
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Power [SD]' )
    ylabel( 'Count' )
    [ cc, pp ]                      = calc_spearman( [ nspksInt( : ) nspksPyr( : ) ], eprops( :, 1 ) * [ 1 1 ], 100 );
    title( sprintf( 'INT: %0.2g (%0.2g); PYR: %0.2g (%0.2g)', cc( 1 ), pp( 1 ), cc( 2 ), pp( 2 ) ) );

    subplot( 3, 3, 6 )
    ph                              = plot( eprops( :, 2 ), nspksInt, '.b', eprops( :, 2 ), nspksPyr, '.r' ); 
    set( ph( 1 ), 'markersize', MS1, 'color', colors( 1, : ) ), 
    set( ph( 2 ), 'markersize', MS2, 'color', colors( 2, : ) )
    axis tight
    axis square 
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Freq [Hz]' )
    ylabel( 'Count' )
    [ cc, pp ]                      = calc_spearman( [ nspksInt( : ) nspksPyr( : ) ], eprops( :, 2 ) * [ 1 1 ], 100 );
    title( sprintf( 'INT: %0.2g (%0.2g); PYR: %0.2g (%0.2g)', cc( 1 ), pp( 1 ), cc( 2 ), pp( 2 ) ) );
    
    subplot( 3, 3, 7 )
    ph                              = plot( nspksInt, mRint, '.b', nspksPyr, mRpyr, '.r' ); 
    set( ph( 1 ), 'markersize', MS1, 'color', colors( 1, : ) ), 
    set( ph( 2 ), 'markersize', MS2, 'color', colors( 2, : ) )
    axis tight
    axis square 
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( '#spikes' )
    ylabel( 'Sync' )
    [ cc, pp ]                      = calc_spearman( [ mRint( : ) mRpyr( : ) ], [ nspksInt( : ) nspksPyr( : ) ], 100 );
    title( sprintf( 'INT: %0.2g (%0.2g); PYR: %0.2g (%0.2g)', cc( 1 ), pp( 1 ), cc( 2 ), pp( 2 ) ) );
    
    subplot( 3, 3, 8 )
    ph                              = plot( eprops( :, 1 ), mRint, '.b', eprops( :, 1 ), mRpyr, '.r' ); 
    set( ph( 1 ), 'markersize', MS1, 'color', colors( 1, : ) ), 
    set( ph( 2 ), 'markersize', MS2, 'color', colors( 2, : ) )
    axis tight
    axis square 
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Power [SD]' )
    ylabel( 'Sync' )
    [ cc, pp ]                      = calc_spearman( [ mRint( : ) mRpyr( : ) ], eprops( :, 1 ) * [ 1 1 ], 100 );
    title( sprintf( 'INT: %0.2g (%0.2g); PYR: %0.2g (%0.2g)', cc( 1 ), pp( 1 ), cc( 2 ), pp( 2 ) ) );

    subplot( 3, 3, 9 )
    ph                              = plot( eprops( :, 2 ), mRint, '.b', eprops( :, 2 ), mRpyr, '.r' ); 
    set( ph( 1 ), 'markersize', MS1, 'color', colors( 1, : ) ), 
    set( ph( 2 ), 'markersize', MS2, 'color', colors( 2, : ) )
    axis tight
    axis square 
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Freq [Hz]' )
    ylabel( 'Sync' )
    [ cc, pp ]                      = calc_spearman( [ mRint( : ) mRpyr( : ) ], eprops( :, 2 ) * [ 1 1 ], 100 );
    title( sprintf( 'INT: %0.2g (%0.2g); PYR: %0.2g (%0.2g)', cc( 1 ), pp( 1 ), cc( 2 ), pp( 2 ) ) );

    
end

if ~isempty( savetype ) && ( isa( savetype, 'cell' ) || ~all( isnan( savetype ) ) ) && ~isempty( figname )
    if ~isa( savetype, 'cell' )
        savetype                    = { savetype };
    end
    for j                           = 1 : length( savetype )
        for i                       = 1 : length( fig )
        if fig( i ) == 0
            continue
        end
            fignameI                = [ figname '.part' num2str( i ) '.' savetype{ j } ];
            fprintf( '%s: Saving figure %s\n', upper( mfilename ), fignameI )
            fig_out( fig( i ), 1, fignameI, savetype{ j } );
        end
    end
end
verb( sprintf( '%s: DONE!\n', upper( mfilename ) ), verbose )

return

% EOF

