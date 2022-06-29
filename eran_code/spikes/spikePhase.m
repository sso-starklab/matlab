% spikePhase      compute spiking phase histograms for multiple frequencies
% 
% [ hRate, phsBins, freqBins, H, T, B ] = spikePhase( filebase, shankclu, varargin )
%
% DOES
% computes the spiking phase-frequency maps for multiple units
%
% ARGUMENTS
% filebase          full path + base or par file
% shankclu          list of units (or shank numbers; see get_spikes.m)
%
% ADDITIONAL ARGUMENTS (given as parameter name/value pairs):
% phaseSource       determines the source of the phase (which channel): 
%                   By default, the phase is taken by averaging over all 
%                       eeg channels of the first electrode group in shankclu. 
%                       Stimulus times are excluded. 
%                       This is specified by {'spontaneous'}. 
%                   Alternatives:
%                   -'ripple': take the channel with the highest ripple
%                       power (of the same electrode group). This requires an
%                       *sps* file, see detect_all_hfos.m
%                   -a list of one or more channel numbers. Then the phase
%                       is taken by averaging over those channels.
%                   -a cell array with parameters to be passed to
%                       stim_select.m. The channel to be use is the one
%                       corresponding to the electrode group (or specified
%                       explicitly using an additional argument 'trigchan')
% trichan           {[]}; specifies a trigger channnel (instead of the one
%                       corresponding to the electrode group)
% periods           {[]}; specifies the periods to be used. The default is
%                       to use all relevant periods (e.g. 'spontaneous',
%                       'ripple', or specified - all stim-free). Options:
%                   -a 2-column matrix of ranges (at the spike Fs)
%                   -a brain state string (see get_states.m)
% spectralMethod    determines the spectral analysis parameters. The
%                       foramt is { fMin fMax nfBins method }, and the 
%                       default is { 2 300 65 'wavelet' }. Alternative
%                       methods are 'linear' and 'log', and then a filter
%                       bank is used (see spikePhaseFreq.m)
% phases            {20} the number of phase bins to be used. Alternatively
%                       a vector of edges ([rad], spanning the unit circle) 
%                       may be specified
% powTH             {2} [SD]. Number of SDs for each frequency. Phase is
%                       evaluated for all valid, but spikes are counted
%                       only for those segments of sufficiently high power.
%                       The baseline (mean, SD) is determined separately
%                       for each frequency bin from the SWS data (see
%                       segmentBehavior and computeBaseline)
% suffix            {'eeg'}. Alternatively 'dat' for wide-band phase
% ilevel            {B'}. ignored if shankclu specifies units; if shankclu
%                       specifies an electrode group, class 'B' units will
%                       be used by default
% Overwrite         1 to recompute and overwrite
%                   0 to just compute (with writing but not overwriting)
%                   -1 to just load/compute if a file exists/doesn't (no writing)
%                   {-2} load if exists, compute and save if doesn't
% savef             {1}; defaults to ../../figs/ directory, otherwise
%                       specify the path
% savetype          {'png'}; format to save figure as
%
% OUTPUT
% hRate             [spikes/sec]; 3D array, phase x freq x unit
%                       easily reconstructed from H: 
%                       h( f ) = H( f ) * n / T( f )
%                       where n = length( phsBins )
% phsBins           [rad]; the phase BINS (not edges)
% freqBins          [Hz]; the frequency bins
% H                 [count] the raw count matrix; same dimensions as hRate
% T                 [sec] the total time spent in each frequency
% B                 [spikes/sec]; the "baseline" rate of each unit (defined
%                       as the mean rate during relevant epochs)
%
% FILES
% requires          *eeg/*dat, *clu*, *res*; *.xml, *stm*
% optional:         *sps*
% intermediate:     *wltBL*/*firBL* 
% output            *wlt.sph*/*fir.sph*
%
% NOTE
% This routine is intended to be ran for one electrode group at a time,
%   thus shankclu should include only units from one shank, so the eeg phase
%   is local. If another source is desired/same source for multi-shank units,
%   the source channel/s have to specified explicitly. 
%                   
% calls:            LoadXml, ParseArgPairs, verb, num3str
%                   get_stimchans, LoadStims, stim_select, get_states
%                   determine_units, get_spikes, get_merged_filenum
%                   makeblocks, sortranges, resampleranges, inranges, setdiffranges
%                   spikePhaseFreq, spikePhasePlot, computeBaseline, plotSpectrogram

% 13-mar-13 ES

% revisions:
% 17-mar-13 enable external determination of periods, or specify a brain
%               state to use
% 24-mar-13 added output sourcePeriods
% 04-apr-13 added output shankclu
%           extended support for externally specified states
%           modified save name to include 3-character state string
%           changed powTH default to zero (use segmentBehavior epochs typically)
% 07-apr-13 call to spikePhasePlot
%           added file information (spectrum, total time spent in state)
% 21-oct-14 modified defaultfigdir
% 17-aug-19 modified call to spikePhaseFreq during computations to account for imprecision of firFreqs
% 18-aug-19 cleaned up

% ISSUES:
% (1) there is an issue with the segmentation. If the brain state is not
%   variable, e.g. theta throughout, then the theta freqeuncies will be
%   under-sampled. Thus the SD method is appropriate for high frequencies (where
%   there is much variability) but much lees for theta/spindle. This can
%   generate very different plots for powTH 0 and 2 (at the 0-50 Hz range,
%   mainly; see e.g. m531r1_40 3.35
% (2) statistical significance... (randomly permute the spike times between
% segments, but this should be done for each frequency separately..)
% 
% MINOR:
% (3) keep somewhere the number of spikes ACTUALLY used (right now the
%       sum(h) is inflating the count)
% (4) speed up the computation somehow - too long when many small blocks..

% add duration on figure. Plot the cumulative spectrum (can be used from
% *phs, but that's for the entire file; for the specified epochs, should
% compute locally as the phase is computed)

function [ hRate, phsBins, freqBins, H, T, B, sourcePeriods, shankclu ] = spikePhase( filebase, shankclu, varargin )

hRate                       = [];
phsBins                     = [];
freqBins                    = [];
H                           = [];
T                           = [];
B                           = [];
sourcePeriods               = [];

%------------------------------------------------------------------%
% globals
%------------------------------------------------------------------%
global BatchMode
if isempty( BatchMode )
    BatchMode               = 0;
end

%------------------------------------------------------------------------%
% constants
%------------------------------------------------------------------------%
verbose                     = 1;
nRandSegments               = 1000;               % [segments]
BLOCKSIZE                   = 2^20;                   % [samples]

%------------------------------------------------------------------------%
% arguments
%------------------------------------------------------------------------%
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    error( 'missing arguments' )
end
[ phaseSource, trigchan, periods, spectralMethod, phases, powTH...
    , specMode, NW...
    , suffix, clustr, ilevel, Overwrite, graphics, savef, savetype ] = ParseArgPairs(...
    { 'phaseSource', 'trigchan', 'periods', 'spectralMethod', 'phases', 'powTH'...
    , 'specMode', 'NW'...
    , 'suffix', 'clustr', 'ilevel', 'Overwrite', 'graphics', 'savef', 'savetype' }...
    , { 'spontaneous', [], [], { 2 300 65 'wavelet' }, 20, 0 ...
    , 'welch', 3 ...
    , 'eeg', 'clu', 'B', -2, 1, 1, 'png' }...
    , varargin{ : } );

if isa( filebase, 'struct' ) && isfield( filebase, 'FileName' )
    par                     = filebase;
    filebase = par.FileName;
elseif isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
    par                     = LoadXml( filebase );
end

% get some general parameters
eegFs                       = par.lfpSampleRate;
spkFs                       = par.SampleRate;
nchans                      = par.nChannels;
switch suffix
    case 'eeg'
        Fs                  = eegFs;
    case 'dat'
        Fs                  = spkFs;
    otherwise
        error( 'erroneous specification of suffix' )
end

% parse the spectral parameters
if length( spectralMethod ) ~= 4 || ~isa( spectralMethod, 'cell' )
    error( 'erroneous format of spectralMethod' )
end
fMin                        = spectralMethod{ 1 };
fMax                        = spectralMethod{ 2 };
nfBins                      = spectralMethod{ 3 };
method                      = lower( spectralMethod{ 4 } );
switch method
    case { 'log', 'linear', 'spikewb' }
        specstr             = 'fir';
    case 'wavelet'
        specstr             = 'wlt';
    otherwise
        error( 'unrecognized methog' )
end
        
% parse the phaseSource parameter
if isa( phaseSource, 'numeric' ) && isvector( phaseSource ) && all( inrange( phaseSource, [ 1 nchans ] ) )
    chans                   = phaseSource;
    phaseSource             = 'specified';
elseif isa( phaseSource, 'char' ) && ismember( lower( phaseSource ), { 'spontaneous', 'ripple', 'theta' } )
    phaseSource             = lower( phaseSource );
elseif isa( phaseSource, 'cell' ) && iseven( length( phaseSource ) )
    stimParams              = phaseSource;
    phaseSource             = 'stimulation';
    powTH                   = 0;
    fprintf( '%s: Better use wnAnalysis.m for this purpose!\n', upper( mfilename ) )
    keyboard
else
    error( 'erroneous format for phaseSource' )
end

% determine the shankclu automatically
if isempty( shankclu )
   switch phaseSource
       case { 'theta', 'specified' }
           shankclu         = 1 : length( par.SpkGrps );
       otherwise % 'ripple', 'spontaneous' - must specify egroup
           error( 'shankclu must be specified for phaseSource %s', phaseSource )
   end
end
if ~isa( shankclu, 'numeric' ) || ~isvector( shankclu ) ...
        && ~inrange( size( shankclu, 2 ), [ 2 3 ] ) ...
        || ~ismatrix( shankclu ) 
    error( 'erroneous format of shankclu' )
end

% parse the periods parameter
if isempty( periods )
    statename               = 'all';
else
    err                     = 0;
    forcePeriods            = 0;
    if isa( periods, 'numeric' ) && size( periods, 2 ) 
        statename           = 'ext';
        periods             = sortranges( periods );
        forcePeriods        = 1;
        err                 = 0;
    elseif isa( periods, 'char' ) || isa( periods, 'cell' )
        statename           = periods;
        [ periods, msg ]    = get_states( filebase, statename );
        if ~isempty( msg )
            err             = 1;
        end
        periods             = resampleranges( periods, spkFs, eegFs );
        if isa( statename, 'cell' )
            statename       = statename{ 1 };
        end
    else
        err                 = 1;
    end
    if err
        error( 'erroneous format for periods' )
    end
    if isempty( periods )
        fprintf( '%s: No relevant periods for %s\n', upper( mfilename ), filebase )
        return
    end
end

if ~isempty( trigchan )
    trigchan                = trigchan( 1 );
end

if powTH < 0 || length( powTH ) > 1
    error( 'erroneous specification of powTH' )
end

%------------------------------------------------------------------------%
% preparations
%------------------------------------------------------------------------%

% determine the phases
if length( phases ) == 1
    npBins                  = abs( round( phases ) );
    phases                  = ( -pi : 2 * pi / npBins : pi )';
else
    phases                  = phases( : );
    npBins                  = length( phases ) - 1;
end
if length( unique( mod( phases( [ 1 end ] ), 2 * pi ) ) ) ~= 1
    error( 'Phase edges should cover the entire unit circle' )
end
phases                      = phases( : );

% determine the frequencies 
switch method
    case { 'linear', 'spikewb' }
        freqs               = linspace( fMin, fMax, nfBins );
    case { 'log', 'wavelet' }
        freqs               = logspace( log10( fMin ), log10( fMax ), nfBins );
end

% get the stimulus times
[ vals, allstimchans, allstims ] = LoadStims( filebase );

% map the source file
sourcefname                 = [ filebase '.' suffix ];
a                           = memmapfile( sourcefname, 'Format', 'int16' );
neeg                        = length( a.Data ) / nchans;

% select the channels + periods for analysis
switch phaseSource
    case { 'spontaneous', 'ripple', 'specified' }
        egroup              = min( shankclu( :, 1 ) );
        switch phaseSource
            case 'theta'
                if exist( [ filebase '.phs' ], 'file' )
                    L       = load( [ filebase '.phs' ], 'eegchan', '-mat' );
                    chans   = L.eegchan;
                else
                    error( 'missing phs file' )
                end
                chanstr     = num3str( chans, 3 );
            case 'spontaneous' % average over all channels
                chans       = par.ElecGp{ egroup } + 1;
                chanstr     = num2str( egroup );
            case 'ripple' % take the ripple channel
                if exist( [ filebase '.sps' ], 'file' )
                    L       = load( [ filebase '.sps' ], '-mat' );
                    chans   = L.vote( egroup );
                else
                    error( 'missing sps file' )
                end
                chanstr     = num3str( chans, 3 );
            case 'specified'
                chanstr     = num3str( chans( 1 ), 3 );
        end
        voltagerange        = par.VoltageRange;
        medx                = 0;
        notvalid            = sortranges( vals( :, 1 : 2 ) );
        if isempty( periods )
            periods         = setdiffranges( [ 1 neeg / Fs * spkFs ], notvalid );
        elseif ~forcePeriods
            periods         = setdiffranges( periods, notvalid );
        end
    case 'stimulation'
        [ stimchans, targets, voltageranges ] = get_stimchans( par );
        if isempty( trigchan )
            egroup          = min( shankclu( :, 1 ) );
            cidx            = targets == egroup;
        else
            cidx            = stimchans == trigchan;
            if sum( cidx ) == 0
                error( 'erroneous specification of trigchan' )
            end
        end
        chans               = stimchans( cidx );
        voltagerange        = voltageranges( cidx );
        stim                = allstims( allstimchans == chans );
        medx                = stim.median;
        out                 = stim_select( stim, stimParams );
        if isempty( periods )
            periods         = out.times;
        else
            periods         = intersectranges( periods, out.times );
        end
        chanstr             = num3str( chans, 3 );
end

% determine which units to use and their properties
if size( shankclu, 2 ) < 2 || size( shankclu, 2 ) > 3
    [ shankclu, sst ]       = determine_units( filebase, shankclu, ilevel );
else
    load( [ filebase '.sst' ], '-mat' );
end
if isempty( shankclu )
    fprintf( '%s: No relevant units for %s\n', upper( mfilename ), filebase )
    return
end
sidx                        = ismember( sst.shankclu, shankclu( :, 1 : 2 ), 'rows' );
loc                         = [ sst.shankclu( sidx, 1 ) round( sst.geo_com( sidx ) ) ];
emap                        = get_emap( par );
[ ~, idx ]                  = ismember( loc, emap( :, 1 : 2 ), 'rows' );
shankcluExpanded            = [ shankclu emap( idx, 3 ) ];

% determine saving/ploting
[ pathname, filename, extname ] = fileparts( filebase );
filename                    = [ filename extname ];
fname                       = sprintf( '%s.%s.sph.%s.%s', filename, specstr, statename, chanstr );
savename                    = [ pathname '/' fname ];
if strcmpi( method, 'spikewb' )
    savename                = [ savename '_' method ];
    spectralMethod{ 4 }     = 'linear';
end
homedir                     = strfind( pathname, 'dat' );
if isempty( homedir )
    homedir                 = [ fileparts( pathname ) '/' ];
else
    homedir                 = pathname( 1 : homedir - 1 );
end
defaultfigdir               = sprintf( '%sfigs', homedir );
if ~graphics
    savef                   = 0;
end
if isequal( savef, 1 )
    if ~exist( defaultfigdir, 'dir' )
        mkdir( homedir, 'figs' )
    end
    figname                 = [ defaultfigdir '/' fname ];
elseif isa( savef, 'char' ) && exist( savef, 'dir' )
    figname                 = [ savef '/' fname ];
else
    figname                 = fname;
end
if Overwrite < 0 && exist( savename, 'file' )
    verb( sprintf( '%s: Loading %s..', upper( mfilename ), savename ), -verbose )
    load( savename, '-mat' );
    verb( sprintf( ' done!' ), verbose )
    kidx                    = ismember( shankcluExpanded( :, 1 : 2 ), shankclu( :, 1 : 2 ), 'rows' );
    shankcluExpanded        = shankcluExpanded( kidx, : );
    if graphics > 0
        if ~exist( 'totT', 'var' )
            totT = sum( diff( sourcePeriods, 1, 2 ) + 1  ) / Fs;
        end
        if ~exist( 'fP', 'var' ) || ~exist( 'P', 'var' )
            fP = [];
            P = [];
        end
        spikePhasePlot( hRate, B, phsBins, freqBins, shankcluExpanded, [ fP P ], totT, figname, savetype );
    end
    return
end

% get the filter bank
if ismember( method, { 'log', 'linear', 'spikewb' } )
    [ ~, ~, ~, freqsHat, ~, ~, ~, ~, ~, hBP ] = spikePhaseFreq( 1 ...
        , 1, 1, 'freqs', spectralMethod, 'compute', 0, 'eegFs', eegFs );
    mindur                  = length( hBP{ 1 } ) / Fs;
else
    mindur                  = 1 / fMin;
end

% keep only periods of sufficient duration
durs                        = ( diff( periods, [], 2 ) + 1 ) / spkFs; % [sec]
ridx                        = durs < mindur;
periods( ridx, : )          = [];
durs( ridx, : )             = [];
sourcePeriods               = resampleranges( periods, Fs, spkFs );
verb( sprintf( '%s: %d periods kept for analysis; %d sec'...
    , upper( mfilename ), size( periods, 1 )...
    , ceil( sum( durs ) ) ), verbose )

% get the scaling factor
scalefactor                 = 1 / 2 .^ par.nBits * voltagerange / par.Amplification * 1000; % [ mV ]

% determine the SWS baseline power
if powTH > 0
    mfilebase               = get_merged_filenum( filebase ); % work with the merged file
    basefile                = sprintf( '%s.%sBL.%s', mfilebase, specstr, chanstr );
    verb( sprintf( '%s: Loading baseline data (%s)', upper( mfilename ), basefile ), verbose )
    if exist( basefile, 'file' )
        L                   = load( basefile, '-mat' );
    end
    if ~exist( basefile, 'file' ) || nfBins ~= length( L.freqs )...
            || ~inrange( L.freqs( 1 ), fMin * [ 0.95 1.05 ] )...
            || ~inrange( L.freqs( end ), fMax * [ 0.95 1.05 ] )%...
        if ismember( method, { 'log', 'linear', 'spikewb' } )
            computeBaseline( mfilebase, chans, 'fir', 'hBP', hBP...
                , 'scalefactor', scalefactor, 'nrand', nRandSegments...
                , 'Overwrite', 1, 'chanstr', chanstr );
        else
            computeBaseline( mfilebase, chans, 'wlt', 'fMin', fMin, 'fMax', fMax, 'nfBins', nfBins...
                , 'scalefactor', scalefactor, 'nrand', nRandSegments...
                , 'Overwrite', 1, 'chanstr', chanstr );
        end
        L                   = load( basefile, '-mat' );
    end
    swsM                    = L.mm;
    swsSD                   = L.ss;
else
    swsM                    = [];
    swsSD                   = [];
end

% load the spikes
[ clu, res, map, shankclu ] = get_spikes( par, shankclu, clustr );
kidx                        = inranges( res, periods );
res                         = res( kidx );
clu                         = clu( kidx );
uclu                        = unique( clu );
ridx                        = ~ismember( map( :, 1 ), uclu );
map( ridx, : )              = []; % consider only active cells
shankclu( ridx, : )         = [];
nclu                        = length( uclu );

%------------------------------------------------------------------------%
% load the eeg/dat piecewise and accumulate results
%------------------------------------------------------------------------%

blocksize                   = 2 ^ floor( log2( BLOCKSIZE / length( chans ) ) );
pidx                        = diff( sourcePeriods, [], 2 ) >= blocksize;
if sum( pidx )
    blocks                  = sourcePeriods( ~pidx, : );
    pidx                    = find( pidx );
    for i                   = 1 : length( pidx )
        seg                 = makeblocks( diff( sourcePeriods( pidx( i ), : ) ) + 1, blocksize, 0 );
        blocks              = [ blocks; seg + sourcePeriods( pidx( i ), 1 ) - 1 ];
    end
    blocks                  = sortranges( blocks, 0 );
else
    blocks                  = sourcePeriods;
end
nblocks                     = size( blocks, 1 );
spkBlocks                   = resampleranges( blocks, spkFs, Fs );
verb( sprintf( '%s: Partitioned into %d blocks: BlockNum | nSpikes | nSpikesInMat: ' ...
    , upper( mfilename ), nblocks ), verbose )
H                           = zeros( npBins, nfBins, nclu );
C                           = zeros( nclu, 1 );
T                           = zeros( nfBins, 1 );

% initialize local spectrum estimation
nFFT                        = 2^floor( log2( Fs ) );
pxxBlock                    = zeros( ceil( ( nFFT + 1 ) / 2 ), nblocks );

% go over blocks
for bidx                    = 1 : nblocks
    
    if verbose
        if ~mod( bidx, 10 )
            verb( sprintf( '\t\t%d | ', bidx ), -verbose )
        end
        if bidx == nblocks
            fprintf( 1, '\n' ); 
        end
    end
    
    % load the data
    sidx                    = inranges( res, spkBlocks( bidx, : ) );
    club                    = clu( sidx );
    resb                    = res( sidx ) - spkBlocks( bidx, 1 ) + 1;
    if isempty( sidx )
        continue
    end
    n                       = diff( blocks( bidx, : ) ) + 1;
    m1                      = ones( n, 1 ) * ( ( blocks( bidx, 1 ) * nchans - nchans + chans ) );
    m2                      = ( ( 0 : ( n - 1 ) )' * nchans ) * ones( 1, length( chans ) );
    widx                    = m1 + m2;
    xb                      = ( single( a.data( widx ) ) - medx ) * scalefactor;
    xb                      = mean( xb, 2 );

    % call the routine
    if strcmp( method, 'wavelet' )
        [ h, freqDurs, phsBins, freqBins, b ] = spikePhaseFreq( club, resb, xb...
            , 'M', swsM, 'SD', swsSD, 'phases', phases, 'freqs', spectralMethod, 'powTH', powTH...
            , 'verbose', 0, 'spkFs', spkFs, 'eegFs', Fs );
    elseif ismember( method, { 'log', 'linear' } )
        [ h, freqDurs, phsBins, freqBins, b ] = spikePhaseFreq( club, resb, xb...
            , 'M', swsM, 'SD', swsSD, 'phases', phases, 'freqs', { hBP, freqsHat }, 'powTH', powTH...
            , 'verbose', 0, 'spkFs', spkFs, 'eegFs', Fs );
    else strcmpi( method, 'spikewb' )
        [ s, aux ] = spikeWB( club, resb, xb, 'periods', [ 1 length( xb ) ]...
            , 'xperiods', [ 1 length( xb ) ], 'freqs', spectralMethod, 'phases', phases...
            , 'M', 1, 'spkFs', spkFs, 'Fs', Fs, 'nreps', [ 10 0 10 0 ] );
        h = s.phsHistRate; % need to multiply by nCycles etc (rate->count)
        % but there is a difference between the s.phsHistRate and the
        % hhRate output from spikePhaseFreq - check M/SD clipping
        % Also, not clear to me how to accumulate other parameters, in
        % particular the coherence... can accumulate the complex spectra
        % using a weighted summation, but this requires special programming
        % (i.e. cannot use spikeWB as is :( .. 
        b                   = s.frs( : );
        phsBins             = s.phsBins;
        freqBins            = s.freqBins;
        freqDurs            = length( xb ) / Fs * ones( length( freqBins ), 1 );
        fprintf( '%s: FINISH PROGRAMMING THIS: (1) coh/phs accumulation (2) spike maps accumulation\n', upper( mfilename ) )
        keyboard
    end
    
    % compute the spectrum of the eeg
    eeg                     = xb - mean( xb );
    switch specMode
        case 'welch'
            [ pxxBlock( :, bidx ), fP ] = my_spectrum( eeg, nFFT, Fs, nFFT, nFFT / 2, 'none' );
        case 'mt' % ~nTapers longer (linear)
            [ pxxBlock( :, bidx ), fP ] = mtcsd( eeg, nFFT, Fs,  nFFT, nFFT / 2, NW, '' );
    end
    
    if verbose && ~mod( bidx, 10 )
        verb( sprintf( '%d | %d ', length( resb ), sum( h( : ) ) ), verbose )
    end
    
    % accumulate
    uidx                    = ismember( uclu, unique( club ) );
    T                       = T + freqDurs;
    H( :, :, uidx )         = H( :, :, uidx ) + h;
    C( uidx )               = C( uidx ) + b * n / Fs;
    
    % TO DO: also accumulate the structures from spikeWB
    
end

if nblocks == 0
    freqBins                = freqs; 
    phsBins                 = mean( [ phases( 1 : end - 1 ) phases( 2 : end ) ], 2 );
end
clear a
nCycles                     = T .* freqBins( : );
phaseBinSize                = 1 ./ freqBins( : ) / npBins;
hRate                       = zeros( size( H ) );
for j                       = 1 : nfBins
    hRate( :, j, : )        = H( :, j, : ) / nCycles( j ) / phaseBinSize( j );
end
durs                        = diff( blocks, [], 2 ) + 1;
B                           = C / sum( durs ) * Fs;
if isempty( bidx )
    P                       = [];
    fP                      = [];
else
    P                       = sum( bsxfun( @times, pxxBlock, durs' ), 2 ) / sum( durs );
end
totT                        = sum( durs / Fs );

%------------------------------------------------------------------------%
% summarize: plot and save
%------------------------------------------------------------------------%
% save data
if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( savename, 'file' ) )
    verb( sprintf( '%s: Saving %s', upper( mfilename ), savename ), verbose )
    generator               = { computer, datestr( now, 'ddmmmyy' ) };
    save( savename, 'filebase', 'suffix', 'Fs', 'scalefactor', 'generator'...
        , 'chans', 'shankclu'...
        , 'sourcePeriods', 'powTH'...
        , 'totT', 'fP', 'P' ...
        , 'hRate', 'phsBins', 'freqBins', 'H', 'T', 'B', '-v6' );
end

if graphics
    [ ~, kidx ]         = ismember( shankclu( :, 1 : 2 ), shankcluExpanded( :, 1 : 2 ), 'rows' );
    shankcluExpanded    = shankcluExpanded( kidx, : );
    spikePhasePlot( hRate, B, phsBins, freqBins, shankcluExpanded, [ fP P ], totT, figname, savetype )
end

if BatchMode
    close all
end

return

% EOF

