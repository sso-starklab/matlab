% computeBaseline       spectrum by randomization over SWS segments
%
% [ mm, ss, freqs ] = computeBaseline( filebase, specChans, specMode )
%
% DOES
% computes the spectral power during SWS. this is done by randomly
% selecting NRAND independent segments, each of SNIPDUR, computing the
% spectrum for each, and averaging over segments. 
%
% ARGUMENTS:
% filebase              full path + base or par file
% specChans             channels for which to compute spectrum. if
%                           multiple, data are averaged over channels
%                           BEFORE spectral estimation
% specMode              {'wlt'}; spectral method: wavelet 'wlt',
%                           multi-taper 'spc', or filter bank 'fir'
%                       
% ADDITIONAL ARGUMENTS:
% nrand                 {1000}, [count] max number of independent segments to randomize
% windur                {0.1} [sec]; duration (exact) of each segment
% suffix                {'eeg'} or 'dat'; of source data file
% chans                 {specChans} additional channels to consider (relevant for CSD computation only)
% Overwrite             1 to recompute and overwrite
%                       0 to just compute (with writing but not overwriting)
%                       -1 to just load/compute if a file exists/doesn't (no writing)
%                       {-2} load if exists, compute and save if doesn't
%
% epochs                {[]}; to specify other epochs
% doCSD                 {0}; to do the CSD; requires at least 3 channels
% refchan               {[]}; subtract a reference channel before computations
% whitenFlag            {0}; to prewhiten data (distorts lower frequencies)
% 
% spectral parameters:
%   'wlt' requires fMin, fMax, nfBins; defaulted to 2,300,65
%   'spc' required M, spectral resolution; defaulted to 1 (~approx 1 Hz)
%   'fir' requires hBP, the fir filter bank, defaulted to 65 filters
%           log-spaced between 2 and 300
%
%
% OUTPUT
% mm, ss                mean and SD in each frequency bin
% freqs                 [Hz] the frequency bins
% 
% ADDITIONAL OUTPUT:
% ssq                   sum of squares, to derive accumulated variance:
%                           ss = sqrt( ssq - mm.^2 )
% nrand                 the actual number of segments used
%
% 
% NOTE
% 1. The data size used depends on windur/nrand but also on fMin and specMode.
% For 'fir', the minimum windur is determined by the longest fir length,
% and for other methods by the lowest frequency. In any case, the window is
% doubled to prevent edge effects
% For instance, for 'fMin = 2, the actual segments will be 6 sec for 'fir'
% and 1 sec for other methods; of these, half will be kept for power
% averaging. 
% 2. There is not temporal or spatial smoothing in this routine (relevant
% for both CSD and LFP). This is critical, temporal smoothing will distort
% the estimate of spiking bleed-through and spatial smoothing may dampen
% HFOs considerably
%
% calls                 ParseArgPairs, LoadXml, makefir, local_max
%                       segmentBehavior, LoadStims, setdiffranges
%                       readbin, WhitenSignal, mtcsg, getWavelet, firfilt, ma_rms
%
% see also              pt_avg, spikePhase

% 12-mar-13 ES

% revisions
% 04-apr-13 (1) firFreqs replaces local extraction of frequencies from fir
%           (2) maximal snipdur for fir is now same as the longest hBP
%           (3) ssq saved as well
% 07-apr-13 (1) CSD handling 
%           (2) output info in s
%           (3) nWindow external control
% 20-oct-13 expand edges for CSD
% 24-feb-15 (1) added 'dog' option. this requires a 3-element filter bank:
%               lowpass, highpass, and power (low pass). it assumes the filter
%               is DOG (see detect_hfos)
%           (2) modified condition for re-computation - windur was assumed
%               to be in seconds; nrands/2->nrands/4
% 10-sep-18 adapatation to Stark Lab

% NOTE:
% there is a bug w/ 'fir': at the last frequency bin, there is a sudden
% rise of power... need to resolve this, but not urgent

function [ mm, ss, freqs, s ] = computeBaseline( filebase, specChans, specMode, varargin )

%------------------------------------------------------------%
% constants
%------------------------------------------------------------%
verbose                     = 1;
behaviorDurSEC              = 2;                                                                    % [s]
OverwriteStates             = -2;                                                                   % do not overwrite if already existing

%------------------------------------------------------------%
% preparations - arguments
%------------------------------------------------------------%
% initialize output
mm                          = [];
ss                          = [];
freqs                       = [];
%s = [];

nargs                       = nargin;
if nargs < 2 || isempty( filebase ) || isempty( specChans )
    error( 'filebase and specChans required' )
end
if nargs < 3 || isempty( specMode )
    specMode                = 'wlt';
end
[ nrand, windur, chans, suffix, epochs, scalefactor, Overwrite, chanstr...
    , expandEdges...
    , whitenFlag, refchan, doCSD, fMin, fMax, nfBins, M, nWindow, hBP, padBuffer ] = ParseArgPairs(...
    { 'nrand', 'windur', 'chans', 'suffix', 'epochs', 'scalefactor', 'Overwrite', 'chanstr'...
    , 'expandEdges' ...
    , 'whitenFlag', 'refchan', 'doCSD', 'fMin', 'fMax', 'nfbins', 'M', 'nWindow', 'hBP', 'padBuffer' }...
    , { 1000, 0, specChans, 'eeg', [], [], -2, [] ...
    , 1 ...
    , 0, [], 0, [], [], [], [], [], [], [ -0.01 0.01 ] }...
    , varargin{ : } );

if isa( filebase, 'struct' ) && isfield( filebase, 'FileName' )
    par                     = filebase;
    filebase                = par.FileName;
elseif isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
    par                     = LoadXml( filebase );
end
eegFs                       = par.lfpSampleRate;
spkFs                       = par.SampleRate;
nchans                      = par.nChannels;

if sum( specChans > nchans ) || sum( specChans < 1 )...
        || sum( specChans ~= round( specChans ) )
    error( 'erroneous specification of specChans' )
end
if sum( chans > nchans ) || sum( chans < 1 )...
        || sum( chans ~= round( chans ) )
    error( 'erroneous specification of chans' )
end
if ~isempty( refchan )
    if refchan > nchans || refchan < 1 || refchan ~= round( refchan )
        refchan             = [];
    end
end
switch lower( specMode )
    case { 'wlt', 'wavelet' }
        specMode            = 'wlt';
    case { 'spc', 'mt', 'multitaper' }
        specMode            = 'spc';
    case { 'fir', 'filterbank' }
        specMode            = 'fir';
    case { 'dog', 'bp' }
        specMode            = 'dog';
    otherwise
        error( 'unrecognizd mode' )
end
nrand                       = round( nrand( 1 ) );
nrand0                      = nrand;
if nrand < 1
    error( 'erroneous specification of nrand' )
end
switch suffix
    case 'eeg'
        Fs                  = eegFs;
    case 'dat'
        Fs                  = spkFs;
    otherwise
        error( 'erroneous specification of suffix' )
end
sourcefile                  = [ filebase '.' suffix ];
windur                      = windur * Fs;                                                          % [s]->[samples]
if isempty( scalefactor )
    egroup                  = get_egroup( par, chans( 1 ) );
    scalefactor             = par.AnatGrps( egroup ).VoltageRange( 1 ) ...
        / 2 .^ par.nBits / par.Amplification * 1e3;                                                 % [ mV ]
end
if isempty( chanstr )
    chanstr                 = num3str( specChans( 1 ) );
end
savename                    = sprintf( '%s.%sBL.%s', filebase, specMode, chanstr );

if Overwrite < 0 && exist( savename, 'file' )
    verb( sprintf( '%s: Loading %s..', upper( mfilename ), savename ), -verbose )
    L                       = load( savename, '-mat' );
    s                       = L;
    verb( sprintf( '%d repetitions, %0.3g sec windows', L.nrand, L.snipdur / L.Fs ), verbose )
    if ~isfield( L, 'mm' ) || ~isfield( L, 'ss' ) || ~isfield( L, 'snipdur' )...
            || ~isfield( L, 'Fs' ) || ~isfield( L, 'nrand' )...
            || L.snipdur < windur || L.nrand < nrand / 4
        if Overwrite ~= -1                                                                          % recompute and possibly rewrite 
            Overwrite       = 1;
        end
    else
        mm                  = L.mm;
        ss                  = L.ss;
        nrand               = L.nrand;
        if isfield( L, 'freqs' )
            freqs           = L.freqs;
        else
            freqs           = [];
        end
        if isfield( L, 'ssq' )
            ssq             = L.ssq;
        else
            ssq             = [];
        end
        return
    end
end
s                           = [];

%------------------------------------------------------------%
% preparations - spectral parameters
%------------------------------------------------------------%
switch specMode
    case 'wlt'
        if ~exist( 'fMmin', 'var' ) || ~isempty( fMin ) ...
                || exist( 'fMax', 'var' ) || ~isempty( fMax ) ...
                || exist( 'nfBins', 'var' ) || ~isempty( nfBins )
            fMin            = 2;
            fMax            = 300;
            nfBins          = 65;
        end
    case 'spc'
        if isempty( M )
            M               = 1;                                                                    % ~1 Hz spectral resolution
        end
        if isempty( nWindow )
            nWindow         = windur + mod( windur, 2 );
        end
        mtNW                = 3;
        dflag               = '';
        nFFT                = 2^floor( log2( Fs * M ) );
        fMin                = Fs / nFFT;
    case 'fir'
        if isempty( hBP )
            fMin            = 2;
            fMax            = 300;
            nfBins          = 65;
            method          = 'log';
            freqs           = logspace( log10( fMin ), log10( fMax ), nfBins );
            verb( sprintf( '%s: Constructing filter bank (%s spacing)...'...
                , upper( mfilename ), method ), -verbose )
            hBP             = cell( 1, nfBins );
            for i           = 1 : nfBins
                verb( sprintf( '%3.1f ', freqs( i ) ), -verbose )
                if i > 1 && i < nfBins
                    bw      = [ diff( freqs( i - 1 : i ) ) diff( freqs( i : i + 1 ) ) ] / 2;
                elseif i == 1
                    bw      = diff( freqs( i : i + 1 ) ) * [ -1 1 ] / 2;
                else
                    bw      = diff( freqs( i - 1 : i ) ) * [ -1 1 ] / 2;
                end
                hBP{ i }    = makefir( freqs( i ) + bw, eegFs, [], 'bandpass' );
            end
            verb( sprintf( 'done.' ), verbose )
        else
            % get central frequencies from filter bank
            freqs           = firFreqs( hBP, eegFs );                                               % this assumes a FIR filter!
            nfBins          = length( hBP );
        end
        fMin                = min( freqs );
        maWin               = ceil( pi ./ freqs * Fs );  
    case 'dog'                                                                                      % these are Gaussian filters!
        if length( hBP ) ~= 3
            error( 'for DOG, three fir filters must be provided!' )
        end
        fir1                = hBP{ 1 };
        fir2                = hBP{ 2 };
        fir3                = hBP{ 3 };
        len                 = zeros( 1, 3 );
        for i               = 1 : 3
            len( i )        = length( hBP{ i } );
        end
        afreqs              = Fs ./ len * 2;
        freqs               = mean( afreqs( 1 : 2 ) );

end
if ismember( specMode, { 'fir', 'dog' } )
    windur                  = max( windur, length( hBP{ 1 } ) / 2 );
else
    windur                  = max( windur, 1 / fMin * Fs );
end

%------------------------------------------------------------%
% determine epochs randomly subject to temporal constraints
%------------------------------------------------------------%
verb( sprintf( '\n%s: Determining normalization epochs..'...
    , upper( mfilename ) ), -verbose )
if isempty( epochs )
    swsfname                = [ filebase '.sts.sws' ];
    if ~exist( swsfname, 'file' )
        verb( 'no epochs specified or stored...', verbose )
        segmentBehavior( filebase, 'windur', behaviorDurSEC, 'Overwrite', OverwriteStates );
    end
    epochs                  = load( swsfname );                                                     % SWS epochs
else
    epochs                  = sortranges( epochs );
end
a                           = memmapfile( sourcefile, 'Format', 'int16' );
nsamples                    = length( a.data ) / nchans;
clear a
epochs                      = intersectranges( [ 1 nsamples ], epochs );
Vals                        = LoadStims( filebase );
if ~isempty( Vals )
    vals                    = resampleranges( Vals( :, 1 : 2 ), Fs, par.SampleRate );
    pad                     = [ floor( padBuffer( 1 ) * Fs ) ceil( padBuffer( 2 ) * Fs ) ];
    vals                    = [ vals( :, 1 ) + pad( 1 ) vals( :, 2 ) + pad( 2 ) ];
    epochs                  = setdiffranges( epochs, vals );                                        % SWS w/o stim
end
durs                        = diff( epochs, [], 2 ) + 1;
pad                         = ceil( windur / 2 );
snipdur                     = windur + 2 * pad;
switch specMode
    case 'spc'
        tidx                = 2;
    otherwise
        tidx                = pad + 1 + ( 0 : windur - 1 );
end
rmv                         = durs < snipdur;
epochs( rmv, : )            = [];                                                                   % short epochs removed
durs( rmv, : )              = [];
if isempty( durs )
    fprintf( '%s: no windows of length %d sec left!! Aborting.\n'...
        , upper( mfilename ), snipdur / Fs )
    return
end 
wins                        = floor( durs / snipdur );
wins                        = cumsum( wins );
nwins                       = wins( end );                                                          % number of independent windows
[ ~, idx ]                  = sort( rand( nwins, 1 ) );                                             % randomize and take many small segments
nrand                       = min( nwins, nrand );
swins                       = sort( idx( 1 : nrand ) );                                             % choose a random set of windows
mat                         = [ [ 1; wins( 1 : end - 1 ) + 1 ] wins ];
rperiods                    = zeros( nrand, 2 );
k                           = 1;
for i                       = 1 : size( mat, 1 )
    slct                    = swins( swins >= mat( i, 1 ) & swins <= mat( i, 2 ) ) - mat( i, 1 ) + 1;
    n                       = length( slct );
    if n > 0
        i0                  = epochs( i, 1 ) : snipdur : epochs( i, 2 );
        rperiods( k : ( k + n - 1 ) ) = i0( slct )';
        k                   = k + n;
    end
end
rperiods( :, 2 )            = rperiods( :, 1 ) + snipdur - 1;
verb( sprintf( '%d/%d/%d %0.3g sec segments selected!'...
    , size( rperiods, 1 ), nrand, nrand0, snipdur / Fs ), verbose )

%------------------------------------------------------------%
% load and compute
%------------------------------------------------------------%
verb( sprintf( '%s: Preprocessing... ', upper( mfilename ) ), -verbose )
if doCSD
    fidx                    = find( ismember( chans, specChans ) );
    fidx                    = unique( clipmat( [ fidx minmax( fidx ) + [ -1 1 ] ], 1, length( chans ) ) );
    lchans                  = chans( fidx );
else
    lchans                  = chans;
end
x                           = readbin( sourcefile, lchans, nchans, rperiods ) * scalefactor;
if ~isempty( refchan )
    xref                    = readbin( sourcefile, refchan, nchans, rperiods ) * scalefactor;
    x                       = bsxfun( @minus, x, mean( xref, 1 ) );
end
if doCSD
    x                       = bsxfun( @minus, x, mean( x, 2 ) );
    if expandEdges
        x                   = -diff( x( [ 1 1 : size( x, 1 ) size( x, 1 ) ], : ), 2 );
    else
        x                   = -diff( x, 2 );
    end
end
if whitenFlag
    x                       = WhitenSignal( x', [], 0 )';
end
xmat                        = reshape( x, [ size( x, 1 ) snipdur nrand ] );
if doCSD
    if expandEdges
        xmat                = squeeze( mean( xmat( ismember( lchans, specChans ), :, : ), 1 ) );
    else
        xmat                = squeeze( mean( xmat( ismember( lchans( 2 : end - 1 ), specChans ), :, : ), 1 ) );
    end
else
    xmat                    = squeeze( mean( xmat( ismember( chans, specChans ), :, : ), 1 ) );
end
if nrand == 1
    xmat                    = xmat'; 
end
xmat                        = bsxfun( @minus, xmat, mean( xmat, 1 ) );

verb( sprintf( 'computing normalization factors for segment ' ), -verbose )
pow                         = [];
for j                       = 1 : nrand
    if verbose == 1 && ~mod( j, ceil( nrand / 10 ) )
        fprintf( '%d ', j )
    end
    xj                      = double( xmat( :, j ) );
    switch specMode
        case 'spc'
            [ pow( :, :, j ), freqs, too ] = mtcsg( xj, nFFT, Fs, nWindow, nWindow/2, mtNW, dflag );
        case 'wlt'
            [ pow( :, :, j ), freqs ] = getWavelet( xj, Fs, fMin, fMax, nfBins - 1 );
        case 'fir'
            for i           = 1 : nfBins
                eegf        = firfilt( xj, hBP{ i } );
                pow( i, :, j ) = ma_rms( eegf, maWin( i ) );
            end
        case 'dog'
            xNotHigh        = firfilt( xj, fir2 );
            xLo             = firfilt( xNotHigh, fir1 );
            xBP             = xNotHigh - xLo;
            xPower0         = firfilt( abs( xBP ), fir3 );
            pow( 1, :, j )  = xPower0;
    end
end
pow                         = pow( :, tidx, : );
pow                         = reshape( pow, [ length( freqs ) length( tidx ) * nrand ] );
verb( sprintf( ' done!' ), verbose )

%------------------------------------------------------------%
% summarize and save
%------------------------------------------------------------%
mm                          = mean( pow, 2 );
ss                          = std( pow, [], 2 );
ssq                         = mean( pow.^2, 2 );
generator                   = { computer, datestr( now, 'ddmmmyy' ) };

if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( savename, 'file' ) )
    verb( sprintf( '%s: Saving %s', upper( mfilename ), savename ), verbose )
    save( savename, 'filebase', 'suffix', 'Fs', 'scalefactor', 'generator'...
        , 'specMode', 'specChans', 'snipdur', 'nrand'...
        , 'mm', 'ss',  'ssq', 'freqs', '-v6' );
end

clear s
s.filebase                  = filebase;
s.suffix                    = suffix;
s.Fs                        = Fs;
s.scalefactor               = scalefactor;
s.generator                 = generator;
s.specMode                  = specMode;
s.specChans                 = specChans;
s.chans                     = chans;
s.doCSD                     = doCSD;
s.snipdur                   = snipdur;
s.nrand                     = nrand;
s.mm                        = mm;
s.ss                        = ss;
s.ssq                       = ssq;
s.freqs                     = freqs;

return

% EOF

% to test:

filebase = datenum2filebase( { 'm660r1', 88 } );
eegchans = 41;
[ mm, ss, freqs, ssq, nrand ] = computeBaseline( mfilebase, eegchans( 1 ), 'wlt', 'nrand', 200, 'Overwrite', -2 );

