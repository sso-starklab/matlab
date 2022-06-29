% spikePhaseFreq            compute spike frequency-phase map
%
% [ h t xbins ybins b ] = spikePhaseFreq( clu, res, eeg )
% 
% DOES
% This is a low-level function: it should receive a set of clu/res vectors
% and an eeg signal (vector), plus the phases and frequencies to evaluate
% (see formatting below) and a power threshold. It then filters the eeg in
% the various frequencies and computes the phases (using either a filter
% bank or CWT), and determines periods of high power in each band.
% It does not save anything, but it optionally plots the results. 
%
% ARGUMENTS
% clu/res       equal size vectors
% eeg           simultaneously recorded eeg signal
%
% OPTIONAL ARGUMENTS
% phases        the phases. {20}. Options: 
%                   -Vector of phase edges covering the range -pi:pi,
%                   thus the first and last must be the same (mod 2pi).
%                   -Scalar indicating the number of phase bins
% freqs         the frequencies. { 2, 50, 25, 'linear' }. Options:
%                   -a scalar, interpreted as the nfft (rounded to the next
%                   power of two). The number of frequency bins will be
%                   nfft/2, spaced linearly between DC and Nyquist
%                   -a vector of discrete frequencies
%                   -a 4-element cell array vector, interpreted as a set: 
%                       [ fmin fmax nbins method ], 
%                   where 'method' is either 'linear', 'log', or 'wavelet'
%                   -a cell array, interpreted as a filter bank
% powTH         {2},. [SDs]; either a scalar (interpreted as SD of the eeg power) or a
%                   vector (then should be same length as freqs)
% spkFs         {20000}, [Hz]
% eegFs         {1250}, [Hz]
% bperiods      {[]}, [eeg samples]; matrix indicating ranges to be used
%                   for baseline determination
% M, SD         {[]}, [mV] Optional mean and SD at each frequency bin to be
%                   used for segment selection. If supplied and valid,
%                   bperiods is ignored. 
% compute       {1} flow flag
% graphics      {0} "
% verbose       {0} "
% 
% RETURNS
% h             count histograms. this is 3D: phase x freq x clunum. the k'th
%                   sheet corresponds to the k'th element of unique( clu )
% t             [sec]; occupancy histogram . for powTH 0, this is a uniform vector
% x/ybins       [rad],[Hz]; the bins centers 
% b             [spikes/sec]; the baseline rates, to be used for computing gain maps.
% 
% ADDITIONAL OUTPUTS
% hRate         rate histograms. same dimensions as h and easily
%                   reconstructed from the other output as follows:
%                   hRate( i, j, k ) = h( i, j, k ) / nCycles( j ) / binSize( j );
%                   -binSize( j ) is the duration [sec] of the a phase bin
%                   for the j'th frequency, identical to 1 / ybins( j ) / length( xbins )
%                   -nCycles( j ) is the number of cycles of the
%                   j'th frequency bin, identical to durs( j ) * ybins( j )
%                   for powTH 0, this corresponds to sum( sum( h( :, :, i ) ) )
% M, SD         the actual mean and SD used for segmentation
% ranges        the time ranges used for each frequency bin
% hBP           the filter bank
% eegPhs        the phases (huge array.. length( eeg ) x length( freqs )
% 
% USAGE:
% Since this is a low-level function, the calling routine should take care
% of data handling (blockwise loading), pre-processing (csd, averaging over
% multiple channels...), segmentation (brain states, animal position,
% stimulation times...), and accumulation over segments
% 
% 
% calls         ParseArgPairs, colvec                                                           (general)
%               dilutesegments, getdatainranges, inranges, parse, resampleranges, sortranges    (sets) 
%               imupsample, alines, myjet                                                       (graph)
%               makefir, ma_rms, local_max, getWavelet, firfilt, plotSpectrogram                (ssp)
%
% see also      detect_hfos, wnAnalysis, spikePhase, spikeWB

% 11-mar-13 ES

% revisions
% 14-mar-13 adapated for use with spikePhase (the wrapper)
% 18-mar-13 added support for external eegPhs/ranges (for resampling purposes, see spikeWB)
% 17-aug-19 added option to extact freqs from a list that accompanies the filter bank
% 18-aug-19 cleaned up

function [ h, durs, phsBins, freqs, b, hRate, M, SD, ranges, hBP, eegPhs ] = spikePhaseFreq( clu, res, eeg, varargin )

h                           = [];
durs                        = [];
b                           = [];
hRate                       = [];
ranges                      = [];

% handle the arguments...
nargs                       = nargin;
if nargs < 3 || isempty( clu ) || isempty( res ) || isempty( eeg )
    error( 'all three arguments are compulsory' )
end
eeg                         = colvec( eeg );
neeg = length( eeg );
if size( eeg, 2 ) ~= 1
    error( 'single channel supported only' )
end
if size( res, 2 ) ~= 1 || size( clu, 2 ) ~= 1
    error( 'clu/res should be column vectors' )
end
if size( res, 1 ) ~= size( clu, 1 )
    error( 'clu/res mismatch' )
end

[ phases, freqs, powTH, bperiods, M, SD, spkFs, eegFs, compute, graphics, verbose ] = ParseArgPairs(...
    { 'phases', 'freqs', 'powTH', 'bperiods', 'M', 'SD', 'spkFs', 'eegFs', 'compute', 'graphics', 'verbose' }...
    , { 20, { 2, 50, 25, 'linear' }, 2, [], [], [], 20000, 1250, 1, 0, 0 }...
    , varargin{ : } );
if isempty( bperiods )
    bperiods                = [ 1 neeg ];
end
bperiods                    = sortranges( bperiods );
if spkFs( 1 ) < 0 || length( spkFs ) > 1
    error( 'spkFs should be a sampling rate' )
end
if eegFs( 1 ) < 0 || length( eegFs ) > 1
    error( 'eegFs should be a sampling rate' )
end
if eegFs > spkFs
    error( 'eegFs should not be larger than spkFs' )
end

%------------------------------------------------------------%
% preparations
%------------------------------------------------------------%

%------------------------------------------------------------%
% downsample the spike trains
if eegFs                    ~= spkFs
    res                     = resampleranges( res, eegFs, spkFs );
end

%------------------------------------------------------------%
% determine the frequency bins/filters/method
err                         = 0;
if isnumeric( freqs ) && length( freqs ) == 1 && freqs > 0
    % scalar: assume nfft and generate a linear array of frequencies
    nfft                    = 2 ^ ceil( log2( freqs ) );
    freqs                   = ( 1 : nfft / 2 )' * eegFs / nfft;
    method                  = 'linear';
elseif ~isnumeric( freqs ) && iscell( freqs ) && length( freqs ) == 4
    % 4-element specification array: determine the filter-bank parameters
    fMin                    = freqs{ 1 };
    fMax                    = freqs{ 2 };
    nfBins                  = freqs{ 3 };
    method                  = freqs{ 4 };    
    switch method
        case 'linear'
            if nfBins == 1
                freqs       = ( fMin + fMax ) / 2;
            else
                freqs       = linspace( fMin, fMax, nfBins );
            end
        case 'log'
            freqs           = logspace( log10( fMin ), log10( fMax ), nfBins );
        case 'wavelet'
            if fMin > fMax || fMax > eegFs / 2
                error( 'parameter mismatch' )
            end
        otherwise
            err             = 1;
    end
elseif iscell( freqs ) && isvector( freqs ) && length( freqs ) == 3 ...
    && isa( freqs{ 1 }, 'cell' ) && isa( freqs{ 2 }, 'numeric' ) ...
    && isequal( size( freqs{ 1 }, 1 ), size( freqs{ 2 }, 2 ), numel( freqs{ 3 } ) )
    % assume already processed data: { Ranges, eegPhs, frequencies }
    method                  = 'preprocessed';
    Ranges                  = freqs{ 1 };
    eegPhs                  = freqs{ 2 };
    freqs                   = freqs{ 3 };
    nfBins                  = length( freqs );
    if size( eegPhs, 1 ) ~= neeg
        err                 = 1;
    end
elseif iscell( freqs ) && isvector( freqs ) && isa( freqs{ 1 }, 'numeric' )
    % cell array: asssume a filter bank, determine the filter frequencies
    hBP                     = freqs;
    nfBins                  = length( hBP );
    freqs                   = firFreqs( freqs, eegFs ); % note - imprecise for high frequencies
elseif iscell( freqs ) && length( freqs ) == 2 && isvector( freqs{ 1 } ) && isa( freqs{ 1 }{ 1 }, 'numeric' ) && isvector( freqs{ 2 } )
    % two-element cell array: assume a filter bank and list of frequencies
    hBP                     = freqs{ 1 };
    nfBins                  = length( hBP );
    freqs                   = freqs{ 2 };
elseif isnumeric( freqs ) && isvector( freqs )
else
    err                     = 1;
end
if err
    error( 'unrecognized format for freqs' )
end
freqs                       = freqs( : ).';
if ~exist( 'method', 'var'  ) && isnumeric( freqs ) && isvector( freqs )
    % given array of frequencies: check if linear or log spaced
    logfreq                 = logspace( log10( freqs( 1 ) ), log10( freqs( end ) ), length( freqs ) );
    linfreq                 = linspace( freqs( 1 ), freqs( end ), length( freqs ) );
    if isequal( freqs, linfreq )
        method              = 'linear';
    elseif isequal( freqs, logfreq )
        method              = 'log';
    else
        freqs               = sort( freqs );
        if length( freqs ) > 1
            fprintf( '%s: Assuming log-spacing!!!\n', upper( mfilename ) )
            method          = 'log';
        else
            method          = 'linear';
        end
    end
end

%------------------------------------------------------------%
% if needed, generate a filter bank
if ismember( method, { 'linear', 'log' } )
    nfBins                  = length( freqs );
    maWin                   = ceil( pi ./ freqs * eegFs );       
    if ~exist( 'hBP', 'var' )
        verb( sprintf( '%s: Constructing filter bank (%s spacing)...'...
            , upper( mfilename ), method ), -verbose )
        switch method
            case 'log'
                eF          = freqs( 1 ) / 2;
                eL          = 2 * freqs( end ) - geomean( freqs( end - 1 : end ) );
                edges       = [ eF geomean( [ freqs( 1 : end - 1 ); freqs( 2 : end ) ], 1 ) eL ];
            case 'linear'
                if length( freqs ) == 1
                    df      = ( fMax - fMin ) / 2;
                else
                    df      = diff( freqs( 1 : 2 ) ) / 2;
                end
                if fMin < df
                    e0      = fMin / 2; 
                else 
                    e0      = fMin - df; 
                end
                edges       = [ e0 linspace( fMin + df, fMax + df, nfBins ) ];
        end
        hBP                 = cell( 1, nfBins );
        for i               = 1 : nfBins
            verb( sprintf( '%3.1f ', freqs( i ) ), -verbose )
            band            = edges( i : i + 1 );
            hBP{ i }        = makefir( band, eegFs, [], 'bandpass' );
        end
        verb( sprintf( 'done.' ), verbose )
    end
end

%------------------------------------------------------------%
% determine the phase bins
if length( phases ) == 1
    npBins                  = abs( round( phases ) );
    phases                  = ( -pi : 2 * pi / npBins : pi )';
end
if length( unique( mod( phases( [ 1 end ] ), 2 * pi ) ) ) ~= 1
    error( 'Phase edges should cover the entire unit circle' )
end
phases                      = phases( : );
npBins                      = length( phases ) - 1;
phsBins                     = ( phases( 1 : npBins ) + phases( 2 : npBins + 1 ) ) / 2;
if ismember( method, { 'linear', 'log', 'preprocessed' } )
    phsBinSize              = 1 / npBins ./ freqs; % bin duration in sec
end

%------------------------------------------------------------%
% determine the power threshold
if length( powTH ) == 1
    powTH                   = ones( nfBins, 1 ) * powTH;
elseif length( powTH ) ~= nfBins
    error( 'powTH should be either a scalar or a vector of same length as freqs (%d)', nfBins )
end

%------------------------------------------------------------%
% determine validity of external mean, SD
localStats                  = 1;
if ~isempty( M ) && ~isempty( SD )
    if length( M ) == 1 
        M                   = M * ones( nfBins, 1 );
    end
    if length( SD ) == 1 
        SD                  = SD * ones( nfBins, 1 );
    end
    if length( M ) == nfBins && length( SD ) == nfBins
        localStats          = 0;
    end
end
if localStats
    M                       = zeros( nfBins, 1 );
    SD                      = zeros( nfBins, 1 );
end
if ~compute
    return
end

%------------------------------------------------------------%
% phase and period computation
%------------------------------------------------------------%
ranges                      = cell( nfBins, 1 );
durs                        = zeros( nfBins, 1 );
switch method
    case 'preprocessed'
        % don't do anything
        ranges              = Ranges;
    case { 'linear', 'log' }
        % compute the band-specific phase using simple filtering (alternatively do wavelet)
        verb( sprintf( '%s: Computing eeg phases via filter bank...'...
            , upper( mfilename )  ), -verbose )
        eegPhs              = single( zeros( size( eeg, 1 ), nfBins ) );
        for i               = 1 : nfBins
            verb( sprintf( '%3.1f ', freqs( i ) ), -verbose )
            if neeg < maWin( i )
                continue
            end
            eegf            = firfilt( eeg, hBP{ i } );
            eegPhs( :, i )  = single( angle( hilbert( eegf ) ) );
            if powTH( i ) ~= 0
                eegPow      = ma_rms( eegf, maWin( i ) );
                if localStats 
                    eegPowBL = getdatainranges( eegPow, bperiods );
                    M( i )  = mean( eegPowBL );
                    SD( i ) = std( eegPowBL );
                end
                TH          = M( i ) + powTH( i )* SD( i );
                mat         = parse( find( eegPow > TH ) );
                if isempty( mat )
                    continue
                end
                mat         = dilutesegments( mat, maWin( i ), maWin( i ) );
            else
                mat         = [ 1 neeg ];
            end
            durs( i )       = sum( diff( mat, [], 2 ) + 1 ) / eegFs;
            ranges{ i }     = mat;
        end
    case 'wavelet'
        hBP                 = []; % returned empty
        verb( sprintf( '%s: Computing eeg phases via CWT...'...
            , upper( mfilename )  ), -verbose )
        [ eegPow, freqs, ~, ~, eegPhs ] = getWavelet( eeg...
            , eegFs, fMin, fMax, nfBins - 1 );
        verb( sprintf( 'done  CWT...' ), -verbose )
        eegPow              = single( eegPow' );
        eegPhs              = single( eegPhs' );
        nfBins              = length( freqs );
        phsBinSize          = 1 / npBins ./ freqs;
        maWin               = ceil( pi ./ freqs * eegFs );
        for i               = 1 : nfBins
            if powTH( i ) == 0
                mat         = [ 1 neeg ];
            else
                if localStats
                    eegPowBL = getdatainranges( eegPow( :, i ), bperiods );
                    M( i )  = mean( eegPowBL );
                    SD( i ) = std( eegPowBL );
                end
                TH          = M( i ) + powTH( i )* SD( i );
                mat         = parse( find( eegPow( :, i ) > TH ) );
                if isempty( mat )
                    continue
                end
                mat         = dilutesegments( mat, maWin( i ), maWin( i ) );
            end
            durs( i )       = sum( diff( mat, [], 2 ) + 1 ) / eegFs;
            ranges{ i }     = mat;
        end
end
verb( sprintf( 'determined all phases and ranges!' ), verbose )
clear eegPow eegPowBL

%------------------------------------------------------------%
% generate the histograms
%------------------------------------------------------------%
uclu                        = unique( clu );
nclu                        = length( uclu );
h                           = single( zeros( npBins, nfBins, nclu ) );
hRate                       = h;
b                           = zeros( nclu, 1 );
verb( sprintf( '%s: Computing for unit#: '...
    , upper( mfilename ) ), -verbose )

for k = 1 : nclu
    clunum                  = uclu( k );
    idx                     = clu == clunum;
    resk                    = res( idx );
    b( k )                  = length( resk ) / neeg * eegFs;
    verb( sprintf( '%d ', uclu( k ) ), -verbose )
    for i = 1 : nfBins
        idx                 = inranges( resk, ranges{ i } );
        spkPhs              = eegPhs( resk( idx ), i );
        if phases( 1 ) < -pi
            phsHist         = histc( spkPhs, [ phases; pi ] );
            phsHist( 1 )    = phsHist( 1 ) + phsHist( npBins + 1 );
        else
            phsHist         = histc( spkPhs, phases );
        end
        phsHist             = reshape( phsHist( ( 1 : npBins ) ), [ npBins 1 ] );
        h( :, i, k )        = single( phsHist );
        totcycles           = durs( i ) * freqs( i );
        hRate( :, i, k )    = phsHist / phsBinSize( i ) / totcycles;
    end
end
verb( sprintf( ' determined all spike counts (%d spikes)!', length( clu ) ), verbose )

%------------------------------------------------------------%
% graphics
%------------------------------------------------------------%
if ~graphics
    return
end

[ ~, ah ]                   = plotSpectrogram( phsBins, freqs, hRate, 5, 'imagesc', 'linear', 'linear', 2 );
for k  = 1 : nclu
    subplot( ah( k ) )
    axis square
    title( sprintf( '%d: %0.3g', uclu( k ), max( max( hRate( :, :, k ) ) ) ) )
end

return

% EOF
