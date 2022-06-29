% streconstruct         reconstruct analog waveform from spike trains
% 
% call                  [ yhat, stats, xc, ax ] = streconstruct( x, y )
%
% gets                  x               spike train/s. sparse matrix (m x n) or a vector of times
%                       y               analog signal. m x 1 (or m x n, or p x 1)
%
% returns               yhat            reconstructed analog signal, same dimensions as y
%                       f               reconstruction filters (mf x nf matrix)
%                       xc              classified spike train/s (same dimension as x)
%                       ax              axes handles (on a new figure)
%
% optional arguments (given as name/value pairs)
%   
%                       minSpikes             {5}; the minimal total number of spikes in x
%                       minTrials             {4}; the minimal number of trials
%                       minSpikesPerTrial     {1}; the minimal mean number of spikes/trial
%                       prepy                 {'zs'} (standradize) or 'mean' (remove mean)
%                       prepMode              {'flat'} (ignore spikes at edges) or 'pad' (zero-pad; triangular-bias in STA/auto-corr)
%                       ftype                 {'sta'}, 'wiener' (first order)
%                                             'csta' (2nd order model), 'cwiener', 'xsta' (requires x2)
%                       x2                    {[]}; a reference set of spike trains (same arrangement as x). relevant for 'xsta' only
%                       nfilt                 {50}; [samples]; filter range. a scalar indicates [ -1 1 ] * nfilt
%                       maxISIbins            {10}; maximal number of ISI bins (relevant for 'csta' and 'xsta' only)
%                       minCountPerBin        {10}; minimal number of spikes/ISI bin (relevant for 'csta' and 'xsta' only)
%                       isiedges              {[]}; [samples]; externaly-specified set of ISI bin edges (overrides maxISIbins/minCountPerBin)
%                       lagsr                 {0}; [samples]; range of lags to optimize over. a scalar
%                                                 indicates a single lag, and positive lags indicate causal filter
%                       cctype                {'rank'} or 'pmoment'; CC type used for post-hoc quantification 
%                       nreps                 {10} number of randomizations per trial
%                       rhalfwin              {[]} (shuffle) or a number (jitter half win)
%                       graphics              {1} or 0; to plot or not (one summary figure)
%                       Fs                    {1}; [samples/s]; sampling rate - used ONLY for ftim/graphics (i.e. if is kHz, abscissa in ms etc)
% 
% Notes:
% 1. This is a low level routine, so everything is in samples; x and y must be of the same sampling rate. 
%       If this is not the case, the calling routine should downsample/upsample appropriately
% 2. Only relevant for zero-mean signals (waveform reconstruction)
% 
% calls                 colvec, make_equal_bins, makeblocks, replacetok, scaleto, shift, ParseArgPairs  (general)
%                       alines, calibration, myjet, patch_band, textf                                   (graph) 
%                       sta, stacalc, stahat                                                            (lfp)
%                       isisort, plot_raster, stmix, stwiener                                           (spikes)
%                       calc_fwhh, firfilt, mtcsd1, my_xcorr                                            (ssp)
%                       calc_pearson, calc_spearman, inrange, make_ai, rankcols, zsc                    (stats)
%
% references            Bialek et al, 1991, Science
%                       Rieke et al., 1997, Spikes
%                       van Rossum 2001, Neural Computation
%                       Schreiber et al., 2003, Neurocomputing
%
% see also              streconstructXval

% 09-jun-14 ES

% revisions:
% 10-jun-14 optimization of filter lag for sta
% 15-jun-14 1. optimize lag for csta/xsta
%           2. compute xsta filters
%           3. use external routine calc_isis
%           4. allow external specification of isiedges, lagsr
% 16-jun-14 1. p-value for reconstruction by shuffling/jittering (single trial + global)
%           2. global measures added
%           3. coherence and phase computed
%           4. temporal statistics (lag, width, ACH) added
% 19-jun-14 1. in case of 'csta' with too few ISIs, tries to do an STA/returns
% 20-jun-14 1. allowed effective nibins differnt from requested (can happen with discrete ISIs)
% 29-sep-14 1. minor modifications
%           2. better annotations
% 30-sep-14 1. added widCC2 (FWHH of pairwise CC between reconstructions)
%               and reliability + randomized (van Rossum) [optional]
% 27-oct-14 1. modularize: isisort.m, stacalc.m, stahat.m
%           2. corrected error in doRel (nansum should have been nanmean)
%           3. ax returned
% 28-oct-14 1. minTrials (3) added
%           2. reliability p-value computed
% 29-oct-14 1. minor cleanup, modifications of graphics
%           2. modified default minCountPerBin from 10->5
% 02-nov-14 1. bug in filter bank computation corrected (x->xpad; xt->xc)
% 11-jan-15 1. wiener filtering added
% 21-jun-15 1. wiener filtering generalized
%           2. pre-processing modified to prevent bias in STA/auto-corr
%           3. added support for xwiener, cwiener, wiener
% 14-oct-19 1. cleaned up and documented

% streconstruct( iST( 250 : end, : ), wn( 250 : end ), 'nfilt', 125, 'ftype', 'sta', 'Fs', 2.5 );

% Open points (not critical): 
% 1. determine why [ 0 50 ] gets a different optimum compared to a [ -50 50 ]  lag range
% 2. determine reason for slight difference with cSTA filters (and reconstrction) upon internal/external bedges
% 3. quantify information rate (bits/s, bits/spike)
% 4. PSTH reconstruction?? (can do also for cSTA - according to the tagging)

function [ yhat, stats, xc, ax ] = streconstruct( x, y, varargin )

yhat                    = [];
stats                   = [];
xc                      = [];
ax                      = [];

%----------------------------------------------------------------------
% constants
%----------------------------------------------------------------------
% reliability
BLOCKSIZE               = 100;

% coherence
fROI                    = [ 0 200 ];
M                       = 0.25;
mtNW                    = 3;
dflag                   = '';

% graphics
colors                  = [ 0 0 0.7; 0 0 0; 1 0 0 ]; % x, y, yhat
blackColor              = [ 0 0 0 ];
whiteColor              = [ 1 1 1 ];
grayColor               = [ 1 1 1 ] * 0.6;
purpleColor             = [ 1 0.5 1 ]; 
redColor                = [ 1 0 0 ];
tstr                    = '';

%----------------------------------------------------------------------
% arguments
%----------------------------------------------------------------------
nargs                   = nargin;
if nargs < 2 || isempty( x ) || isempty( y )
    return
end
[ minSpikes, minTrials, minSpikesPerTrial, prepy ...
    , prepMode ...
    , ftype, x2, nfilt ...
    , maxISIbins, minCountPerBin, isiedges ...
    , lagsr ...
    , cctype, nreps, rhalfwin, doRel ...
    , doScatter, doCoherence, Fs, graphics, verbose ] = ParseArgPairs(...
    { 'minSpikes', 'minTrials', 'minSpikesPerTrial', 'prepy' ...
    , 'prepMode' ...
    , 'ftype', 'x2', 'nfilt' ...
    , 'maxISIbins', 'minCountPerBin', 'isiedges' ...
    , 'lagsr' ...
    , 'cctype', 'nreps', 'rhalfwin', 'doRel' ...
    , 'doScatter', 'doCoherence', 'Fs', 'graphics', 'verbose' }...
    , { 5, 4, 1, 'zs' ...
    , 'flat' ...
    , 'sta', [], 50 ...
    , 10, 5, [] ...
    , 0 ...
    , 'rank', 10, [], 0 ...
    , 0, 0, 1, 0, 1 }...
    , varargin{ : } );

% scale analog
y                       = colvec( y );
y0                      = y;
prepy                   = lower( prepy );
switch prepy
    case 'zs'
        y               = zsc( y0, 1 );
        zstr            = 'Z';
    case 'mean'
        y               = bsxfun( @minus, y, nanmean( y, 1 ) );
        zstr            = 'zraw';
    otherwise
        error( 'written ONLY for zero-mean signals' )
end

% input size
y1                      = y;
if issparse( x )
    [ m, nt ]           = size( x );
    [ my, ny ]          = size( y );
    if nt > 1
        if ny == 1
            if m == my
                y       = repmat( y, [ 1 nt ] );
            else
                error( 'input size mismatch: X sparse, Y mismatching vector!' )
            end
        elseif nt ~= ny
            error( 'input size mismatch: X sparse, Y matrix with mismatching columns!' )
        elseif m ~= my
            error( 'input size mismatch: X sparse, Y matrix with mismatching rows!' )
        end
    elseif m ~= my
        error( 'input size mismatch: X sparse, Y matrix with mismatching rows!' )
    end
else
    error( 'not supported yet: X must be sparse' )
end
if nt < minTrials
    if verbose
        fprintf( '%s: Too few trials!! (%d)\n', upper( mfilename ), nt )
    end
    return
end

% filter type
ftype                   = lower( ftype );
if ~ismember( ftype, { 'wiener', 'cwiener', 'xwiener', 'sta', 'csta', 'xsta', 'corr', 'sta1', 'sta2', 'sta3' } )
    error( 'mismatching filter type' )
end
switch ftype
    case { 'sta', 'sta1' }
        ftype           = 'sta';
    case { 'csta', 'sta2' }
        ftype           = 'csta';
    case { 'xsta', 'sta3' }
        ftype           = 'xsta';
end
if isequal( ftype, 'xsta' ) && ~isequal( [ m nt ], size( x2 ) )
    error( 'input size mismatch: X and X2 must have identical dimensions' )
end

% filter length
if length( nfilt ) == 1
    nfilt               = nfilt * [ -1 1 ];
end
if length( nfilt ) ~= 2
    error( 'NFILT must be a range' )
end
nfilt                   = round( sort( nfilt ) );
maxlag                  = max( abs( nfilt ) );

% isis bins
if ismember( ftype, { 'cwiener', 'csta', 'xsta' } ) 
    if isempty( isiedges )
        if ( isempty( maxISIbins ) || maxISIbins < 1 )
            error( 'MAXISIBINS must be non-negative' )
        end
        maxISIbins      = round( maxISIbins );
    else
        [ ~, tmp ]      = size( isiedges );
        if tmp ~= 1
            error( 'ISIBINS must be a vector of edges' )
        end
        isiedges        = sort( isiedges );
    end
end

% randomization scale
rhalfwin                = abs( round( rhalfwin ) );

%----------------------------------------------------------------------
% prepare for filter computation
%----------------------------------------------------------------------
% pad to prevent edge effects
switch prepMode
    case 'pad'
        padbuffer       = zeros( maxlag, nt );
        ypad            = [ y; padbuffer ];                                 % does not remove any spikes
        xpad            = [ x; padbuffer ];                                 % does not modify the analog
    case { 'flat', 'unbiased' }
        x( 1 : maxlag, : ) = 0;
        x( end - maxlag + 1 : end, : ) = 0;
        xpad            = x;                                                % removed data from the spikes
        ypad            = y;                                                % no padding of the analog
        padbuffer       = [];                                               % to support c-sta/c-wiener etc
end

% initialize
nspikes                 = full( sum( x( : ) ) );                            % number of spikes
nf                      = diff( nfilt ) + 1;
f                       = NaN * ones( nf, 1 );
fsem                    = f;
ftim                    = ( nfilt( 1 ) : nfilt( 2 ) )' / Fs;                % [s]
bedges                  = [];
cpb                     = [];
sidx                    = sparse( nspikes, 1 );

% make sure enough spikes:
if nspikes < minSpikes || nspikes / nt < minSpikesPerTrial
    if verbose
        fprintf( '%s: Too few spikes!! (%d trials, %d spikes)\n', upper( mfilename ), nt, nspikes )
    end
    return
end

%----------------------------------------------------------------------
% compute filters
%----------------------------------------------------------------------

xc                      = x;
switch ftype
    
    case 'corr'                                                             % reverse correlation filter
        f               = my_xcorr( ypad( : ), full( xpad( : ) ), maxlag, 0 );
        f               = f / nspikes;
        fsem            = ones( size( fsem ) ) / sqrt( nspikes );
        
    case 'sta'                                                              % STA filter (mathematically identical, just faster for sparse trains)
        xt              = find( xpad( : ) ); 
        [ f, ss, nn ]   = sta( xt( : ), ypad( : ), nfilt );
        fsem            = ss / sqrt( nn );

    case 'wiener'                                                           % WIENER filter (point process implementation)
        xt              = find( xpad( : ) );
        [ f, fsem ]     = stwiener( xt, ypad( : ), nfilt );
        
    case { 'csta', 'xsta', 'cwiener', 'xwiener' }                           % conditional STA/Wiener filter (depending on the preceding ISI)
        
        % sort the spikes by the ISIs:
        nspikesEff      = sum( max( full( sum( x, 1 ) ), 1 ) - 1 );
        if isempty( isiedges )
            nibins      = min( maxISIbins, floor( nspikesEff / minCountPerBin ) );
            if nibins < 2 
                if verbose
                    fprintf( '%s: Too few ISIs!! (%d trials, %d valid ISIs)\n', upper( mfilename ), nt, nspikesEff )
                end
                return
            end
            isiedges    = nibins;
        end
        switch ftype
            case { 'csta', 'cwiener' }
                [ xc, bedges, cpb ] = isisort( x, isiedges );
            case { 'xsta', 'xwiener' }
                [ xc, bedges, cpb ] = isisort( x, isiedges, x2 );
        end
        nibins          = length( bedges ) - 1;
        sidx            = full( xc( xc > 0 ) );
        
        % compute filter bank:
        xt              = find( [ xc; padbuffer ] );
        switch ftype
            case { 'csta', 'xsta' }
                cmethod = 'sta';
            case { 'cwiener', 'xwiener' }
                cmethod = 'wiener';
        end
        [ f, fsem ]     = stacalc( xt, ypad( : ), nfilt, sidx, 1 : nibins, cmethod );
        
end

% compute the 1st order STA filter anyhow
if ismember( ftype, { 'corr', 'sta', 'wiener' } )
    f1                  = f;
    f1sem               = fsem;
else
    [ f1, ss, nn ]      = sta( xt( : ), ypad( : ), nfilt );
    f1sem               = ss / sqrt( nn );
end

%----------------------------------------------------------------------
% reconstruction
%----------------------------------------------------------------------

% determine the lag
if length( lagsr ) == 1                                                     % fixed lag
    optlag              = lagsr;
    f                   = shift( f, -optlag, 0 );
    if optlag > 0
        f( ftim > 0, : ) = 0;
    end        
else
    % optimize over all possible lags:
    fprintf( '%s: Optimizing STA lag: ', upper( mfilename ) )
    lags                = lagsr( 1 ) : lagsr( 2 );
    nlags               = length( lags );
    fx                  = ftim .* Fs;
    fidx                = fx <= 0;
    cci                 = zeros( 2 * maxlag + 1, nlags );
    for i               = 1 : nlags
        fprintf( '.' )
        lag             = lags( i );
        fhat            = shift( f, -lag, 0 );                              % shift
        fhat( ~fidx, : )= 0;                                                % make causal
        yhat            = stahat( xc, fhat, size( y, 1 ), 0 );
        cci( :, i )     = mean( my_xcorr( y, yhat, maxlag, -1 ), 2 );
    end
    % keep the causal filter with the best performance:
    [ ~, lagidx ]       = find( cci == max( cci( : ) ) );
    optlag              = lags( lagidx );
    f                   = shift( f, -optlag, 0 );
    f( ~fidx, : )       = 0;
    fprintf( 'Optimal lag: %d\n', optlag )
end

% reconstruct + shift back
yhat                    = stahat( xc, f, size( y, 1 ), optlag );

%----------------------------------------------------------------------
% quantification + statistics
%----------------------------------------------------------------------

% initialize
pvals                   = NaN * ones( 1, nt );
pvalAll                 = [];
rAll                    = [];

% reconstruction quality (Q := ccAll)
switch cctype
    case 'rank'
        yr              = rankcols( y );
        cct             = my_xcorr( yr, rankcols( yhat ), maxlag, -1 );
        act             = my_xcorr( yr( :, 1 ), yr( :, 1 ), maxlag, -1 );
        cctAll          = my_xcorr( yr( : ), rankcols( yhat( : ) ), maxlag, -1 );
        ccstr           = 'Rank CC';
    case 'pmoment'
        cct             = my_xcorr( y, yhat, maxlag, -1 );
        act             = my_xcorr( y( :, 1 ), y( :, 1 ), maxlag, -1 );
        cctAll          = my_xcorr( y( : ), yhat( : ), maxlag, -1 );
        ccstr           = 'CC';
end
cc                      = cct( maxlag + 1, : );                             % zero lag cc
ccAll                   = nanmean( cc );                                    % mean over trials
[ ccmax, ccmaxlag ]     = max( cct, [], 1 );                                % max cc - may be noisy
[ ccmaxAll, ccmaxAllLag ] = max( cctAll, [], 1 );

% randomize the zero-lag cc (Q) for p-value
if nreps > 0
    fprintf( '%s: Randomizing (%d trials, %d times) ', upper( mfilename ), nt, nreps )
end
yhatRand                = zeros( size( yhat, 1 ), nt );
for i                   = 1 : nt
    if nreps > 0
        fprintf( '.' )
    else
        continue
    end
    xm                  = stmix( x( :, i ), rhalfwin, nreps );
    switch ftype
        case { 'corr', 'sta', 'wiener' }                                    % 1st order
            yhatR       = firfilt( full( xm ), f );
        case { 'csta', 'xsta', 'cwiener', 'xwiener' }                       % 2nd order
            xm          = isisort( xm, bedges, x2 );
            yhatR       = stahat( xm, f, size( y, 1 ), optlag );
    end
    yhatRand( :, i )	= yhatR( :, 1 );
    switch cctype
        case 'rank'
            r           = calc_pearson( yr( :, i ) * ones( 1, nreps ), rankcols( yhatR ) );
        case 'pmoment'
            r           = calc_pearson( y( :, i ) * ones( 1, nreps ), yhatR );
    end
    rAll                = [ rAll r ];
    if nreps < 100                                                          % assume Gaussian distribution of the shuffled data
        mm              = nanmean( r );
        ss              = nanstd( r );
        pvals( i )      = 1 - normcdf( cc( i ), mm, ss );
    else                                                                    % compute empirical p-value
        pvals( i )      = ( sum( r >= cc( i ) ) + 1 ) / ( nreps + 1 );
    end
end
if nreps > 0 && nreps < 100
    mm                  = nanmean( rAll );
    ss                  = nanstd( rAll );
    pvalAll             = 1 - normcdf( ccAll, mm, ss );
elseif nreps > 100
    pvalAll             = ( sum( rAll >= ccAll ) + 1 ) / ( nreps + 1 );
end
if nreps > 0
    fprintf( '\n' )
end
if isempty( rAll )
    rAll                = NaN; 
    pvalAll             = NaN;
end 

% temporal statistics (T := lagSta)
[ widSta, lagSta ]      = calc_fwhh( f1 );                                  % from the 1st order filter
lagSta                  = lagSta - maxlag - 1;
wid                     = calc_fwhh( cct );                                 % from the reconstruction (one per trial)
widAll                  = calc_fwhh( cctAll );                              % from the reconstruction (concatenated data)
widAct                  = calc_fwhh( act );                                 % upper bound: the FWHH of the input auto-correlation

% reliability + width by pair-wise CC between reconstructions (R := ccRel); Schreiber et al., 2003
widCC2                  = NaN;
ccRel                   = [ NaN NaN NaN ];
ccRelRand               = [ NaN NaN NaN ];
if doRel
    
    aidx                = make_ai( nt );
    blocks              = makeblocks( size( aidx, 1 ), BLOCKSIZE, 0 );
    nblocks             = size( blocks, 1 );
    cc2                 = zeros( 2 * maxlag + 1, 1 );
    for i               = 1 : nblocks
        bidx            = blocks( i, 1 ) : blocks( i, 2 );
        cc2i            = my_xcorr( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ), maxlag, -1 );
        cc2             = cc2 + nanmean( cc2i, 2 );
    end
    cc2                 = cc2 / nblocks;
    widCC2              = calc_fwhh( mean( cc2, 2 ) );                      % width of the reliability cross-correlation
    rel                 = cc2( maxlag + 1, : );
    relPval             = signrank( rel, 0 );
    if nanmedian( rel ) > 0
        relPval         = relPval / 2;
    else
        relPval         = 1 - relPval;
    end
    ccRel               = [ nanmean( rel ) calc_sem( rel ) relPval ];
    
    if nreps > 0
        cc2rand         = zeros( 1 );
        for i           = 1 : nblocks
            bidx        = blocks( i, 1 ) : blocks( i, 2 );
            switch cctype
                case 'rank'
                    cc2i = calc_spearman( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ) );
                case 'pmoment'
                    cc2i = calc_pearson( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ) );
            end
            cc2rand     = cc2rand + nansum( cc2i, 2 );
        end
        cc2rand         = cc2rand / nblocks;
        relRand         = cc2rand;
        relRandPval     = signrank( relRand, 0 );
        if nanmedian( relRand ) > 0
            relRandPval = relRandPval / 2;
        else
            relRandPval = 1 - relRandPval;
        end
        ccRelRand       = [ nanmean( relRand ) calc_sem( relRand ) relRandPval ];
    end
    
end

% rescale (for scatter, for plotting)
if isequal( prepy, 'zs' )
    yhat                = zsc( yhat, 1 );
end

% scatter - plot <y>|yhat in [a,b) vs. yhat in [a,b); Rieke et al., 1997
if doScatter || graphics
    nybins              = min( floor( m / 20 ), 50 );
    [ ~, ~, yidx ]      = make_equal_bins( yhat( : ), nybins );
    ysct                = zeros( nybins, 2 );
    for i               = 1 : nybins
        idx             = yidx == i;
        ysct( i, : )    = [ mean( yhat( idx ) ) mean( y( idx ) ) ];
    end
    switch cctype
        case 'rank'
            ccLinear    = calc_spearman( ysct );
        case 'pmoment'
            ccLinear    = calc_pearson( ysct );
    end
else
    ccLinear            = NaN;
end

% coherhence and phase between y and yhat
if doCoherence || graphics
    nFFT                = 2^floor( log2( M * Fs * 1000 ) ); 
    nWindow             = nFFT;
    [ yo, frq ]         = mtcsd1( [ y( :, 1 ) yhat ], nFFT, Fs * 1000, nWindow, nWindow/2, mtNW, dflag );
    frqidx              = inrange( frq, fROI );
    yo                  = yo( frqidx, :, : );
    frq                 = frq( frqidx );
    p1                  = repmat( yo( :, 1, 1 ), [ 1 nt ] );
    p2                  = squeeze( yo( :, 1, 2 : ( nt + 1 ) ) );
    c12                 = squeeze( yo( :, 2, 2 : ( nt + 1 ) ) );
    cohs                = single( abs( c12 .^ 2 ) ./ ( p1 .* p2 ) );
    phs                 = single( atan2( imag( c12 ), real( c12 ) ) );
else
    frq                 = [];
    cohs                = [];
    phs                 = [];
end

%----------------------------------------------------------------------
% summarize results
%----------------------------------------------------------------------
% input:
stats.nspikes           = nspikes;              % [count] number of spikes
stats.nt                = nt;                   % [count] number of trials
stats.Fs                = Fs;                   % [samples/s]
stats.rhalfwin          = rhalfwin;             % [samples]

% filters:
stats.f                 = f;                    % [Z] actual (STA/cSTA/xSTA) filter
stats.fsem              = fsem;                 % [Z] SEM of actual (1st or 2nd order) filter
stats.ftim              = ftim;                 % [s] for filters
stats.f1                = f1;                   % [Z] 1st order (STA) filter
stats.f1sem             = f1sem;                % [Z] SEM of 1st order filter
stats.isiedges          = bedges;               % [samples] ISI bin edges
stats.cpb               = cpb;                  % [count] per ISI bin
stats.sidx              = sidx;                 % [class] ISI class of each spike
stats.optlag            = optlag;               % [samples] optimized (or fixed) lag (how much the filter+reconstruction was shifted) 

% reconstruction quality:
stats.ccLinear          = ccLinear;             % quantifies how bad are saturation effects (linear regression of yhat vs. y)
stats.frq               = frq;                  % [Hz]; frequency bins for the coherence/phase estimates
stats.cohs              = cohs;                 % nfbins x nt; coherence between analog and its reconstruction
stats.phs               = phs;                  % [rad]

% summary stats
stats.cc                = cc;                   % [cc] @ zero-lag, 1 x nt
stats.pvals             = pvals;                % per trial (by randominzation - {shuffle} or jitter)
stats.ccAll             = [ ccAll calc_sem( cc ) ];                % [cc] @ zero-lag, mean and SEM over trials
stats.pval              = pvalAll;              % same, for entire dataset
stats.ccRand            = [ nanmean( rAll ) calc_sem( rAll ) ]; % randomized

stats.ccRel             = ccRel;                % zero-lag CC between reconstructions, mean and SEM
stats.ccRelRand         = ccRelRand;            % randomized

stats.lagsta            = lagSta;               % [samples], here negative are causal
stats.widsta            = widSta;               % [samples]; from the filter
stats.wid               = wid;                  % [samples] FWHH: from the reconstruction (for each trial)
stats.widAll            = widAll;               % [samples] FWHH: from the reconstruction (concatenated data)
stats.widCC2            = widCC2;               % [samples] FWHH: from the CC between reconstructions
stats.widAct            = widAct;               % [samples] from the input

stats.ccmax             = ccmax;                % [cc] either rank or pearson (for each trial)
stats.ccmaxlag          = ccmaxlag - maxlag;    % [samples]; lag of max cc (for each trial)
stats.ccmaxAll          = ccmaxAll;                
stats.ccmaxAllLag       = ccmaxAllLag - maxlag;

% the "most important" results are:
% 1.1. ccAll        Q:  CC (trial-averaged mean and SEM) at zero-lag, for the selected filter (STA/cSTA etc)
% 1.2. pval             pval of the ccAll
% 2. lagsta         T:  [samples] time lag, based on the STA filter
% 3. ccRel          R:  CC (trial-averaged mean, SEM, and p-value) for the selected filter
% 4. halfwid        [samples] half width, based reconstructed concatenated data + input auto-correlation

%----------------------------------------------------------------------
% summarize graphically
%----------------------------------------------------------------------
if graphics 
    ytim                = ( 1 : m )' / Fs;
    
    % layout
    figure
    ax( 1 )             = axes( 'position', [ 0.05  0.75 0.9 0.2 ] );       % spike trains
    ax( 2 )             = axes( 'position', [ 0.05  0.55 0.9 0.2 ] );       % analog + yhat
    ax( 3 )             = axes( 'position', [ 0.05  0.3  0.2 0.2 ] );       % STA
    ax( 4 )             = axes( 'position', [ 0.3   0.3  0.2 0.2 ] );       % (cSTA)
    ax( 5 )             = axes( 'position', [ 0.525 0.3  0.2 0.2 ] );       % y-yhat xcc
    ax( 6 )             = axes( 'position', [ 0.75  0.3  0.2 0.2 ] );       % y-yhat scatter
    ax( 7 )             = axes( 'position', [ 0.05  0.05 0.2 0.2 ] );       % 2D optimization
    ax( 8 )             = axes( 'position', [ 0.3   0.05 0.2 0.2 ] );       % 1D optimization
    ax( 9 )             = axes( 'position', [ 0.525 0.05 0.2 0.2 ] );       % coherence
    ax( 10 )            = axes( 'position', [ 0.75 0.05 0.2 0.2 ] );        % pvals
    
    % output - spikes
    subplot( ax( 1 ) )
    plot_raster( x, ytim, [], [], [], colors( 1, : ) );

    % reconstruction example
    subplot( ax( 2 ) )
    [ maxval, maxidx ]  = max( cc );
    line( ytim, yhat( :, maxidx ), 'color', colors( 3, : ) );
    xpos                = min( xlim ) + 0.9 * diff( xlim );
    ypos                = min( ylim ) + 0.9 * diff( ylim );
    th                  = text( xpos, ypos, sprintf( 'tr%d; cc=%0.3g', maxidx, maxval ) ); 
    set( th, 'color', colors( 3, : ) )
    
    % input - analog
    subplot( ax( 2 ) );
    line( ytim, y1( :, maxidx ), 'color', colors( 2, : ) );
    axis tight
    calibration( [ 100 1 ], [], [], { 'ms', zstr } );
    alines( 0, 'y', 'color', blackColor, 'linestyle', '--' );
  
    % filters
    subplot( ax( 3 ) )
    if ismember( ftype,  { 'cwiener', 'xwiener', 'csta', 'sta2', 'xsta', 'sta3' } )
        patch_band( ftim, f1, f1sem );
        title( sprintf( 'STA (%d spks; W=%0.3g ms)', stats.nspikes, stats.widsta / stats.Fs / 2 ) )
    else
        patch_band( ftim, f, fsem );
        title( sprintf( '%s (%d spks; W=%0.3g ms)', upper( ftype ), stats.nspikes, stats.widsta / stats.Fs / 2 ) )
    end
    axis tight
    if ismember( ftype,  { 'cwiener', 'xwiener', 'csta', 'sta2', 'xsta', 'sta3' } )
        subplot( ax( 4 ) )
        imagesc( ftim, 1 : nibins, f' ), axis xy
        set( ax( 4 ), 'ytick', ( 1 : nibins ) - 0.5, 'yticklabel', bedges( 1 : nibins ) / Fs )
        hold on
        plot( -mean( bedges( [ 1 : end - 1; 2 : end ] ) ) / Fs * 2, ( 1 : nibins ), '.-k' ) % twice the mean ISI
        alines( 0, 'x', 'color', whiteColor, 'linestyle', '--' );
        title( sprintf( '%s (%d-%d spks)', [ ftype( 1 ) upper( ftype( 2 : end ) ) ] ...
            , min( stats.cpb ), max( stats.cpb ) ) )
    else
        ax( 4 ) = -ax( 4 );
    end
    
    % filter optimization (time lags)
    if length( lagsr ) == 1
        ax( [ 7 8 ] ) = -ax( [ 7 8 ] );
    else
        subplot( ax( 7 ) )
        imagesc( lags / Fs, ftim, cci ), axis xy
        title( sprintf( 'CC; optimal lag: %0.2g', optlag / Fs ) )
        alines( optlag / Fs, 'x', 'color', redColor, 'linestyle', '--' );
        
        subplot( ax( 8 ) )
        plot( lags / Fs, max( cci, [], 1 ) )
        title( 'max CC/lag' )
    end

    % input-reconstruction CC
    subplot( ax( 5 ) )
    line( ftim, cct, 'color', grayColor )
    axis tight
    line( ftim, scaleto( act, ylim ), 'color', purpleColor, 'linewidth', 2 )
    line( ftim, cct, 'color', grayColor )
    line( ftim, cct( :, maxidx ), 'color', colors( 3, : ) );
    line( ftim, nanmean( cct, 2 ), 'color', blackColor, 'linewidth', 2 );
    xlim( minmax( ftim ) )
    title( sprintf( '%s: %0.2g (p=%0.2g)', ccstr, ccAll, stats.pval ) )
    
    % scattergram
    subplot( ax( 6 ) )
    plot( ysct( :, 1 ), ysct( :, 2 ), '.k' )
    axis tight
    title( sprintf( '%0.3g', ccLinear ) )
    xlabel( 'yhat', 'color', colors( 3, : ) )
    ylabel( '<y>|yhat', 'color', colors( 2, : ) )
    set( ax( 6 ), 'yAxisLocation', 'right' )

    % coherence
    subplot( ax( 9 ) )
    imagesc( frq, 1 : nt, cohs' ), axis xy
    mcoh                = nanmean( cohs, 2 );
    [ maxcoh, maxidx ]  = max( mcoh ); 
    line( frq, scaleto( mcoh, ylim ), 'color', colors( 2, : ), 'linewidth', 2 );
    xlabel( 'Freq [Hz]' )
    title( sprintf( '<Coh (%0.2gHz)> %0.2g', frq( maxidx ), maxcoh ) )

    % pvalues (shuffles)
    subplot( ax( 10 ) )
    [ hh, bb ]          = hist( pvals, 0.005 : 0.01 : 1 );
    bar( bb, hh, 'EdgeColor', colors( 2, : ), 'FaceColor', colors( 2, : ) );
    xlim( [ 0 1 ] + [ -1 1 ] * 0.01 )
    ylim( [ 0 nt ] )
    xlabel( 'p-value (CC)' )
    alines( 0.05, 'x', 'color', blackColor, 'linestyle', '--' );
    
    % organize
    for i               = 1 : length( ax )
        if ax( i ) < 0
            subplot( -ax( i ) )
            axis off
            continue
        end
        subplot( ax( i ) );
        set( gca, 'tickdir', 'out', 'box', 'off' )
        if i == 1 || i == 2
            set( gca, 'ytick', [] )
            xlim( minmax( ytim ) )
            axis off
        end
        if ismember( i, [ 3 5 6 7 ] )
            alines( 0, 'y', 'color', blackColor, 'linestyle', '--' );
            alines( 0, 'x', 'color', blackColor, 'linestyle', '--' );
        elseif i == 8
            alines( 0, 'x', 'color', blackColor, 'linestyle', '--' );
        end
        
    end
    textf( 0.5, 0.975, replacetok( tstr, '\_', '_' ) );
    colormap( myjet )
end

return

% EOF
