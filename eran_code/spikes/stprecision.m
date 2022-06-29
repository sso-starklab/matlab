% stprecision               estimate reliability (and precision) of single-cell spike trains
%
% call                      [ rel, prec, cc, ccSD ] = stprecision( x, yhat )
%
% gets                      x               multi-trial spike trains (sparse matrix, n x nt)
%                           yhat            a specific reconstruction (same dimensions as x)
% 
% optional arguments (given as name/value pairs):
%
%                           s               {[ 1 2 4 8 16 32 64 ]}  m-element vector of Gaussian SDs [samples]
%                           ftype           { 'gauss' }             filter type; 'exp', 'alpha', and 'rect' are also supported
%                           graphics        {0}                     flag
%
% returns                   rel             mean, SEM, and p-value of all cc (for all possible yhat pairs); p-value by one-sided signed-rank test
%                           prec            precision, [samples]; computed only when s includes multiple SDs
%                           cc              the individual n*(n-1)/2 values (for all pairs in the reconstruction)
%                           ccSD            for all pairs, and all jitters (m by n*(n-1)/2 matrix)
%
% calls                     makeblocks, minmax, ParseArgPairs           (general)
%                           calc_sem, calc_spearman, make_ai            (stats)
%                           firfilt, makegaussfir                       (ssp) 
%                           alines, patch_band                          (graph)
%
% does                      answers the following three questions:
%                               1. how reliable is the reconstruction yhat
%                               2. how reliable are the spike trains x (in terms of reconstruction with a non data-dependent filter)
%                               3. how precise temporally are the spike trains (computed by comparing the reconstructions of yhat and x)
%
% definitions:              reliability := reproducibility of the same spike train
%                                   this is similar to the Van Rossum metric (Neural Computation 2001) 
%                                   and identical to the Schreiber et al. metric (Neurocomputing, 2003).
%
%                           precision   := the temporal jitter between spike trains
%                                   to determine the precision, we compute reliability using yhat, which 
%                                   may be generated by reconstruction with a data-dependent filter which 
%                                   is not necessarily well-defined in terms of precision, but gives a fixed reliability.
%
% NOTE:                     can use the cc to classify the spike trains into subsets (e.g. using hierarchical clustering). 

% 27-oct-14 ES

% revisions
% 28-oct-14 reliability pvalue computed
% 07-jan-15 also supports quick cimputation if s is NaN (no smoothing)
% 21-jun-15 graphics( 1 ) == 1 for plotting
% 14-oct-19 cleaned up and documented

function [ rel, prec, cc, ccSD ] = stprecision( x, yhat, varargin )

% block-wise to prevent memory exhaustion:
BLOCKSIZE                   = 50;

% arguments
nargs                       = nargin;
if nargs < 2 || isempty( x ) || isempty( yhat )
    return
end
[ s, ftype ...
    , graphics ]            = ParseArgPairs(...
    { 's', 'ftype' ...
    , 'graphics' }...
    , { 2 .^ ( 0 : 6 ), 'gauss' ...
    , 0 }...
    , varargin{ : } );
nt                          = size( x, 2 );
if ~isequal( size( x ), size( yhat ) )
    error( 'input size mismatch' )
end

% preps
aidx                        = make_ai( nt );
npairs                      = size( aidx, 1 );
blocks                      = makeblocks( size( aidx, 1 ), BLOCKSIZE, 0 );
nblocks                     = size( blocks, 1 );

% compute reliability for the reconstruction:
cc                          = zeros( npairs, 1 );
for i                       = 1 : nblocks
    bidx                    = blocks( i, 1 ) : blocks( i, 2 );
    cc( bidx )              = calc_spearman( yhat( :, aidx( bidx, 1 ) ), yhat( :, aidx( bidx, 2  ) ) );
end

% compute p-value
relPval                     = signrank( cc, 0 );
if nanmedian( cc ) > 0
    relPval                 = relPval / 2;
else
    relPval                 = 1 - relPval / 2;
end
rel                         = [ nanmean( cc ) calc_sem( cc ) relPval ];

% compute reliability for various temporal jitters:
nSD                         = length( s );
if all( isnan( s ) )
    nSD                     = 0;
end
ccSD                        = NaN * ones( nSD, npairs );
xf                          = full( x );

for k                       = 1 : nSD

    % build the filter
    tau                     = s( k );
    switch lower( ftype )
        case 'gauss'
            iwin            = makegaussfir( tau, 1 );
        case 'alpha'
            t               = 0 : 1 : 5 * tau;
            iwin            = t .* exp( -t / tau );
        case 'exp'
            t               = 0 : 1 : 5 * tau; % 3 is >95% of the support, 5 is >99% of the support
            iwin            = exp( -t / tau );
        case { 'rect', 'ma', 'box', 'boxcar' }
            iwin            = ones( tau, 1 );
        otherwise
            error( 'unsupported filter type' )
    end
    iwin                    = iwin / sum( iwin );
    
    % convolve the spike trains:
    switch lower( ftype )
        case { 'exp', 'alpha' }
            fg              = filter( iwin, 1, xf ); % exponential tail
        case { 'gauss', 'rect', 'ma', 'box', 'boxcar' }
            fg              = firfilt( xf, iwin ); % symmetric around spike
    end
    
    % compute the cc:
    for i                   = 1 : nblocks
        bidx                = blocks( i, 1 ) : blocks( i, 2 );
        ccSD( k, bidx )     = calc_spearman( fg( :, aidx( bidx, 1 ) ), fg( :, aidx( bidx, 2  ) ) );
    end
    
end

% compute the precision (closest reliability to the reconstruction):
[ ~, minidx ]               = min( abs( mean( ccSD, 2 ) - rel( 1 ) ) );
if minidx > 1
    minidx                  = minidx + [ -1 0 ];
end
prec                        = geomean( s( minidx ) ); % [samples]

if graphics
    
    newplot
    mm                      = mean( ccSD, 2 );
    ss                      = calc_sem( ccSD, 2 );
    if nSD ~= 0
        plot( s, ccSD, 'color', [ 1 1 1 ] * 0.7 )
        patch_band( s, mm, ss );
        xlim( minmax( s ) )
        set( gca, 'xscale', 'log' )
        set( gca, 'xtick', s, 'xticklabel', s )
        xlabel( 's [samples]' )
    end
    ylabel( 'CC-vr' )
    alines( rel( 1 ), 'y', 'color', [ 1 0 0 ], 'linestyle', '-' );
    alines( rel( 1 ) + [ -1 1 ] * rel( 2 ), 'y', 'color', [ 1 0 0 ], 'linestyle', '--' );
    axis tight
    alines( prec, 'x', 'color', [ 1 0 0 ], 'linestyle', '--' );
    set( gca, 'tickdir', 'out', 'box', 'off' ),
    title( sprintf( 'Precision=%0.3g samples', prec  ) )
end

return

% EOF
