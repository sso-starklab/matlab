% optclust      optimal clustering from N-dim data by enumeration
% 
% [ clu, centroids, dispersion ] = optclust( x, mode, graphics )
%
% ARGUMENTS
% x             data; n x p
% mode          optimization mode (also determines clustering method):
%                   {'SIC'} uses a Gaussian mixture model
%                   'AIC' (not recommended), also uses GMM
%                   'ratio' uses kmeans and comnputes the ratio between 
%                       mean intra-cluster distance and min inter-cluster distance
%                   'max' (kmeans) is more robust but cannot detect k<4
% graphics      {0}; plots the measure + the clustered data
%                   for >2D, only the first 2 dimensions are plotted
%
% OUTPUT
% clu           the labels
% centroids     the cluster centers (means)
% dispersion    either covariance matrix (GMM) or within-cluster
%               point-to-centroid squared euclidean distances (KMEANS)
%
% ALGORITHM
% step 1: search
%       for each k, clusters the data using a GMM (or k-means) and computes a
%       goodness of fit measure (e.g. BIC)
% step 2: optimization
%       chooses k that minimizes the goodness of fit measure
% step 3: post processing
%       the resulting cluster are then sorted according to the distance of the
%       centroid from the origin
%
% ADDITIONAL PARAMETERS allow some control over the internal behavior:
% nc            range to consider; default is [ 1 kmax ], where
%                   kmax = min( K, length( unique( x, 'rows' ) ), 
%                   where K = ceil( 0.1 * n ) (if n<100),
%                   K=10 (if n<1000), or 20 (if n>=1000)
% ni            number of iterations for each k; default is [ 10 100 ] for
%                   k-means and [ 2 10 ] for EM, for the search and summary steps
% addnull       flag for adding a null cluster; without it, k-means cannot detect k<2
%                   (inter-cluster distance is undefined for a singleton)
%
% ADVICE
% (1) For a Gaussian mixture, this algorithm when used with the SIC yields
% errors which are close to the Bayes error. The AIC (or the K-means based
% method) tend to overcluster. Thus although when data are non-gaussian or 
% the Bayes error is large (clusters overlap) the k-means based methods may
% be more appropriate, for most purposes the SIC mode provides the best
% results
% (2) There is a random initialization step in the clustering; to get
% replicable results, set rng('default'); alternatively to check
% robustness, set rng('shuffle'); before running this routine
%
% PERFORMANCE
% This is a slow routine, it is useful for small n (up to ~1000 elements)
% or for small nc( 2 ) (internally limited to 20 by default)
%
% CALLS         MATLAB functions
% 
% see also      gmdistribution, kmeans

% 20-jan-13 ES

% revisions:
% 21-jan-13 (1) added noise cluster option
%           (2) sort the centroids by distance from 0
%           (3) integrate GMM and BIC
% 01-aug-19 (1) added warning off for EmptyClusterRep (necessary in R2018a)
% 17-aug-19 cleaned up
% 13-sep-10 warning off for FailedToConvergeReps

function [ clu, centroids, dispersion ] = optclust( x, mode, graphics, nc, ni, nul, vflag )

%----------------------------------------------------------------------%
% initialization
%----------------------------------------------------------------------%

mfname                      = upper( mfilename );
warning( 'off', 'stats:kmeans:EmptyCluster' )
warning( 'off', 'stats:gmdistribution:FailedToConverge' )
warning( 'off', 'stats:gmdistribution:IllCondCov' )
warning( 'off', 'stats:kmeans:FailedToConverge' )
warning( 'off', 'stats:kmeans:FailedToConvergeReps' )
warning( 'off', 'stats:kmeans:EmptyClusterRep' )

clu                         = [];
centroids                   = [];
dispersion                  = [];

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( x )
    return, 
end
[ m, n ]                    = size( x );
if m == 1 || length( unique( x( : ) ) ) == 1
    clu                     = ones( m, 1 );
    centroids               = x;
    dispersion              = 0;
    return
end
if nargs < 2 || isempty( mode )
    mode                = 'bic';
end
mode                    = lower( mode );
switch mode
    case { 'ratio', 'max' }
        core            = 'KMEANS';
    case { 'aic', 'bic', 'sic' }
        core            = 'GMM';
    otherwise
        return
end
if nargs < 3 || isempty( graphics )
    graphics            = 0;
end

if nargs < 4 || isempty( nc ) || length( nc ) ~= 2
    nc                  = [];
end
if nargs < 5 || isempty( ni ) || length( ni ) ~= 2
    ni                  = [];
end
if nargs < 6 || isempty( nul )
    nul                 = 1;
end
if nargs < 7 || isempty( vflag )
    vflag               = 0;
end

% determine enumeration parameters
if isempty( nc ) 
    if nul
        nc( 1 )         = 1;
    else
        nc( 1 )         = 2;
    end
    if m < 100
        nc( 2 )         = max( ceil( 0.1 * m ), 10 );
    elseif m < 1000
        nc( 2 )         = 10;
    else
        nc( 2 )         = 20;
    end
end
nc                      = round( nc );
nc( 1 )                 = max( nc( 1 ), 1 );
nc( 2 )                 = min( nc( 2 ), length( unique( x, 'rows' ) ) );
if isempty( ni )
    switch mode
        case { 'ratio', 'max' }
            ni          = [ 10 100 ];
        case { 'aic', 'bic', 'sic' }
            ni          = [ 10 20 ];
    end
end
ni                      = round( ni );
if nul
    z                   = zeros( 1, n );
else
    z                   = [];
end

%----------------------------------------------------------------------%
% search stage
%----------------------------------------------------------------------%
verb( sprintf( '%s: Initiating search stage (%s): m=%d; n=%d; k=%d:%d; Using '...
    , mfname, core, m, n, nc( 1 ), nc( 2 ) ), -vflag )
r                       = NaN * ones( nc( 2 ), 1 );
for k                   = nc( 1 ) : nc( 2 )
    verb( sprintf( '%d', k ), -vflag )
    switch core
        case 'KMEANS'
            try
                [ clu, centroids, d1 ] = kmeans( x, k...
                    , 'Start', 'Uniform'...
                    , 'Distance', 'sqEuclidean'...
                    , 'EmptyAction', 'drop'...
                    , 'Display', 'Off'...
                    , 'Replicates', ni( 1 ) );
            catch ME
                verb( '! ', -vflag )
                continue
            end
            %d2 = pdist( centroids, 'euclidean' ) .^ 2;
            d2          = pdist( [ centroids; z ], 'euclidean' ) .^ 2;
            r( k )      = mean( d1( ~sum( isnan( centroids ), 2 ) ) ) / min( d2 );
        case 'GMM'
            try
                % three options:
                % -initialize randomly with k samples, many times, keep the best (LL) model
                % -initialize with a random subset; then compute the mean and initialize once with this
                if k == 1
                    gm          = gmdistribution.fit( x, k, 'Replicates', ni( 1 ) );
                else
                    mukeep      = [];
                    sigkeep     = [];
                    NlogL       = [];
                    for i       = 1 : ni( 1 )
                        ilabels     = ceil( rand( m, 1 ) * k );
                        gm          = gmdistribution.fit( x, k, 'Start', ilabels );
                        [ mukeep( :, :, i ), sidx ] = sortrows( gm.mu, 1 );
                        sigkeep( :, :, :, i )       = gm.Sigma( :, :, sidx );
                        NlogL( i )                  = gm.NlogL;
                        bic( i )                    = gm.BIC;
                    end
                    s.mu        = mean( mukeep, 3 );
                    s.Sigma     = mean( sigkeep, 4 );
                    gm          = gmdistribution.fit( x, k, 'Start', s );
                end
                gmm{ k }        = gm;
            catch ME
                verb( '! ', -vflag )
                continue
            end
            if gm.Converged
                switch mode
                    case 'aic'
                        r( k )  = gm.AIC;
                    case { 'bic', 'sic' }
                        r( k )  = gm.BIC;
                end
            end
    end
    verb( sprintf( ' ' ), -vflag )
end

%----------------------------------------------------------------------%
% optimization
%----------------------------------------------------------------------%
k                           = [];
if strcmp( mode, 'max' ) % the minimum following the first local maximum
    d                       = diff( sign( diff( r ) ) );
    [ row, col ]            = find( d < -1 );
    if ~isempty( row )
        rm                  = row( 1 ) + 1;
        [ minval, khat ]    = min( r( rm : end ) );
        k                   = rm( 1 ) + khat - 1;
    end
end
if isempty( k ) % global minimum
    [ minval, k ]           = min( r );
end

%----------------------------------------------------------------------%
% post-processing (recompute and sort centroids)
%----------------------------------------------------------------------%
verb( sprintf( 'clusters... %d clusters chosen! Summarizing...', k ), -vflag )

switch core
    case 'KMEANS'
        [ clu, centroids, dispersion ] = kmeans( x, k...
            , 'Start', 'Uniform'...
            , 'Distance', 'sqEuclidean'...
            , 'EmptyAction', 'drop'...
            , 'Display', 'Off'...
            , 'Replicates', ni( 2 ) );
    case 'GMM'
        % two options: 
        % initialize with the previous mean gmm (may lead to overfitting)
        % or just use the selected model dimension (k) and intialize
        % randomly
        gm                  = gmdistribution.fit( x, k, 'Replicates', ni( 2 ) );
        % the following may lead to overfitting:
%         s.mu = gmm{ k }.mu; 
%         s.Sigma = gmm{ k }.Sigma; 
%         gm2 = gmdistribution.fit( x, k, 'Start', s );
        centroids           = gm.mu;
        dispersion          = gm.Sigma;
        %clu = cluster( gm, x );
        p                   = posterior( gm, x ); 
        [ minp, clu ]       = max( p, [], 2 );
end

[ sc, sidx ]                = sort( sum( centroids .^ 2, 2 ) );
reorder                     = [ sidx ( 1 : k )' ];
clu0                        = clu;
for i                       = 1 : k
    clu( clu0 == reorder( i, 1 ) ) = reorder( i, 2 );
end
centroids                   = centroids( sidx, : );
switch core
    case 'KMEANS'
        dispersion          = dispersion( sidx, : );
    case 'GMM'
        dispersion          = dispersion( :, :, sidx );
end
verb( sprintf( 'Done!' ), vflag )

%----------------------------------------------------------------------%
% plot
%----------------------------------------------------------------------%
if graphics 
    
    ah( 1 )                 = subplot( 2, 1, 1 );
    cla
    switch mode
        case { 'ratio', 'max' }
            label           = '1 / ratio';
            plot( 1 : nc( 2 ), 1./r, 'o-' )
        case { 'bic', 'sic' }
            label           = 'BIC';
            plot( 1 : nc( 2 ), r, 'o-' )
        case 'aic' 
            label           = 'AIC';
            plot( 1 : nc( 2 ), r, 'o-' )
    end
    ylabel( label )
    xlabel( 'k' )
    set( ah( 1 ), 'yscale','log', 'box', 'off', 'tickdir','out' ),
    
    ah( 2 )                 = subplot( 2, 1, 2 );
    cla
    map                     = lines( k );
    idx                     = 1 : m;
    for ki                  = 1 : k
        hold on,
        if n == 1
            ph              = plot( idx( clu == ki ), x( clu == ki, 1 ), 'o' );
            if ki == 1
                xlabel( 'Sample number' )
                ylabel( 'x(1)' )
            end
        elseif n == 2
            ph              = plot( x( clu == ki, 1 ), x( clu == ki, 2 ), 'o' );
            if ki == 1
                xlabel( 'x(1)' )
                ylabel( 'x(2)' )
                axis square
            end
        end
        set( ph, 'color', map( ki, : ) );
    end
    title( sprintf( 'Selection: %d clusters; %s: %0.3g', k, label, r( k ) ) );
    set( ah( 2 ), 'box', 'off', 'tickdir','out' ),
    
end

warning( 'on', 'stats:kmeans:EmptyCluster' )
warning( 'on', 'stats:gmdistribution:FailedToConverge' )
warning( 'on', 'stats:gmdistribution:IllCondCov' )
warning( 'on', 'stats:kmeans:FailedToConverge' )

return

% EOF

%----------------------------------------------------------------------%
% testing:
%----------------------------------------------------------------------%

% generate a Gaussian mixture:
n = 100;
mu1 = [1 3]; Sigma1 = [1 .5; .5 2]; 
%mu2 = [7 8 ]; Sigma2 = [1 .5; .5 2];
mu2 = [3 2 ]; Sigma2 = [1 .5; .5 2];
mu3 = [5 1]; Sigma3 = [1 .5; .5 2];
% z1 = repmat(mu1,n,1) + randn(n,2)*chol(Sigma1);
% z2 = repmat(mu2,n,1) + randn(n,2)*chol(Sigma2);
% z3 = repmat(mu3,n,1) + randn(n,2)*chol(Sigma3);
% z = [ z1; z2; z3 ]; 
% cluGround = [ ones( n, 1 ); 2 * ones( n, 1 ); 3 * ones( n, 1 ) ];

% compute the Bayes error for this mixture:
mu = [mu1;mu2;mu3];
Sigma = cat(3,Sigma1,Sigma2,Sigma3);
mixp = ones(1,3)/3;
gm = gmdistribution(mu,Sigma,mixp);
[ z cluGround ] = random( gm, 300 );
cluOpt = cluster( gm, z );

centroids = gm.mu;
k = size( centroids, 1 );
[ sc sidx ] = sort( sum( centroids .^ 2, 2 ) );
reorder = [ sidx [ 1 : k ]' ];
clu0 = cluOpt;
clu1 = cluGround;
for i = 1 : k
    cluOpt( clu0 == reorder( i, 1 ) ) = reorder( i, 2 );
    cluGround( clu1 == reorder( i, 1 ) ) = reorder( i, 2 );
end

figure( 1 ),
map = lines( k );
for ki = 1 : k
    sidx = cluOpt == ki;
    hold on
    ph = plot( z( sidx, 1), z( sidx, 2 ), '.' );
    set( ph, 'color', map( ki, : ) );
end
fprintf( 1, 'BAYES error: %0.3g\n', sum( cluOpt ~= cluGround ) / length( cluOpt ) )

% now cluster the data and compute the true error:
figure( 2 ), clu = optclust( z, 'ratio', 1 );
fprintf( 1, 'Clustering error: %0.3g\n', sum( clu ~= cluGround ) / length( clu ) )
figure( 3 ), clu = optclust( z, 'bic', 1 );
fprintf( 1, 'Clustering error: %0.3g\n', sum( clu ~= cluGround ) / length( clu ) )

% 
% figure, Z = linkage( z, 'ward', 'euclidean' );  dendrogram( Z )

