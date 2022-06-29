% calc_mr               calculate multiple regression
%
% call                  [ pcc, cc, c, beta, R2 ] = calc_mr( y, x, rankf, nreps  )
%
% gets                  y           system output, n x 1
%                       x           system input, n x m
%                       rankf       {1}; flag for rank correlation
%                       nreps       {100}; number of bootstrap for SEM of contribution and R2
%
% returns               pcc         partial cc      m x 3: [ pcc    SD  p-values ]
%                       cc          cc              m x 3: [ cc     SD  p-values ]
%                       c           contribution    m x 2: [ c      SD ]
%                       beta        beta            m x 3: [ beta   SD  p-values ]
%                       R2          R2              1 x 3: [ R2     SD  p-value  ]
%
% calls                 calc_pcc, calc_pearson, rankcols, mixmat, calc_sem

% 15-nov-20 ES & HS

% revisions
% 16-nov-20 modified output to include SD

function [ pcc, cc, c, beta, R2 ] = calc_mr( y, x, rankf, nreps )

% arguments
nargs = nargin;
if nargs < 2 || isempty( y ) || isempty( x )
    error( 'missing inputs' )
end
[ n, m ]                        = size( x );
if n ~= size( y, 1 )
    error( 'input size mismatch' )
end
if nargs < 3 || isempty( rankf )
    rankf                       = 1;
end
if nargs < 4 || isempty( nreps )
    nreps                       = 100;
end

% rank
if rankf
    y                           = rankcols( y );
    x                           = rankcols( x );
end

% partial correlation, standardized regression 
[ pcc, beta ]                   = calc_pcc( y, x );

% correlation 
cc                            	= calc_pearson( y * ones( 1, m ), x )';

% significance (correlations, partial correlations)
pvals_cc                        = NaN( m, 1 );
pvals_pcc                       = NaN( m, 1 );
for i                           = 1 : m
    [ ~, pvals_cc( i) ]         = partialcorr( y, x( :, i ), ones( n, 1 ) );
    cidx                        = setdiff( 1 : m, i );
    [ ~, pvals_pcc( i ) ]       = partialcorr( y, x( :, i ), x( :, cidx ) );
end

% contribution and R2
c                               = beta .* cc; % C may be negative for some of the regressors
R2                              = sum( c );

% compute SEM for the contribution and R2
mix                             = mixmat( ( 1 : n )' * ones( 1, nreps ), 1, 2 );
cc_bs                           = NaN( m, nreps );
pcc_bs                          = NaN( m, nreps );
beta_bs                         = NaN( m, nreps );
for i                           = 1 : nreps 
    ridx                        = mix( :, i );
    [ pcc_bs( :, i ) ...
        , beta_bs( :, i ) ]     = calc_pcc( y( ridx, : ), x( ridx, : ) );
    cc_bs( :, i )               = calc_pearson( y( ridx, : ) * ones( 1, m ), x( ridx, : ) )';
end
sd_cc                           = calc_sem( cc_bs, 2 );
sd_pcc                          = calc_sem( pcc_bs, 2 );
sd_beta                         = calc_sem( beta_bs, 2 );
c_bs                            = pcc_bs .* beta_bs;
sd_c                            = calc_sem( c_bs, 2 );
sd_R2                           = calc_sem( sum( c_bs, 1 ) );

% p-value for full model
[ ~, ~, ~, ~, stts ]            = regress( y, [ ones( n, 1 ) x ] );

% organize output
cc                              = [ cc sd_cc pvals_cc ];
pcc                             = [ pcc sd_pcc pvals_pcc ];
beta                            = [ beta sd_beta pvals_pcc ];
c                               = [ c sd_c ];
R2                              = [ R2 sd_R2 stts( 3 ) ];

return

% EOF
