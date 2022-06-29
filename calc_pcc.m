% CALC_PCC          calculate all maxmial order partial correlation coefficients and standardized regression coefficients
%
% call              [ PCC, BETA ] = CALC_PCC( Y, X )
%
% gets              Y           system output, n x 1
%                   X           system input, n x m
%
% returns           PCC         all (m-1)th order partial CCs, m x 1
%                   BETA        all standardized regression coefficients, m x 1
%
% calls             nothing
%
% notes             1. no error checking. assume input size as above; no NaNs
%                   2. To compute contribution: 
%
%                           Ci = beta_i * rho_i 
%
%                       where rho_i is correlation coefficient, obtained using
%
%                           rho = calc_pearson( y * ones( 1, size( x, 2 ) ), x )
%
%                       and the sum of all contributions is the model R2:
%
%                           R2 = sum( Ci )
%
%                   3. to compute Spearman (rank) correlations, apply ranking to the inputs: 
%
%                           calc_pcc( rankcols( y ), rankcols( x ) )
%
% see also calc_pearson, calc_spearman, rankcols

% 22-mar-06 ES

% 10-nov-20 (1) cleaned up
%           (2) added beta to output

function [ rp, beta ] = calc_pcc( y, X )

[ n, p ]       	= size( X );
y               = y - sum( y ) / n;
s               = sqrt( sum( y .* y, 1 ) / ( n - 1 ) );
if s == 0
    y           = zeros( n, 1 ); 
else
    y           = y / s; 
end
m               = sum( X, 1 ) / n;
d               = ones( n, 1 );
X               = X - d * m;
s               = sqrt( sum( X .* X, 1 ) / ( n - 1 ) );
zidx            = s == 0;
den             = d * s;
Z               = zeros( n, p );
Z( :, ~zidx )   = X( :, ~zidx ) ./ den( :, ~zidx );
I               = eye( p );
[ Q, R ]        = qr( Z, 0 );
beta            = ( R \ ( Q' * y ) );
Tstat           = beta ./ sqrt( diag( R \ I * ( R \ I )' ) ) / norm( y - Z * beta ) * sqrt( n - p );
rp              = Tstat ./ ( sqrt( Tstat .^ 2 + n - p ) );

return

% EOF
