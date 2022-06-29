% BINOMIAL_INLINES        exact and approximate (normal) binomial confidence limits.
%
% [ bino_ci_exact, bino_ci_norm, bino_se_norm ] = binomial_inlines
%
% first two inlines get ALPHA, K, N (tails, successes, total)
% and return [lower point upper] estimates of prob
% 
% the third gets K, N and returns only the SEM
%
% all assume that ALPHA is a scalar (0-1) and K/N are integers (vectors)

% 24-mar-13 ES

function [ bino_ci_exact, bino_ci_norm, bino_se_norm ] = binomial_inlines

bino_ci_exact = inline( '[ binoinv( alpha( 1 ) / 2, n( : ), k( : ) ./ n( : ) ) k( : ) binoinv( 1 - alpha( 1 ) / 2, n( : ), k( : ) ./ n( : ) ) ] ./ ( n( : ) * [ 1 1 1 ] )'...
    , 'alpha', 'k', 'n' );
bino_ci_norm = inline( '( k( : ) ./ n( : ) ) * [ 1 1 1 ] + ( sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : ) ) * norminv( 1 - alpha( 1 ) / 2 ) * [ -1 0 1 ]'...
    , 'alpha', 'k', 'n' );
bino_se_norm = inline( 'sqrt( k( : ) .* ( 1 - k( : ) ./ n( : ) ) ) ./ n( : )', 'k', 'n' );

return