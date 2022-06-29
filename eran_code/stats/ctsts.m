% CTSTS         calculate time-series trend and significance
%
% call          [ r0, rn, p ]  = ctsts( x, n, g )
%
% gets          x       data vector
%               n       number of resampling trials
%               g       graphics flag
%
% returns       r0      correlation coefficient
%               rn      resampled (bootstraped) correlation coefficient
%               p       p-value
%
% calls         calc_p, calc_pearson, mixmat

% 10-jun-06 ES

% revisions
% 14-feb-07 corrections - mixmat called without replacement
% 28-oct-14 superficial modifications
% 16-sep-19 cleaned up

function [ r0, rn, p ]  = ctsts( x, n, g )

r0                      = [];
rn                      = [];
p                       = [];

nargs                   = nargin;
if nargs < 1 || isempty( x )
    return
end
x                       = x( : );

d1                      = ( 1 : length( x ) )'; 
r0                      = calc_pearson( x( : ), d1 );
if exist( 'n', 'var' ) && n( 1 ) > 0
    d2                  = ones( 1, n( 1 ) ); 
    xn                  = mixmat( x * d2, 1, 1 ); 
    rn                  = calc_pearson( xn, d1 * d2 );
    p                   = calc_p( rn.^2, r0.^2, 0 );
else
    rn                  = NaN;
    p                   = NaN;
end

if exist( 'g', 'var' ) && g
    figure
    subplot( 2, 1, 1 )
    plot( x )
    title( texlabel( sprintf( 'rho = %0.3g', r0 ) ) )
    if n > 0
        subplot( 2, 1, 2 )
        hist( rn, 20 )
        separators( r0 ); 
        xlim( [ -1 1 ] )
        title( sprintf( 'p: %0.3g', p ) )
    end
end

return

% EOF
