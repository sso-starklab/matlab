% CALC_COM          center of mass of matrix columns
%
% CALL              [ c, sd ] = calc_com( x, w )
% 
% GETS              x     data matrix (vector or matrix)
%                   w     weights (vector or matrix)
%                   if both are matrices, must have same number of columns
%                   nans are supported. 
%
% RETURNS           c     center of mass (1xn)
%                   sd    from center of mass (1xn)
%
% DOES              calc_com( x ) is like [ mean( x ) std( x ) ]
%
% non-weighted:
% c = sum( x ) / n
% sd = sqrt( sum( ( x - c ) ^2 ) )
% 
% weighted:
% c = sum( x * w ) / sum( w )
% sd = sqrt( sum( w * ( w * x - c )^2 ) ) ) )
%
% see also: var

% 10-jun-11 ES

% revisions
% 10-jan-12 also supported vector x, matrix of w
% 17-aug-19 cleaned up
% 19-jul-20 modified sd computation to support NaNs

function [ c, sd ] = calc_com( x, w, flag )

% initialize output
c                   = [];
sd                  = [];

% arguments 
nargs               = nargin;
if nargs < 1 || isempty( x )
    return
end
if nargs < 2 || isempty( w )
    w               = [];
end
if nargs < 3 || isempty( flag )
    flag            = 0;
end

% organize input
[ m, n ]            = size( x );
[ m1, n1 ]          = size( w );
if m == m1 && n == n1
elseif m1 == 1 || n1 == 1
    if m1 == 1 && n1 == m
        w           = w';
        ex          = 1;
    elseif m1 == m && n1 == 1
        ex          = 1;
    elseif m1 == n && n1 == m
        x           = x';
        ex          = 0;
    else
        error( 'input size mismatch' )
    end
    if ex
        w           = w * ones( 1, n );
    end
elseif m == 1 || n == 1
    if m == 1 && n == m1
        x           = x';
    elseif m == m1 && n == 1
    else
        error( 'input size mismatch' )
    end
    x               = x * ones( 1, n1 );
elseif isempty( w )
    w               = ones( m, n );
else
    error( 'input size mismatch' )
end

% compute
w                   = bsxfun( @rdivide, w, nansum( w ) );
c                   = nansum( x .* w );
xhat                = bsxfun( @minus, x, c );
sd                  = sqrt( nansum( bsxfun( @times, xhat .^2, w ) ) );
if ~flag
    sd              = sd * sqrt( m / ( m - 1 ) );
end

return

% EOF
