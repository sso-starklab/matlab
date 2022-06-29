% INRANGE           indices of data in range.
%
% call              Y = INRANGE( X, B, DIM, DTYPE )
%
% gets              X       data
%                   B       bounds (matrix of 2 vectors, each same size
%                           as one of X's dimensions)
%                   DIM
%                   DTYPE   {'linear'} or 'circular'
%
% returns           Y       indices of data in range
%
% calls             nothing
%
% example           
%                   x = rand( 10, 4 )
%                   b = sort( rand( 10, 2 ), 2 )
%                   y = inrange( x, b )
%
%                   returns indices of x that are within the range of b
%                   (rows by row)

% 08-mar-04 ES

% revisions
% 17-apr-04 2nd output added
% 02-may-04 4th input
% 15-jul-15 empty input case handled
% 31-aug-19 cleaned up

function [ y, xin ] = inrange( x, b, dim, datatype )

y                       = [];
xin                     = [];

nargs                   = nargin;
if nargs < 2, error( '2 arguments' ), end
if nargs < 4 || isempty( datatype ), datatype = 'linear'; end

[ m, n ]                = size( x );
[ mb, nb ]              = size( b );
eflag = 0;
if m == 0 || mb == 0
    return
end

if nargs < 3 || isempty( dim )
    if mb == 2
        dim             = 1;
    elseif nb == 2
        dim             = 2;
    else
        dim             = 0;
    end
end
    
if dim == 1
    if n == nb
        bmin            = repmat( b( 1, : ), m, 1 );
        bmax            = repmat( b( 2, : ), m, 1 );
    elseif m == nb
        bmin            = repmat( b( 1, : )', 1, n );
        bmax            = repmat( b( 2, : )', 1, n );
    else
        eflag           = 1;
    end
elseif dim == 2
    if m == mb
        bmin            = repmat( b( :, 1 ), 1, n );
        bmax            = repmat( b( :, 2 ), 1, n );
    elseif n == mb
        bmin            = repmat( b( :, 1 )', m, 1 );
        bmax            = repmat( b( :, 2 )', m, 1 );
    else
        eflag           = 1;
    end
else
    eflag               = 1;
end

if eflag
    error( 'input size mismatch' )
end

if strcmp( lower( datatype( 1 : 4 ) ), 'circ' )
    x                   = mod( x, 2 * pi );
    y                   = ~true( m, n );
    idx                 = bmin <= bmax;
    y( idx )            = x( idx ) >= bmin( idx ) & x( idx ) <= bmax( idx );
    idx                 = bmin > bmax;
    y( idx )            = ( x( idx ) >= bmin( idx ) & x( idx ) <= 2 * pi ) |...
        ( x( idx ) >= 0 & x( idx ) <= bmax( idx ) );
else
    y                   = x >= bmin & x <= bmax;
end

if nargout > 1
    xin                 = x( y );
end

return

% EOF
