% MAKE_AI           make auto-indices for vector of n elements.
%
% call              [ I J ] = MAKE_AI( N )
%
% gets              N           number of elements in each group
%
% returns           I, J        columns vectors of group indices
%                               if only I is requested, then a matrix of
%                               two columns is returned.
%
% calls             nothing
%
% see also MAKE_CI (for two groups), MAKE_MI (for more than two groups).

% 08-apr-04 ES

function [ i, j ] = make_auto_indices( n )

if nargin < 1 |...
        ~isa( n, 'double' ) | prod( size( n ) ) ~= 1 | any( n ~= round( n ) )
    error( '1 integer argument' )
end

if ~n
    i = [];
    j = [];
    return
end

m = n * ( n - 1 ) / 2;
i = zeros( m, 1 );
j = ones( m, 1 );

p = ( n - 1 ) : -1 : 2;
i( cumsum( [ 1 p ] ) ) = 1;
i = cumsum( i );

j( cumsum( p ) + 1 ) = 2 - p;
if n > 1
    j( 1 )=2;
end
j = cumsum( j );

if nargout <= 1
    i = [ i j ];
end

return

