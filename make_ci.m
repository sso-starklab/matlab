% MAKE_CI           make cross-indices for vectors of n1 and n2 elements.
%
% call              [ I J ] = MAKE_CI( N1, N2 )
%
% gets              N1, N2      number of elements in each group
%
% returns           I, J        columns vectors of group indices
%                   I           of the group with N1 elements
%                   J           of the group with N2 elements
%
%                               if only I is requested, then a matrix of
%                               two columns is returned.
%
% calls             nothing
%
% see also MAKE_AI (for one group), MAKE_MI (for more than two groups).

% 08-apr-04 ES

function [ i, j ] = make_ci( n1, n2 )

if nargin < 2 |...
        ~isa( n1, 'double' ) | prod( size( n1 ) ) ~= 1 | any( n1 ~= round( n1 ) )...
        ~isa( n2, 'double' ) | prod( size( n2 ) ) ~= 1 | any( n2 ~= round( n2 ) )
    error( '2 integer arguments' )
end

if ~n1 | ~n2
    i = [];
    j = [];
    return
end

i = zeros( n1 * n2, 1 );
pi = 1 : n2 : n1 * n2;
i( pi ) = 1;
i = cumsum( i );

pj = ( 1 : n2 )';
j = pj( :, ones( 1, n1 ) );
j = j( : );

if nargout <= 1
    i = [ i j ];
end

return
