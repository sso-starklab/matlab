% MAKE_MI           make multiple-indices.
%
% call              I = MAKE_MI( N )
%
% gets              N       a vector of integers, each indicating the
%                           number of members in each group.
%
% returns           I       a matrix of integers, prod( N ) by length( N ),
%                           generated sorted by rows
%
% calls             nothing
%
% see also MAKE_AI (for one group), MAKE_CI (for 2 groups only).

% 09-apr-04 ES

function y = make_mi( n )

% exception handling

if nargin < 1 | ~isa( n, 'double' ) | any( n ~= round( n ) )
    error( '1 vectorial integer argument' )
end
n = n( : );
m = length( n );
if m <= 1
    if m
        y = [ 1 : n ]';
    else
        y = [];
    end
    return
end

% preparations

b( 2 : m + 2 ) = [ 1; cumprod( n ) ];

if b( end ) > 1e6
    error( 'matrix will be unacceptably large' )
end

a( m + 1 : -1 : 2 ) = [ 1; cumprod( n( m : -1 : 2 ) ) ];
n = [ 1; n; 1 ];

% build by columns

for k = 2 : m + 1
    v = ones( a( k ), 1 ) * [ 1 : n( k ) ];         % block of group indices
    t = v( : );
    tt = t( :, ones( 1, b( k ) ) );                 % replicate this block
    y( :, k - 1 ) = tt( : );
end

return