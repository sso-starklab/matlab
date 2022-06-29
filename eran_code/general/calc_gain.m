% calc_gain         ratio between two numbers
%
% call              g = calc_gain( a, b, d )
% 
% gets              a, b            numbers
%                   d               {1}, exponent of 10 for limiting the
%                                   range (use inf to get standard behavior)
%
% does              computes ratio between two numbers (g=a/b), with corrections for zero cases
%
%                   a=0, b=0:     gain := 1           (matlab: 0/0=NaN)
%                   a~=0, b=0:    gain := 10^(d)      (matlab: x/0=inf)
%                   a=0, b~=0:    gain := 10^(-d)     (matlab: 0/x=0)
%
%                   default d=1, i.e. x/0 -> g=10, 0/x -> g=0.1
%
% calls             nothing

% 12-may-13 ES

% revisions
% 14-apr-20 cleaned up

function g = calc_gain( a, b, d )

% arguments
nargs                   = nargin;
if ( nargs < 2 || isempty( b ) ) && ~isempty( a ) && size( a, 2 ) == 2
    b                   = a( :, 2 );
    a                   = a( :, 1 );
end
if isempty( b )
    g                   = [];
    return
end
if nargs < 3 || isempty( d )
    d                   = 1;
end
d                       = real( round( d( 1 ) ) );

% computations
g                       = a ./ b;
g( a == 0 & b == 0 )    = 1;
g( a == 0 & b ~= 0 )    = 10.^( -d );
g( a ~= 0 & b == 0 )    = 10.^( d );

return

% EOF

