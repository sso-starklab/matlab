% makeedges         from a vector of bin centers
%
% y = makeedges( x, mode )
%
% x                 vector of bin centers (use makebins to generate from data)
% mode              can be {'linear'} or 'log'
%
% see also          makebins

% 20-jan-13 ES

function y = makeedges( x, mode )

y = [];
if isempty( x )
    return
end
if ~isvector( x )
    x = x( : );
end
if isrow( x )
    trans = 1;
    x = x( : );
else
    trans = 0;
end
m = length( x );
if sum( diff( x ) < 0 )
    x = sort( x );
end
if exist( 'mode', 'var' ) && strcmpi( mode( 1 : 3 ), 'log' )
    v = [ x( 1 ) / 2; x; x( m ) * 2 ];
    y = geomean( [ v( 1 : end - 1 ) v( 2 : end ) ], 2 );
else
    if m == 1
        v = [ x( 1 ) / 2; x; x( m ) * 2 ];
    else
        v = [ x( 1 ) - abs( x( 2 ) ); x; x( m - 1 ) + abs( x( m ) ) ];
    end
    y = mean( [ v( 1 : end - 1 ) v( 2 : end ) ], 2 );
end
if trans
    y = y';
end
%z = [ y( 1 : end - 1 ) y( 2 : end ) ];

return