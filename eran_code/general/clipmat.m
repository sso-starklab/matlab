% CLIPMAT       clip array values to between minv and maxv

% 26-mar-11 ES

function mat = clipmat( mat, minv, maxv )

if length( minv ) >= 2
    maxv = minv( 2 );
    minv = minv( 1 );
end

mat( mat < minv ) = minv;
mat( mat > maxv ) = maxv;

return