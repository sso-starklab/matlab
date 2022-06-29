% DERIVE            numerical difference with same number of samples
%
% CALL              DX = DERIVE( X, DT, DN )
%
% x - data (matrix/vector)
% dt - time base (if ~=1, then returns dx/dt)
% dn - n-th derivative (always 1st order though)

% 03-dec-10 ES

% revisions
% 30-oct-14 dn added
% 29-jan-15 dn default to 1
% 14-oct-19 cleaned up

function dx = derive( x, dt, dn )

nargs           = nargin;
if nargs < 1 || isempty( x )
    dx          = zeros( size( x ) );
    return
end
if nargs < 2 
    dt          = [];
end
if nargs < 3 || isempty( dn )
    dn          = 1;
end

dx              = x( ( dn + 1 ) : end, : ) - x( 1 : ( end - dn ), : );
dx              = [ dx( 1 : dn, : ); dx; dx( end - dn + 1, : ) ]; 
fdx             = filter( ones( 2, 1 ) / 2, 1, dx ); 
dx              = fdx( 2 : end, : );
if ~isempty( dt ) && ( length( dt ) == 1 || isequal( size( dt ), size( dx ) ) )
    dx          = dx / dt;
end

return

% EOF
