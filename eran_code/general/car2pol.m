% car2pol       convert data from Cartesian to polar coordinates
%
% call          [ r, phi ] = car2pol( x, y )
%
% gets          x, y            Cartesian
%
% returns       r, phi          polar
% 
% calls         nothing
%
% see also      pol2car

% 11-dec-03 ES

% revisions
% 18-aug-19 updated logic 
% 12-dec-19 cleaned up

function [ r, phi ] = car2pol( x, y )

nargs               = nargin;
if nargs < 2
    error( '2 arguments' )
end
if ~isequal( size( x ), size( y ) )
    error( 'input size mismatch' )
end

r                   = sqrt( x.^2 + y.^2 );
phi                 = mod( atan2( y, x ), 2*pi );

return

% EOF
