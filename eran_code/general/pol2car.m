% pol2car       convert data from polar to Cartesian coordinates
%
% call          [ x, y ] = pol2car( r, ph )
%
% gets          r, phi          polar
%
% returns       x, y            Cartesian
% 
% calls         nothing
%
% see also      car2pol

% 11-dec-03 ES

% revisions
% 12-dec-19 cleaned up

function [ x, y ]   = pol2car( r, ph )

nargs               = nargin;
if nargs < 2
    error( '2 arguments' )
end
if ~isequal( size( r ), size( ph ) )
    error( 'input size mismatch' )
end

x                   = r .* cos( ph ); 
y                   = r .* sin( ph );

return

% EOF
