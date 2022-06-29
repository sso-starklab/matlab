% textf      text on a figure
%
% th = textf( x, y, str, params )
%
% x, y      in figure coordinates
% str       the string
% params    parameter/value pairs, passed as is to text.m

% 26-feb-13 ES

function [ th, ah ] = textf( x, y, str, varargin )

nargs = nargin;
if nargs < 1 || isempty( x ), x = 0.95; end
if nargs < 2 || isempty( y ), y = 0.5; end
if nargs < 3 || isempty( str ), str = ''; end

if isempty( str ) 
    th = [];
    ah = [];
    return
end

figure( gcf );
ah0 = gca;
ah = axes( 'position', [ x y 0.01 0.01 ] );
axis off;
th = text( 0.5, 0.5, str );
set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' )
subplot( ah0 );

return

% EOF