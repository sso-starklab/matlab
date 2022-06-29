% PATCH_BAND            plot a line with a patch around it
%
% CALL                  ph = patch_band( x, y, s, c1, c2, sy, nanflag )
% 
% GETS
%   x, y              parameters for line
%   s                 error in y
%   sy                error in x
%   c1                color for line (default: [ 0 0 1 ])
%   c2                color for patch (default: every zero in c1 is brightened to 0.3)
%   nanflag           remove NaNs before plotting
%
% RETURNS 
%   ph                  handles to line and to patch
%
% see also          sleeve

% 03-apr-13 ES

% revisions
% 12-sep-19 cleaned up and expanded help

function ph = patch_band( x, y, s, c1, c2, sy, nanflag )

ph = [];

nargs = nargin;
if nargs < 2 || ~isequal( numel( x ), numel( y ) ) || isempty( x )
    return
else
    x = x( : );
    y = y( : );
end
if isempty( s )
    s = zeros( size( x ) );
else
    s = s( : );
end
if nargs < 4 || isempty( c1 ) || length( c1 ) ~= 3
    c1 = [ 0 0 1 ];
else
    c1 = c1( : )';
end
if nargs < 5 || isempty( c2 ) || length( c2 ) ~= 3
    c2 = c1;
    c2( c2 == 0 ) = 0.3;
else
    c2 = c2( : )';
end
if nargs < 6 || isempty( sy )
    sy = zeros( size( y ) );
else
    sy = sy( : );
end
if nargs < 7 || isempty( nanflag )
    nanflag = 0;
end

if nanflag
    nans = isnan( x ) | isnan( y );
    x( nans ) = [];
    y( nans ) = [];
    s( nans ) = [];
    sy( nans ) = [];
end
if isempty( x )
    return
end

h0 = ishold;
if ~h0
    hold on
end
ph( 2 ) = patch( [ x + sy; flipud( x - sy ) ], [ y + s; flipud( y - s ) ], c2 );
set( ph( 2 ), 'edgecolor', c2 )
ph( 1 ) = line( x, y, 'color', c1, 'linewidth', 2 );
if ~h0
    hold off
end

return

% EOF
