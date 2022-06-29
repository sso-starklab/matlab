% ELLIPSE       plot a tilted elliplse (x,y) with major axis a and minor b.
%
% call          H = ELLIPSE( X, Y, A, B, TILT, C, MC );
%
% get           X, Y        origin
%               A, B        major (X) and minor (Y) axes lengths
%               TILT        {0}; angle (radians)
%               C           {[1 1 1]}; perimiter color
%               MC          {[]}; fill color; if not specified - line is drawn
%
% returns       H           handle to the plot.
%
% calls         nothing

% 28-jan-04 ES

% revisions
% 10-oct-05 TILT added
% 30-dec-19 cleaned up

% reference
% http://mathworld.wolfram.com/Ellipse.html

function h = ellipse( x, y, a, b, tilt, C, MC )

nargs           = nargin;
if nargs < 4
    error( 'usage: ellipse( x, y, a, b )')
end
if nargs < 5 || isempty( tilt )
    tilt        = 0;
end
if nargs < 6 || isempty( C )
    C           = [1 1 1];
end
if nargs < 7 || isempty( MC )
    MC          = []; 
end

if ~isnan(C)
    theta       = ( 0 : 360 ) * pi / 180;
    rot         = [ cos( tilt ) -sin( tilt ); sin( tilt ) cos( tilt ) ];
    xy          = rot * [ a * cos( theta ); b * sin( theta ) ];
    x           = x + xy( 1, : );
    y           = y + xy( 2, : );
    if ~isempty(MC)
        h       = patch( x, y, C );
        set( h, 'EdgeColor', MC );
    else
        h       = line( x, y );
        set( h, 'Color', C );
    end
end

return

% EOF
