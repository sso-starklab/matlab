% BARWERROR         plot a bar graph w/ error bars
%
% call              [ bh, eh ] = barwerror( x, y, e, barcolor, barwidth, ebarcolor, ebarlinewidth, ebarwidth, orientation )
%
% gets              x       (optional!)
%                   y       values
%                   e       error, the error bars will be 2e
%
% returns           bh      handle to the bar
%                   eh      handle to the error bar
%
% orientation: { 'vertical' }
%
% calls             make_cl_on_bar

% 17-jul-12 ES

% revisions
% 05-aug-12 vectorized
% 24-apr-13 replaced errorbar w/ make_cl_on_bar
% 15-dec-19 cleaned up

function [ bh, eh ] = barwerror( x, y, e, barcolor, barwidth, ebarcolor, ebarlinewidth, ebarwidth, orientation )

% initialize output
bh                      = [];
eh                      = [];

% arguments
nargs                   = nargin;
if nargs < 2 || isempty( y )
    return
end
x                       = x( : );
y                       = y( : );
n                       = length( y );
if isempty( x )
    x                   = 1 : n;
end
if length( x ) ~= n
    return
end
if nargs < 3 || isempty( e )
    e                   = [];
end
e                       = e( : );
if ~isempty( e )
    if length( e ) ~= n
        return
    end
end
if nargs < 4 || isempty( barcolor ) || numel( barcolor ) ~= 3
    barcolor            = [ 0 0 1 ];
end
if nargs < 5 || isempty( barwidth )
    barwidth            = 0.8;
end
if nargs < 6 || isempty( ebarcolor )
    ebarcolor           = [ 0 0 0 ];
end
if nargs < 7 || isempty( ebarlinewidth )
    ebarlinewidth       = 1;
end
if nargs < 8 || isempty( ebarwidth )
    if n == 1
        ebarwidth       = barwidth / 2;
    else
        ebarwidth       = mean( diff( x ) ) * barwidth / 4;
    end
end
if nargs < 9 || isempty( orientation )
    orientation = 'vertical';
end
orientation             = lower( orientation( 1 ) );

% plot
if orientation == 'h'
    bh                  = barh( x, y, barwidth );
else
    bh                  = bar( x, y, barwidth );
end
set( bh, 'EdgeColor', barcolor, 'FaceColor', barcolor )
if length( x ) > 1
    if orientation == 'h'
        ylim( x( [ 1 end ] )' + ( barwidth + diff( x( 1 : 2 ) ) ) / 2 * [ -1 1 ] )
    else
        xlim( x( [ 1 end ] )' + ( barwidth + diff( x( 1 : 2 ) ) ) / 2 * [ -1 1 ] )
    end
end

if ~isempty( e )
    hstate = ishold;
    if ~hstate
        hold on
    end
    eh = make_cl_on_bar( x, y + e, y - e, ebarcolor, ebarwidth, ebarlinewidth, orientation );
    if ~hstate
        hold off
    end
end

return

% EOF 
