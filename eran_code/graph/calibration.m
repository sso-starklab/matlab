% calibration    adds x-y calibration bars
%
% calibration   without any arguments will add black calibration bars at 
%               the lower left part of the current axis and label with the 
%               values actually used
%
% [ lh th ] = calibration( len, xy, h, s, params )
%
% len           absolute xy lengths. 
%                   negative to specify fractions (0-1)
%                   defaults to 10% of the current limits
% xy            position of bottom left corner of the calibration bar.
%                   imaginary to specify in frations
%                   defaults to 10% from the bottom left corner of the axes
% h             handle to relevant axis
%                   defualts to current axis
% s             1 x 2 cell array of strings to append to text
%                   -logical zero to omit all text (empty string will be
%                    printed, can then e.g. 
%                    >> set( th( 1 ), 'string', 'xval' )
%                   defaults to {'' '' }
% params        parameter/value pairs, passed as is to line 
%
% lh            handle to lines
% th            1 x 2 handles to text
%
% note          assumes 2D plot

% 24-mar-13 ES

% revisions
% 20-nov-13 dual-format support (absolute + fractions)
% 16-apr-19 made sure len and xy are double
% 18-aug-19 cleaned up

function [ lh, th ] = calibration( len, xy, h, s, varargin )

% constants
loc0                        = [ 0.1 0.1 ]; % default xy in relative coordinates

% arguments
lh                          = [];
th                          = [];
nargs                       = nargin;
if nargs < 1 || isempty( len )
    len                     = -[ 0.1 0.1 ];
end
if length( len( : ) ) ~= 2
    fprintf( '%s: input size mismatch (len)', upper( mfilename ) )
    return
end
len                         = double( len );
if nargs < 2 || isempty( xy )
    xy                      = [];
end
xy                          = double( xy );
if nargs < 3 || isempty( h )
    h                       = gca;
end
if nargs < 4 || isempty( s ) || ~isa( s, 'cell' ) || length( s ) ~= 2
    if ~exist( 's', 'var' ) || ~isequal( s, 0 )
        s                   = { '' '' };
    end
end
if ~isempty( varargin )
    lineparams              = varargin;
else
    lineparams              = { 'color', [ 0 0 0 ], 'linewidth', 2 };
end

% get limits, xy
try
    subplot( h )
catch
    fprintf( '%s: input type mismatch (h)', upper( mfilename ) )
    return
end
xlims                       = xlim;
ylims                       = ylim;
dx                          = diff( xlims );
dy                          = diff( ylims );
if isempty( xy )
    xy                      = [ min( xlims ) + loc0( 1 ) * dx min( ylims ) + loc0( 2 ) * dy ];
end
if length( xy( : ) ) ~= 2
    fprintf( '%s: input size mismatch (xy)', upper( mfilename ) )
    return
end
if all( imag( xy ) > 0 )
    f                       = abs( imag( xy ) );
    xy                      = double( [ min( xlims ) min( ylims ) ] + [ dx dy ] .* f( : ).' );
end

% get len
dxy0                        = [ dx dy ];
dxy( len >= 0 )             = len( len >= 0 );
dxy( len < 0 & len >= -1 )  = dxy0( len < 0 & len >= -1 ) .* abs( len( len < 0 & len >= -1 ) );
dxy( len < -1 )             = NaN;
if sum( isnan( dxy ) )
    fprintf( '%s: input type mismatch (len)', upper( mfilename ) )
    return
end

% draw the line
x                           = xy( 1 ) + [ 0 0 dxy( 1 ) ];
y                           = xy( 2 ) + [ dxy( 2 ) 0 0 ];
lh                          = line( x, y, lineparams{ : } );

% add the text
if isequal( s, 0 )
    str1                    = '';
    str2                    = '';
else
    str1                    = sprintf( '%0.3g %s', dxy( 1 ), s{ 1 } );
    str2                    = sprintf( '%0.3g %s', dxy( 2 ), s{ 2 } );
end
if dxy( 1 ) == 0
    dxy( 1 )                = dx * 0.1;
    str1                    = '';
end
if dxy( 2 ) == 0
    dxy( 2 )                = dy * 0.1;
    str2                    = '';
end
th( 1 )                     = text( xy( 1 ) + dxy( 1 ) / 2, xy( 2 ) - dxy( 2 ) / 2, str1 );
set( th( 1 ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', 0 )
th( 2 )                     = text( xy( 1 ) - dxy( 1 ) / 2, xy( 2 ) + dxy( 2 ) / 2, str2 );
set( th( 2 ), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'rotation', 90 )

return

% EOF
