% imagescbar        scaled image and bar graph
%
% CALL              [ hh ah ] = imagescbar( x, y, z, usf, xy )
%
% GETS              x/y       coordinates (assumed sorted)
%                   z         data matrix, data in columns
%                   usf       (opional) interpolation factor, 
%                               useful when x/y are not evenly spaced
%
% RETURNS           hh        handles to image, bar
%                   ah        handles to the axes
%
% DOES              splits the current axis 3:1 along the x-axis
%                   on the left plots z with each column scaled 0-1
%                   on the right plots a bar graph of the max values of z
%                   this is useful for visualizing dynamics and absolute values together
% 
% EXAMPLE:
% [ hh ] = imagescbar( x, y, z );
% colormap( flipud( gray  ) ), 
% set( hh( 2 ), 'EdgeColor', [ 0 0 1 ], 'FaceColor', [ 0 0 1 ] )
%
% CALLS             imupsample, scale

% 18-mar-13 ES

% revisions
% 28-oct-14 modified trans so that sizes of x,y,z would always match
% 17-aug-19 cleaned up

function [ hh, ah ] = imagescbar( x, y, z, usf, xy )

% arguments
x                           = unique( x( : ) );
y                           = unique( y( : ) );
m                           = length( x ); 
n                           = length( y );
if nargin < 4 || isempty( usf )
    usf                     = 0;
end
if nargin < 5 || isempty( xy )
    xy                      = 1;
end
if isequal( size( z ), [ m n ] )
    trans                   = 1;
    z                       = z';
elseif isequal( size( z ), [ n m ] )
    trans                   = 0;
else
    % size mismatch
    return
end

% upsample and rescale
[ x1, y1, z1 ]              = imupsample( x, y, z, usf );
if size( z1, 2 ) == m %trans
    maxz1                   = max( z1, [], 2 )';
    z1                      = scale( z1' )';
else
    maxz1                   = max( z1 );
    z1                      = scale( z1 );
end

% create a subplot
ah1                         = gca;
pos                         = get( ah1, 'position' );
f                           = 0.75;
pos1                        = [ pos( [ 1 2 ] ) pos( 3 ) * f pos( 4 ) ];
pos2                        = [ pos( 1 ) + pos( 3 ) * f pos( 2 ) pos( 3 ) * ( 1 - f ) pos( 4 ) ];
set( ah1, 'position', pos1 )
ah2                         = axes( 'position', pos2 );

% plot 
subplot( ah1 )
h1                          = imagesc( x1, y1, z1 ); 
if xy
    axis xy
end
ylims                       = ylim;

subplot( ah2 )
if xy
    h2                      = barh( y1( : ), maxz1( : ) );
else
    h2                      = barh( y1( : ), flipud( maxz1( : ) ) );
end
set( h2, 'FaceColor', [ 0 0 1 ], 'EdgeColor', [ 0 0 1 ] );
ylim( ylims )

set( ah1, 'YAxisLocation', 'left',  'tickdir', 'out', 'box', 'off' )
set( ah2, 'YAxisLocation', 'right', 'tickdir', 'out', 'box', 'off' )
set( ah2, 'yticklabel', [] )

% return these values
ah                          = [ ah1 ah2 ];
hh                          = [ h1 h2 ];

return

% EOF
