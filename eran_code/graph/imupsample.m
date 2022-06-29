% IMUPSAMPLE            with x,y vectors and z a 2D image
%
% CALL                  [ x1, y1, z1 ] = imupsample( x0, y0, z0, USF, iNaN, method )
%
% GETS                  x0, y0, z0      x- and y- vectors of the data in z0
%                       USF             {10},   up sampling factor
%                       iNaN            {0},    flag: interpolate over NaNs (default: ignore NaNs)
%                       methods         {'im'}  determines interpolation resolution:
%                                                   'im': mean size
%                                                   'i1': by first interval
% 
% DOES                  linear interpolation between z0 values on the grid
%                           defined by x0, y0
%
% RETURNS               x1, y1, z1      upsampled
% 
% usage
% (1) if data are sampled on a linear evenly-spaced grid, this routine just
%       upsamples z0 so plots are easier to read
% (2) if the data are sampled non-evenly/non-linearly, this routine
%       converts the sampling to linear; this is useful for e.g. log-spaced
%       spectrograms
%
% Notes:
% (1) z0: must be length( y0 ) x length( x0 )
% (2) the interpolation is based on the mean intervals of x0/y0 (relevant for
%   non-uniformly spaced data, e.g. log spaced)
% (3) NaNs are are propagated by interp2. This routine (imupsample) circumvents 
%   the propagation by replacing NaNs with the global mean and then replacing 
%   back with NaNs, which is OK for some purposes. 
%   For other behavior, control externally
%
% CALLS                 nothing.

% 17-jul-12 ES

% revisions
% 08-apr-13 NaN support
% 25-nov-13 method changed to mean sampling interval (1st interval still
%               supported)
% 12-aug-19 type cast input arguments to double
% 18-aug-19 cleaned up

function [ x1, y1, z1 ] = imupsample( x0, y0, z0, USF, iNaN, method )

% initialize output
x1                          = [];
y1                          = [];
z1                          = [];

% arguments
nargs                       = nargin;
[ m, n ]                    = size( z0 );
v                           = find( size( x0 ) == n );
if isempty( v )
    z0                      = z0';
    [ m, n ]                = size( z0 );
    v                       = find( size( x0 ) == n );
    if isempty( v )
        return
    end
end
if v == 1
    x0                      = x0';
end
if numel( x0 ) ~= n
    return
end
v = find( size( y0 ) == m );
if isempty( v )
    return
end
if v == 2
    y0                      = y0';
end
if numel( y0 ) ~= m
    return
end
if nargs < 4 || isempty( USF )
    USF                     = 10;
end
x0                          = double( x0 );
y0                          = double( y0 );
z0                          = double( z0 );
if length( x0 ) == 1 || length( y0 ) == 1 || USF == 0
    x1                      = x0;
    y1                      = y0;
    z1                      = z0;
    return
end
if nargs < 5 || isempty( iNaN )
    iNaN                    = 0;
end
if nargs < 6 || isempty( method )
    method                  = 'im';
end

% interpolate
switch method
    case 'i1'                                                   % 1st interval
        dx                  = diff( x0( 1 : 2 ) ) ;
        dy                  = diff( y0( 1 : 2 ) );
    case 'im'                                                   % mean interval
        dx                  = diff( x0( [ 1 end ] ) ) / ( length( x0 ) - 1 );
        dy                  = diff( y0( [ 1 end ] ) ) / ( length( y0 ) - 1 );
end
x1                          = x0( 1 )  : dx / USF : x0( end );
y1                          = ( y0( 1 ) : dy / USF : y0( end ) )';
[ xi, yi ]                  = meshgrid( x1, y1 );
nans                        = isnan( z0 );
if ~sum( nans( : ) )
    z1                      = interp2( x0, y0, z0, xi, yi, 'linear' );
else
    newval                  = nanmean( z0( : ) );               % ideally, NaNs between non-NaNs should be interpolated over
    z0( nans )              = newval;
    z1                      = interp2( x0, y0, z0, xi, yi, 'linear' );
    if ~iNaN
        x0i                 = repmat( x0, [ length( y0 ) 1 ] );
        y0i                 = repmat( y0, [ 1 length( x0 ) ] );
        v                   = find( nans( : ) ).';
        for i               = v
            xlims           = x0i( i ) + diff( x0( 1 : 2 ) ) / 2 * [ -1 1 ];
            ylims           = y0i( i ) + diff( y0( 1 : 2 ) ) / 2 * [ -1 1 ];
            roi             = xi >= xlims( 1 ) & xi <= xlims( 2 ) & yi >= ylims( 1 ) & yi <= ylims( 2 );
            z1( roi )       = NaN;
        end
    end
end

return

% EOF
