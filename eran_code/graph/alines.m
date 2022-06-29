% alines        plots seperators on current axis.
%
% lh = alines( x, ax, varargin )
%
% x         points to place separators
% ax        which axis to separate {'x'}
% the rest of the arguments are passed to line
%
% lh        handle to the lines

% 15-jan-13 ES

function lh = alines( x, ax, varargin )

if nargin < 2 || isempty( ax ), ax = 'x'; end
ax = lower( ax( 1 ) );

x = x( : );
nx = length( x );
ah = gca;
if ax == 'x'
    X = [ x x NaN * ones( nx, 1 ) ]';
    X = X( : );
    Y = repmat( [ get( ah, 'ylim' ) NaN ]', nx, 1 );
elseif ax == 'y'
    Y = [ x x NaN * ones( nx, 1) ]';
    Y = Y( : );
    X = repmat( [ get( ah, 'xlim' ) NaN ]', nx, 1 );
else
    return
end
lh = line( X, Y, varargin{:} );

return