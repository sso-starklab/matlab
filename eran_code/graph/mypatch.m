% MYPATCH       calls patch.
%
% call          PH = MYPATCH( X, Y, C, BS, BAR )
%
% gets          X, Y        vectors
%               C           color triplet {[ 0 0 1 ]}
%               BS          histogram bin size {[]}
%               BAR         flag; disconnects bars {0}
%
% returns       PH          handle to the patch
%
% does          plots Y (patched) as a function of X
%               sets the edge to C
%
% calls         PATCH
%
% see also      SLEEVE

% 14-dec-03 ES

% note
% mypatch trick (adding 0 on yhat) is problematic for plotting
% 'analog' data (e.g. SEM sleeve), use sleeve or
% patch directly & set(ph, 'edgecolor',get(ph,'facecolor'));    

function ph = mypatch( x, y, c, bs, bar )

if nargin < 2
    error( '2 argumetns' )
end
if nargin < 3 || isempty( c )
    c = [ 0 0 1 ];
end
if nargin < 4, bs = []; end
if numel( c ) ~= 3 || any( c < 0 ) || any( c > 1 )
    error( 'c must be a color triplet' )
end
if nargin < 5, bar = []; end
sx = size( x );
sy = size( y );
if any( sx ~= sy ) && ( sx( 1 ) ~= sy( 2 ) || sx( 2 ) ~= sy( 1 ) )
    error( 'input length mismatch' )
end
if prod( sx ) ~= max( sx )
    error( 'vector input required' )
end

x = x(:);
y = y(:);
sx = length( x );

% histogram
if ~isempty( bs )
    if bar
        y = [ y'; y'; zeros( 2, sx ) ];
        x = [ x' - bs/2; x' + bs/2; x' + bs/2; [ x( 2 : end ); x( end ) + bs/2 ]' - bs / 2 ];
    else
        y = [ y'; y' ];
        x = [ x' - bs/2; x' + bs/2 ];
    end
    y = y( : );
    x = x( : );
end

% patch
xhat = [ x(1) - eps; x; x(end)+eps ];
%yhat = [ 0; y; 0 ];
yhat = [ y( 1 ) - eps; y; y( end ) + eps ];
ph = patch( xhat, yhat, c );
set( ph, 'edgecolor', c );

return