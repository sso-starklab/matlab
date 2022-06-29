% SCALETO       the range 0-1.
%
% call          Y = SCALETO( X, NEWRANGE, GFLAG )
%
% gets          X           matrix to scale
%               NEWRANGE    {[0 1]}
%               GFLAG       {0} scales each column separately,
%                               1 scales the entire matrix together
%               OLDRANGE    {[min(x) max(x)]}, what to scale by
%
% calls                     nothing

% 21-mar-13 ES

% revisions
% 14-oct-21 cleaned up

function y                      = scaleto( x, newrange, g, oldrange )

nargs                           = nargin;
if nargs < 2 || isempty( newrange ) || numel( newrange ) ~= 2 
    newrange                    = [ 0 1 ];
end
newrange                        = sort( newrange );
if nargs < 3 || isempty( g )
    g                           = 0;
end
g                               = logical( g( 1 ) );
if nargs < 4 || isempty( oldrange )
    if g
        oldrange                = [ min( x( : ) ); max( x( : ) ) ];
    else
        oldrange                = [ min( x, [], 1 ); max( x, [], 1 ) ];
    end
end
if isequal( size( oldrange ), [ 1 2 ] )
    oldrange                    = oldrange';
end

% scale to 0-1 (or any other, determined by the oldrange) by global- or column-wise min/max
if g
    y                           = ( x - oldrange( 1 ) ) / diff( oldrange );
else
    d                           = ones( size( x, 1 ), 1 );
    num                         = x - d * oldrange( 1, : );
    den                         = d * diff( oldrange, 1, 1 );
    y                           = num ./ den;
end

% shift to new range
if isequal( newrange, [ 0 1 ] )
    return
end
y                               = newrange( 1 ) + y * diff( newrange );

return

% EOF
