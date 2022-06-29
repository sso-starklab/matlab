% MAKE_CL_ON_BAR        at position X.
%
% call                  lh = make_cl_on_bar( x, ytop, ybot, color, width, linewidth, orientation )
%
% calls                 nothing
%
% see also              barwerror

% 25-apr-13 ES

% 24-may-13 orientation added
% 15-dec-19 cleaned up

function lh = make_cl_on_bar( x, ytop, ybot, color, width, linewidth, orientation )

% argument handling
nargs                   = nargin;
if nargs < 4 || isempty( color )
    color               = [ 0 0 1 ];
end
if nargs < 5 || isempty( width )
    width               = 0.4;
end
if nargs < 6 || isempty( linewidth )
    linewidth           = 1;
end
if nargs < 7 || isempty( orientation )
    orientation         = 'vertical';
end
orientation             = lower( orientation( 1 ) );

% prepare for plotting
x                       = x( : );
ytop                    = ytop( : );
ybot                    = ybot( : );
len                     = [ length( x ), length( ytop ), length( ybot ) ]; 
n                       = max( len );
if length( unique( length( len > 1 ) ) ) > 1
    error( 'mismatch' )
end
if len( 1 ) == 1
    x                   = repmat( x, [ n 1 ] );
end
if len( 2 ) == 1
    ytop                = repmat( ytop, [ n 1 ] );
end
if len( 3 ) == 1
    ybot                = repmat( ybot, [ n 1 ] );
end

% actually plot
lh                      = zeros( n, 1 );
for i                   = 1 : n
    if orientation == 'h'
        lh( i )         = line( [ ytop( i ) * [ 1 1 1 ] ybot( i ) * [ 1 1 1 ] ]...
            , [ x( i ) + [ -1 1 ] * width / 2 [ x( i ) x( i ) ] x( i ) + [ -1 1 ] * width / 2 ] );    
    else
        lh( i )         = line( [ x( i ) + [ -1 1 ] * width / 2 [ x( i ) x( i ) ] x( i ) + [ -1 1 ] * width / 2 ]...
            , [ ytop( i ) * [ 1 1 1 ] ybot( i ) * [ 1 1 1 ] ] );
    end
    set( lh( i ), 'color', color, 'linewidth', linewidth );
end

return

% EOF

