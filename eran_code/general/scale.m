% SCALE                     to the range 0-1.
%
% call                      Y = SCALE( X, GFLAG )
%
% gets                      X       matrix to scale
%                           GFLAG   {0}, scale each column
%                                   1, scaling by global min and max
%                                   [ min max ], scaling by preset min and max
%
% called by                 PARTIAL_CROSS_3D_CONCAT_LAGS
%
% calls                     nothing

% 05-apr-07 ES

% revisions
% 27-oct-13 scaling by zero -> zero (not NaN anymore)
% 10-jul-18 [ min max ] modified

function y = scale( x, g )

if exist( 'g', 'var' )
    lg = length( g );
    if lg == 1
        y = ( x - min( x( : ) ) ) / ( max( x( : ) ) - min( x( : ) ) );
    elseif lg == 2
        %y = ( x - g( 1 ) ) / ( g( 2) - g( 1 ) );
        y = scale( x ) * ( g( 2 ) - g( 1 ) ) + g( 1 );
    end
else
    d = ones( size( x, 1 ), 1 );
    num = x - d * min( x );
    den = d * max( x ) - d * min( x );
    y =  num ./ den;
    y( :, sum( den ) == 0 ) = 0;
end

return

% EOF