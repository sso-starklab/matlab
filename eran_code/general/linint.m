% LININT            linear intepolation of scalar values (rows of the input matrix)
%
% pointer is the first column of the matrix
%
% note - works if the missing values are NaN but then isolated samples are
% not interpolated; to override this behavior replace NaNs with a number

% 03-dec-10 ES

% revisions
% 24-nov-12 modified so that interpolation is NOT done at ends

function x = linint( x, val )
%mat = parse( find( x( :, 1 ) == val ) );
mat = parse( find( sum( x == val, 2 ) >= 1 ) );
for i = 1 : size( mat, 1 )
    n = diff( mat( i, : ) ) + 1;
    if mat( i, 1 ) ~= 1 && mat( i, 2 ) ~= size( x, 1 )
        p0 = x( mat( i, 1 ) - 1, : );
        p1 = x( mat( i, 2 ) + 1, : );
        xm = ones( n, 1 ) * p0 + [ 1 : n ]' * ( p1 - p0 ) / ( n + 1 );
        %disp( round( [ i n p0 abs( [ p0 - p1 ] ) ] ) )
        x( mat( i, 1 ) : mat( i, 2 ), : ) = xm;
    end
    %mi( i, : ) = mean( xm ); % keep for debugging purposes
end
return
