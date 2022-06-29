function [xx, yy, medi] = CDF_for_GWN (x,graphics) 

[ ycdf, xcdf ] = cdfcalc( x );

xx = ( xcdf * ones( 1, 2 ) )'; 
xx = [ xx( 1 ); xx( : ) ];

yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
yy( end ) = [];

medi = xx( find( yy >= 0.5, 1, 'first' ) ); 
if graphics
plot( xx, yy, '-k');
end

 return