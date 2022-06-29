shirly =  load('rez_for_lamir.mat');

a1 =shirly.rez{1,1}{4,1}.sst1.ach_com;
idx1 =shirly.rez{1,1}{4,1}.ispos;
a2 =shirly.rez{2,1}{4,1}.sst1.ach_com;
idx2 =shirly.rez{2,1}{4,1}.ispos;
a3 =shirly.rez{3,1}{4,1}.sst1.ach_com;
idx3 =shirly.rez{3,1}{4,1}.ispos;

pyridx = shirly.rez{1,1}{4,1}.sst1.pyr;

ach_com_punits1 = a1(idx1);
ach_com_punits2 = a2(idx2);
ach_com_punits3 = a3(idx3);

ach_com_nunits = a1(~idx1);
ach_com_nunitsPYR = a1(~idx1&pyridx);
ach_com_nunitsINT = a1(~idx1&~pyridx);

color1 = [0 163 126]/255;
color2 = [145 47 114]/255;
color3 = [42 101 176]/255;

figdiff = figure;
ylim([0 1.1])

[ y1 x1 ] = cdfcalc( ach_com_punits1);
x                               = [ x1 ]; 
y                               = y1(2:end);
xq                              = x( 1 ) : 0.1 : x( end ); 
yq                              = interp1( x, y, xq );
med00i                          = x( find( y >= 0.5, 1, 'first' ) );
xq00                            = x;
yq00                            = y;

ph                              = plot( xq00, yq00, '-','color',color1);
alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' )
text(med00i,1,sprintf('%0.2g',med00i));
axis square
hold on


[ y1 x1 ] = cdfcalc( ach_com_punits2);
x                               = [ x1 ]; 
y                               = y1(2:end);
xq                              = x( 1 ) : 0.1 : x( end ); 
yq                              = interp1( x, y, xq );
med00i                          = x( find( y >= 0.5, 1, 'first' ) );
xq00                            = x;
yq00                            = y;

ph                              = plot( xq00, yq00, '-','color',color2);
alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
text(med00i,1,sprintf('%0.2g',med00i));


[ y1 x1 ] = cdfcalc( ach_com_punits3);
x                               = [ x1 ]; 
y                               = y1(2:end);
xq                              = x( 1 ) : 0.1 : x( end ); 
yq                              = interp1( x, y, xq );
med00i                          = x( find( y >= 0.5, 1, 'first' ) );
xq00                            = x;
yq00                            = y;

ph                              = plot( xq00, yq00, '-','color',color3);
alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
text(med00i,1,sprintf('%0.2g',med00i));


[ y1 x1 ] = cdfcalc( ach_com_nunitsPYR);
x                               = [ x1 ]; 
y                               = y1(2:end);
xq                              = x( 1 ) : 0.1 : x( end ); 
yq                              = interp1( x, y, xq );
med00i                          = x( find( y >= 0.5, 1, 'first' ) );
xq00                            = x;
yq00                            = y;

ph                              = plot( xq00, yq00, '-');
alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
text(med00i,1,sprintf('%0.2g',med00i));

[ y1 x1 ] = cdfcalc( ach_com_nunitsINT);
x                               = [ x1 ]; 
y                               = y1(2:end);
xq                              = x( 1 ) : 0.1 : x( end ); 
yq                              = interp1( x, y, xq );
med00i                          = x( find( y >= 0.5, 1, 'first' ) );
xq00                            = x;
yq00                            = y;

ph                              = plot( xq00, yq00, '-');
alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
text(med00i,1,sprintf('%0.2g',med00i));

% [ y1 x1 ] = cdfcalc( ach_com_nunits);
% x                               = [ x1 ]; 
% y                               = y1(2:end);
% xq                              = x( 1 ) : 0.1 : x( end ); 
% yq                              = interp1( x, y, xq );
% med00i                          = x( find( y >= 0.5, 1, 'first' ) );
% xq00                            = x;
% yq00                            = y;

% ph                              = plot( xq00, yq00, '-');
% alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
% alines( 0.5, 'y', 'color','k', 'linestyle', '--' );

xlim([10 35])
ylabel('Columative distribution')
xlabel('ACH - COM [ms]')
legend
