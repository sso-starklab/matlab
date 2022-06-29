% spikePhasePlot        plot spike phases
%
% CALL                  spikePhasePlot( hRate, B, phsBins, freqBins, shankclu, spec, totDur, figname, savetype )
% 07-apr-13 ES

% revisions
% 18-aug-19 cleaned up

function spikePhasePlot( hRate, B, phsBins, freqBins, shankclu, spec, totDur, figname, savetype )

% constants
USF                         = 3;
plotmode                    = 'imagesc';
zscaling                    = 'linear';
yscaling                    = 'linear';
nCycles                     = 2;

tColor                      = [ 1 0 1 ];                 % fig title text color
tFW                         = 'normal';                     % title font weight
tFS                         = 12;                           % title font size
tXY                         = [ 0.5 0.985 ];                % title location

xFS                         = 10;                           % text font size
xColor                      = [ 1 1 1 ] * 0;                 % text color
xXY                         = [ 0.5 0.9 ];                  % text location

% argument handling
[ pathname, filename, extname ] = fileparts( figname );
fname                       = [ filename extname ];
shankclu                    = [ shankclu NaN * ones( size( shankclu, 1 ), 4 - size( shankclu, 2 ) ) ];

% plot
if isempty( hRate )
    return
end
hRate                       = cat( 3, hRate, NaN * ones( size( hRate( :, :, 1 ) ) ) );
[ fig, ah ]                 = plotSpectrogram( phsBins, freqBins, hRate, USF, plotmode, zscaling, yscaling, nCycles );
for k                       = 1 : size( shankclu, 1 )
    subplot( ah( k ) )
    ashankclu               = shankclu( k, : );
    title( '' )
    txt1                    = sprintf( '%d.%d(%d;%d)'...
        , ashankclu( 1 ), ashankclu( 2 ), ashankclu( 3 ), ashankclu( 4 ) );
    txt2                    = sprintf( '%0.3g/BL=%0.3g'...
        , max( max( hRate( :, :, k ) ) ), B( k ) );
    th1                     = text( min( xlim ) + diff( xlim ) * 0, min( ylim ) + diff( ylim ) * 1, txt1 );
    th2                     = text( min( xlim ) + diff( xlim ) * 1, min( ylim ) + diff( ylim ) * 0, txt2 );
    set( th1, 'color', xColor, 'fontsize', xFS...
        , 'fontweight', tFW, 'HorizontalAlignment', 'left'...
        , 'VerticalAlignment', 'top' )
    set( th2, 'color', xColor, 'fontsize', xFS...
        , 'fontweight', tFW, 'HorizontalAlignment', 'right'...
        , 'VerticalAlignment', 'bottom' )
    
    axis square
end
tstr                        = replacetok( fname, '\_', '_' );
th                          = textf( tXY( 1 ), tXY( 2 ), tstr ); 
set( th, 'color', tColor, 'fontsize', tFS, 'fontweight', tFW )

subplot( ah( k + 1 ) )
axis square
pos                         = get( ah( k + 1 ), 'position' );
axis( ah( k + 1 ), 'off' )
title( '' )
cla
newax                       = axes( 'position', pos );
if ~isempty( spec )
    lh                      = line( log10( spec( :, 2 ) ), spec( :, 1 ), 'color', [ 0 0 1 ] );
    ylim( [ 0 max( freqBins ) ] )
    axis square
    xlims                   = xlim;
    alines( -2 : 2, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    set( gca, 'xlim', xlims, 'XAxisLocation', 'top', 'tickdir', 'out', 'box', 'off' )
end
txt                         = sprintf( '%d sec', ceil( totDur ) );
th                          = text( min( xlim ) + diff( xlim ) * 0.5, min( ylim ) + diff( ylim ) * 1, txt );
set( th, 'color', xColor, 'fontsize', xFS...
    , 'fontweight', tFW, 'HorizontalAlignment', 'center'...
    , 'VerticalAlignment', 'top' )
    axis( newax, 'off' )

% save figure
if ~isempty( pathname )
    if k >= 4
        subplot( ah( 1 ) )
        xlabel( '' )
        ylabel( '' )
        set( ah( 1 ), 'xticklabel', '' );
    end   
    fig_out( fig, 1, [ figname '.' savetype ], savetype );
    fprintf( '%s: Saved %s.%s\n', upper( mfilename ), figname, savetype )
end

return

% EOF
