% plot_spk_waveforms        of selected/subset of spikes
% 
% CALL                      plot_spk_waveforms(spk,ntoplot,color,zoomfactor,scalefactor,whattoplot)
% 
% GETS                      spk           3D array (channels x samples x instances )
%                           ntoplot       {100} number of waveforms to plot
%                           color         {[0.7 0.7 1]} color to plot the waveforms
%                           zoomfactor    {1} for magnification
%                           scalefactor   {0.12207} for conversion
%                           whattoplot    {[1 1 0]} logical vector: [ waveforms mean SD ]
%                           Invert        {0}, if 1, then plot Channel 1 at top of panel
% 
% RETURNS                   handle to plot
% 
% DOES                      newplot; channel 1 is plotted at bottom of panel
%
% CALLS                     nothing

% 08-jan-12 ES

% revisons
% 25-jan-12 uninverted plots; clibration bar; support of mean,SD input
% 17-aug-19 cleaned up
% 01-sep-19 added default value for Invert

function lh = plot_spk_waveforms( spk, ntoplot, color, zoomfactor, scalefactor, whattoplot, Invert )

% constants (can move to arguments)
Fs                          = 20000;            % just for scaling the time axis as 

% arguments
nargs                       = nargin;
if nargs < 2 || isempty( ntoplot )
    ntoplot                 = 100; 
end
if nargs < 3 || isempty( color )
    color                   = [ 0.7 0.7 1 ]; 
end
mcolor                      = zeros( 1, 3 );
mcolor( color == max( color ) ) = 1;
if isequal( mcolor, [ 1 1 1 ] )
    mcolor                  = [ 0.5 0.5 0.5 ]; 
end
if nargs < 4 || isempty( zoomfactor )
    zoomfactor              = 1; 
end
if nargs < 5 || isempty( scalefactor )
    scalefactor             = 1 / 2^15 * 4 * 1e6 / 1e3; 
end
if nargs < 6 || isempty( whattoplot )
    whattoplot              = [ 1 1 0 ]; 
end
if ntoplot >= 1
    whattoplot( 1 )         = 1; 
end
if ntoplot == 0
    whattoplot( 1 )         = 0; 
end
if length( whattoplot ) < 3
    whattoplot              = [ whattoplot( : ); zeros( 3 - length( whattoplot ), 1 ) ]; 
end
if nargs < 7 || isempty( Invert )
    Invert = 0;
end

% statistics
if isequal( size( ntoplot ), size( spk ) )
    spk                     = cat( 3, spk, ntoplot );
end
[ nchans, nsamps, nspk ]    = size( spk );
ntoplot                     = min(ntoplot,nspk);
mspk                        = mean( spk, 3 ) * scalefactor;
sspk                        = std( spk, [], 3 ) * scalefactor;
if nspk == 2
    mspk                    = spk( :, :, 1 );
    sspk                    = spk( :, :, 2 );
    ntoplot                 = 1;
    whattoplot              = [ 0 1 1 ];
end
t                           = ( -nsamps / 2 : ( nsamps / 2 - 1 ) ) / Fs * 1e3; % ms
Delta                       = max( ( max( mspk, [], 2 ) - min( mspk, [], 2 ) ) / zoomfactor );

if isempty( spk )
    lh                      = NaN;
    return
end

% plot
newplot
for i = 1 : nchans
    if whattoplot( 1 ) % waveforms
        wf                  = shiftdim(spk(i,:,1 : floor( nspk / ntoplot ) : nspk, 1 ) );
        if Invert
            line( t, wf * scalefactor + ( nchans - i ) * Delta, 'Color', color );
        else
            line( t, wf * scalefactor + i  *Delta, 'Color', color );
        end
    end
    if whattoplot( 2 ) % mean
        if Invert
            lh( i, 1 )      = line( t, mspk( i, : ) + ( nchans - i ) * Delta, 'Color', mcolor, 'LineWidth', 2 );
        else
            lh( i, 1 )      = line( t, mspk( i, : ) + i * Delta, 'Color', mcolor, 'LineWidth', 2 );
        end
    end
    if whattoplot( 3 ) % SD
        if Invert
            lh( i, 2 )      = line(t,mspk(i,:)+sspk(i,:)+(nchans-i)*Delta,'Color',mcolor,'LineWidth',1,'LineStyle', '--' );
            lh( i, 3 )      = line(t,mspk(i,:)-sspk(i,:)+(nchans-i)*Delta,'Color',mcolor,'LineWidth',1,'LineStyle', '--' );
        else
            lh( i, 2 )      = line(t,mspk(i,:)+sspk(i,:)+i*Delta,'Color',mcolor,'LineWidth',1,'LineStyle', '--' );
            lh( i, 3 )      = line(t,mspk(i,:)-sspk(i,:)+i*Delta,'Color',mcolor,'LineWidth',1,'LineStyle', '--' );
        end
    end
end
set( gca, 'XLim', t( [ 1 end ] ) );
if Invert
    set( gca, 'Ylim', [ -Delta * zoomfactor nchans * Delta + Delta * ( zoomfactor - 1 ) ] )
else
    set( gca, 'Ylim', [ 0 nchans * Delta + Delta * ( zoomfactor ) ] )
end

% scale bar
lh( i + 1, 1 )              = line( [ -0.5 -0.5 -0.25 ], [ 75 25 25 ] ); %0.25 ms, 50 uV
set( lh( i + 1, 1 ), 'color', [ 0 0.7 0 ], 'linewidth', 2 )

return

% EOF
