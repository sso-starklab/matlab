% CALC_SPATIAL_WAVEFORM_FEATURES        bi-polarity and wavefront
%
% call                  [ vA, vT, bpi ] = calc_spatial_waveform_features( w, TH, graphics )
%
% gets                  w           waveform over multiple channels of a
%                                   single unit, [ nsites x nsamples ]
%
% returns               vA          [uV] vector of ( nsites x 1 ), Amplitude of extremum on each site
%                       vT          [ms] same, Time lag of local extermum relative to global extremum
%                       bpi         scalar describing vA, index of abs(max pos)-abs(max neg) divided by their sum
%
% calls                 alines, plot_spk_waveforms, imupsample, myjet

% 08-dec-19 SSo & ES

% revisions
% 27-dec-19 set vT to NaN for uninformative amplitudes (vA==0)

function [ vA, vT, bpi ] = calc_spatial_waveform_features( w, TH, spkFs, graphics )

nargs                   = nargin;
if nargs < 1 || isempty( w )
    return
end
if nargs < 2 || isempty( TH )
    TH                  = 10;
end
if nargs < 3 || isempty( spkFs )
    spkFs               = 20000;
end
if nargs < 4 || isempty( graphics )
    graphics            = 0;
end

if ~ismatrix( w )
    return
end
[ m, n ]                = size( w );

% get the extermum value for each site (value and sign):
[ ~, extidx ]           = max( abs( w ), [], 2 );
wt                      = w'; 
vA                      = wt( ( 0 : n : n * ( m - 1 ) )' + extidx );

% determine time lag of local extermum relative to global extremum
[ ~, maxidx ]           = max( abs( vA ) );                 % get the origin
vT                      = extidx - extidx( maxidx );        % shift origin
vT                      = vT / spkFs * 1000;                % convert to ms
% clip to remove noise
vA( abs( vA ) < TH )    = 0;
vT( vA == 0 )           = NaN;

% derive an index defined as the difference between the global minimum
% (negative) and the global maximum (positive)
m1                      = min( vA );            % define as the global minimum 
m2                      = max( vA );            % define as the global max
m1                      = min( m1, 0 );              % force m1 to be non-positive
m2                      = max( m2, 0 );              % force m2 to be non-negative
m1                      = abs( m1 );                 % take the absolute value of the negative extremum
m2                      = abs( m2 );                 %                   " positive "
bpi                     = ( m1 - m2 ) ./ ( m1 + m2 );

if graphics
    
    newplot
    
    % raw waveform
    subplot( 1, 3, 1 )
    plot_spk_waveforms( w )
    alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' );
    
    % same, as an upsampled color map:
    USF                 = 50;
    x0                  = ( ( -n / 2 + 1 ) : ( n / 2 ) ) / spkFs * 1000;
    y0                  = 1 : m;
    [ x1, y1, w1 ]      = imupsample( x0, y0, w, USF );
    subplot( 1, 3, 2 )
    imagesc( x1, y1, w1 )
    axis xy
    colormap( myjet )
    set( gca, 'tickdir', 'out', 'box', 'off' );

    % plot vA
    subplot( 3, 3, 3 )
    plot( vA, y0, '.-b' )
    alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'Amplitude [{\mu}V]' )
    title( sprintf( 'BPI=%0.3g', bpi ) )
    
    % plot vT
    subplot( 3, 3, 6 )
    plot( vT, y0, '.-b' )
    alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'Time lag [ms]' )
    
end

return

% EOF

