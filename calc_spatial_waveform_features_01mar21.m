% CALC_SPATIAL_WAVEFORM_FEATURES        bi-polarity and wavefront
%
% call                  [ vA, vT, bpi, tV, vB, pA, nA, pT, nT ] = calc_spatial_waveform_features( w, TH, graphics )
%
% gets                  w           waveform over multiple channels of a
%                                   single unit, [ nsites x nsamples ]
%
% returns               vA          [uV] vector of ( nsites x 1 ), signed amplitude of extremum on each site
%                       vT          [ms] same, Time lag of local extermum relative to global extremum
%                       bpi         global BPI, the vB of the extremal channel, 
%                                       the closer 0, the more bipolar
%                       tV          [ms] SD of vT
%                       vB          vector of ( nsites x 1), channel level BPI, calculated as: 1 - abs( p - n ) / ( p + n )
%                                       with the condition of peak appearing before trough
%                       pA          [uV] vector of ( nsites x 1 ), peak value on each site
%                       nA          [uV] vector of ( nsites x 1 ), trough value on each site
%                       pT          [ms] same, Time lag of peak relative to global extremum
%                       nT          [ms] same, Time lag of trough relative to global extremum
%
% calls                 alines, plot_spk_waveforms, imupsample, myjet, calibration

% 08-dec-19 SSo & ES

% revisions
% 27-dec-19 set vT to NaN for uninformative amplitudes (vA==0)
% 23-feb-21 (1) extended vA and vT to {pA,nA} and {pT,nT}
%           (2) modified graphics to include {pA,nA} and {pT,nT}
%           (3) changed BPI to be channel-based, and include a condition on
%                   the sequence (positive then negative)
%           (4) added channel level BPI (vB)
%           (5) compute tV (SD of vT) - unit level metric
% 01-mar-21 updated documentation

% to do:
% (1) update documentation (help) - DONE
% (2) consider TH as a fraction (e.g. -0.5 will take 50% of vT( maxidx ) )
% (3) extend BPI to AUC (compute pA and nA using AUC)
% (4) allow a window for local extremum (samples 7-24)

function [ vA, vT, bpi, tV, vB, pA, nA, pT, nT ] = calc_spatial_waveform_features( w, TH, spkFs, graphics )

nargs                           = nargin;
if nargs < 1 || isempty( w )
    return
end
if nargs < 2 || isempty( TH )
    TH                          = 10;
end
if nargs < 3 || isempty( spkFs )
    spkFs                       = 20000;
end
if nargs < 4 || isempty( graphics )
    graphics                    = 0;
end

if ~ismatrix( w )
    return
end
[ m, n ]                        = size( w );

% get the peak value for each site (value, sign, and lag):
[ ~, pidx ]                     = max( w, [], 2 );
wt                              = w'; 
pA                              = wt( ( 0 : n : n * ( m - 1 ) )' + pidx );

% get the trough value for each site (value, sign, and lag):
[ ~, nidx ]                     = min( w, [], 2 );
wt                              = w'; 
nA                              = wt( ( 0 : n : n * ( m - 1 ) )' + nidx );

% get the extermum value for each site (value and sign):
[ ~, extidx ]                   = max( abs( w ), [], 2 );
wt                              = w'; 
vA                              = wt( ( 0 : n : n * ( m - 1 ) )' + extidx );

% determine time lag of local exterma relative to global extremum
[ ~, maxidx ]                   = max( abs( vA ) );                         % [channels] get the origin
orig                            = extidx( maxidx );                         % [samples]
vT                              = extidx - orig;                         	% shift origin
vT                              = vT / spkFs * 1000;                        % convert to ms
nT                              = nidx - orig;
pT                              = pidx - orig;
nT                              = nT / spkFs * 1000;
pT                              = pT / spkFs * 1000;

% clip to remove noise
vA( abs( vA ) < TH )            = NaN;
vT( isnan( vA ) )               = NaN;
tV                              = nanstd( vT, 0, 1 );                      % estimate the probability of multi-compartmental recording by the variance in the peak lag

% derive BPI for each channel
% if pT < nT (i.e. positive peak before negative trough)
% 1 - abs( a - b ) / ( a + b )
vB                              = 1 - abs( pA - abs( nA ) ) ./ ( pA + abs( nA ) );
vB( pT >= nT )                  = NaN;

% derive a global BPI as the vB of the extremal channel
bpi                             = vB( maxidx );

% % derive an index defined as the difference between the global minimum
% % (negative) and the global maximum (positive)
% m1                              = min( vA );            % define as the global minimum 
% m2                              = max( vA );            % define as the global max
% m1                              = min( m1, 0 );              % force m1 to be non-positive
% m2                              = max( m2, 0 );              % force m2 to be non-negative
% m1                              = abs( m1 );                 % take the absolute value of the negative extremum
% m2                              = abs( m2 );                 %                   " positive "
% bpi                             = ( m1 - m2 ) ./ ( m1 + m2 );

if graphics
    
    newplot
    
    % raw waveform
    xx                          = ( ( -n / 2 + 1 ) : 1 : n / 2 ) / spkFs * 1000; % [ms]\
    ylims                       = [ min( nA ) max( pA ) ];
    xlims                       = xx( [ 1 n ] );
    for i = 1 : m
        ch                      = m - i + 1;
        %subplot( m, 3, ( i - 1 ) * 3 + 1 )
        axes( 'position', [ 0.13 0.1 * ch 0.21341 0.1 ] )
        plot( xx, w( ch, : ), 'k', nT( ch ), nA( ch ), '.b', pT( ch ), pA( ch ), '.r' )
        set( gca, 'tickdir', 'out', 'box', 'off', 'xlim', xlims, 'ylim', ylims );
        alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
        axis off
        if i == m
            calibration( [ 0.5 100 ], [], [], { 'ms', '\muV' } );
        end
    end
     
%     subplot( 1, 3, 1 )
%     plot_spk_waveforms( w )
%     alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
%     set( gca, 'tickdir', 'out', 'box', 'off' );
    
    % same, as an upsampled color map:
    USF                         = 50;
    x0                          = ( ( -n / 2 + 1 ) : ( n / 2 ) ) / spkFs * 1000;
    y0                          = 1 : m;
    [ x1, y1, w1 ]              = imupsample( x0, y0, w, USF );
    subplot( 1, 3, 2 )
    imagesc( x1, y1, w1 )
    axis xy
    colormap( myjet )
    set( gca, 'tickdir', 'out', 'box', 'off' );

    % plot vA
    subplot( 3, 3, 3 )
    plot( vA, y0, '.-b' )
    ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
    xlim( ylims + [ -1 1 ] * 10 )
    alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'Amplitude [{\mu}V]' )
    title( sprintf( '%0.3g \\muV', vA( maxidx ) ) )
    
    % plot vT
    subplot( 3, 3, 6 )
    plot( vT, y0, '.-b' )
    ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
    xlim( xlims )
    alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlabel( 'Time lag [ms]' )
    title( sprintf( '%0.3g ms, SD = %0.3g', vT( maxidx ), tV ) )

    % plot vB
    subplot( 3, 3, 9 )
    plot( vB, y0, '.-b' )
    ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
    xlim( [ 0 1 ] )
    alines( bpi, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' );
    title( sprintf( 'BPI=%0.3g', bpi ) )
    
end

return

% EOF

