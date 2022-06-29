% CALC_SPATIAL_WAVEFORM_FEATURES        bi-polarity and wavefront
%
% call                  [ vA, vT, bpi, tV, vB, pA, nA, pT, nT ] = calc_spatial_waveform_features( w, TH, delta, spkFs, graphics )
%
% gets                  w           waveform over multiple channels of a
%                                   single unit, [ nsites x nsamples ]
%                       TH          {10}; threshold of SNR 
%                       delta       {10}; "significant" difference between amplitudes of different samples
%                       spkFs       {20000}; sampling rate [Hz]
%                       graphics    {0}; plot
%
% returns               vA          [uV] vector of ( nsites x 1 ), signed amplitude of extremum on each site
%                                       with the condition of peak above TH 
%                       vT          [ms] same, Time lag of local extermum relative to global extremum
%                                       with the condition of peak above TH 
%                       bpi         global BPI, the vB of the extremal channel, 
%                                       the closer to 1, the more symmetric ('bipolar')
%                       tV          [ms] SD of vT
%                       vB          vector of ( nsites x 1), channel level BPI, calculated as: 1 - abs( p - n ) / ( p + n )
%                                       with the condition of peak appearing before trough
%                                       the closer to 1, the more symmetric
%                       pA          [uV] vector of ( nsites x 1 ), peak value on each site
%                       nA          [uV] vector of ( nsites x 1 ), trough value on each site
%                       pT          [ms] same, Time lag of peak relative to global extremum
%                       nT          [ms] same, Time lag of trough relative to global extremum
%
% calls                 alines, plot_spk_waveforms, imupsample, myjet, calibration
%                       local_max

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
% 30-mar-21 commented
% 06-apr-21 modified BP definition to rely on local extrema rather than
%               global extrema 
% 07-apr-21 (1) added delta parameter: "significant" difference between amplitudes of different samples
%           (2) added BPI conditions: global maximum is larger by delta from last sample
%                                     global minimum is smaller by delta from first sample

% to do:
% (1) update documentation (help)
% (2) consider TH as a fraction (e.g. -0.5 will take 50% of vT( maxidx ) )
% (3) extend BPI to AUC (compute pA and nA using AUC)
% (4) allow a window for local extremum (samples 7-24), or use local
% maximum followed by local minimum
% (5) fix red and blue dots 
% (6) consider reference of baseline instead of zero

function [ vA, vT, bpi, tV, vB, pA, nA, pT, nT ] = calc_spatial_waveform_features( w, TH, delta, spkFs, graphics )

nargs                           = nargin;
if nargs < 1 || isempty( w )
    return
end
if nargs < 2 || isempty( TH )
    TH                          = 10;                                       % same units as w (e.g. micro-V)
end
if nargs < 3 || isempty( delta )
    delta                       = 10;                                       % same units as w (e.g. micro-V)
end
if nargs < 4 || isempty( spkFs )
    spkFs                       = 20000;
end
if nargs < 5 || isempty( graphics )
    graphics                    = 0;
end

if ~ismatrix( w )
    return
end
[ m, n ]                        = size( w );

% deslope and shift first sample to zero
%w                               = deslope( w', 1 )';

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
% 1 - abs( a - b ) / ( a + b )
vB                              = 1 - abs( pA - abs( nA ) ) ./ ( pA + abs( nA ) );

% assign NaNs to non-classical BPI
% (1) first condition: positive peak before negative trough
% if pT < nT (i.e. positive peak before negative trough)
ridx1                           = pT >= nT;
% vB( pT >= nT )                  = NaN;
% (2) second condition: 
%   local maximum before local minimum (removes maxima at first sample),
%   without additional local extrema (i.e. remove tripolar waveforms)
[ eidx, ~, etype ]              = local_max( w', 'ext' );
uc                              = unique( eidx( :, 2 ) );
kidx                            = false( m, 1 );
for i = 1 : length( uc )
    cidx                        = eidx( :, 2 ) == uc( i );
    kidx( uc( i ) )             = isequal( etype( cidx, : ), [ 1; -1 ] );
end
ridx2                           = ~kidx;
% (3) third condition: 
%   no global maximum at last sample
ridx3                           = pT == ( n - orig ) / spkFs * 1000;
% ridx                            = ridx2 | ridx3;
% (4) forth condition: 
%  global maximum is larger by delta from last sample
wEnd                            = w(:,32);
ridx4                           = (pA - wEnd) < delta;
% (5) fifth condition: 
%  global minimum is smaller by delta from first sample
wStart                          = w(:,1);
ridx5                           = (wStart - nA) < delta;

% ridx                            = ridx2 | ridx4 | ridx5;

ridx                            = ridx2 | ridx3 | ridx5; % this is currently (07-apr-21) the correct one 

vB( ridx )                      = NaN;

% derive a global BPI as the vB of the extremal channel
bpi                             = vB( maxidx );

% derive a global BPI as the weighted mean of the channel-specific index
%bpi                             = calc_com( vB, abs( vA ) );

if graphics
    
    newplot
    
    % raw waveform
    xx                          = ( ( -n / 2 + 1 ) : 1 : n / 2 ) / spkFs * 1000; % [ms]
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

