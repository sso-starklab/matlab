% calc_asg          compute ASG from CCH

% 14-may-20 ES

% revisions
% 15-may-20 (1) expanded i/o to support various formats (partial)
%           (2) added deconvolution of ACH (partial)
% 17-may-20 (1) removed deconvolution of ACH (temporarily)

% input configurations:
% (1) s2s, n1, n2 (assumes s2s has fields: ccg, shankclu, t, nspks, pred, t_ROI)
% (2) the actual CCH, N1, dt, predictor (or W), t_ROI, (ACH)
% (3) spike trains ... and parameters (dt, W, t_ROI, periods)

% to add:
% (1) convolution (to compute the predictor)
% (2) deconvolution (to correct CCH by the ACH)

% Notes: 
% In spikes2spikes, the ACH are computed from the raw spike trains, whereas
% the CCH are computed from the diluted ("burst-filtered", with a cutoff of
% 5 ms - i.e. all spikes which follow ISIs shorter than 5 ms are removed;
% note that the removal in that routine is very aggressive - it is as is
% done in burstf, which will remove in a vector such as [ 5 10 15 18 22 ]
% all spikes except the first - rather than e.g. the second and fourth). 
%
% Also, in spikes2spikes, for same-shank pairs, I simply ignore the
% zero-lag bin in the computations (remove it), then compute everything
% normally (as if all samples are there). Then add a NaN at the center 
% this is probably the best way to do things

% convType = 'gauss';
% BinSizeMS = 1; halfWidthMS = 50; jitWindowMS = 5; SpikesFs = 20000;
% BinSize                 = round( BinSizeMS * SpikesFs / 1000 ); % ms -> samples
% W                       = ceil( 2 * jitWindowMS / ( BinSize / SpikesFs * 1000 ) + 1 );
% nBins                   = halfWidthMS/BinSizeMS; % number of bins on each side of the CCH
% n = 2 * nBins + 1; 
% % thus BinSize is 20 samples, W is 11 bins, and n is 101 bins. 
% cch0 = cch;
%         cch0( nBins + 1, : )                = [];
%         [ pvalsUpper, pred, pvalsLower ]    = cch_conv( cch0, W, convType );
% 
function [ g1, g2, act, sil, s, fig ] = calc_asg( st1, st2, varargin )

[ BinSizeMS, halfWidthMS, jitWindowMS, SpikesFs ...
    , roiMS, convType, alpha0, deadTimeMS, CutOffMS ...
    , ach, t, periods, cmode, graphics ] = ParseArgPairs( ...
    { 'BinSizeMS', 'halfWidthMS', 'jitWindowMS', 'SpikesFs' ...
    , 'roiMS', 'convType', 'alpha0', 'deadTimeMS', 'CutOffMS' ...
    , 'ach', 't', 'periods', 'cmode', 'graphics' } ...
    , { 1, 50, 5, 20000 ...
    , [ 0 5 ], 'gauss', 0.001, 0.0004, 5 ...
    , [], [], [], 1, 0 }...
    , varargin{ : } );

% calc_asg( s2s, n12, mode, varargin )
% calc_asg( cch, n1, mode, varargin )
% calc_asg( s1, s2, mode, varargin )

% 1. get basic parameters (cch, cchbins, N1, dt, alfa, t_ROI
switch cmode
    case 1 % use the data in s2s (cch, pred, ach)
        s2s                     = st1;
        n12                     = st2;                                          % [ shankclu1; shankclu2 ];
        n1                      = find( s2s.shankclu( :, 1 ) == n12( 1, 1 ) & s2s.shankclu( :, 2 ) == n12( 1, 2 ) );
        n2                      = find( s2s.shankclu( :, 1 ) == n12( 2, 1 ) & s2s.shankclu( :, 2 ) == n12( 2, 2 ) );
        
        t                       = s2s.t;
        nBins                   = ( length( t ) - 1 ) / 2; 
        t_ROI                   = s2s.t_ROI;
        
        cch                     = s2s.ccg( :, n1, n2 );
        ach                     = s2s.ccg( :, n1, n1 );
        N1                      = s2s.nspks( n1 );
    case 2 % use the cch as given (and if given, t and ach)
        cch                     = st1;
        N1                      = st2;
        nBins                   = ( length( cch ) - 1 ) / 2; 
        if isempty( t )                                                         % assume default bin size (e.g. 1 ms)
            t                   = ( -nBins : nBins )' * BinSizeMS;
        else
            BinSizeMS             = diff( t( 1 : 2 ) );                         % assume fixed bin size
        end
        W                       = ceil( 2 * jitWindowMS / BinSizeMS + 1 );
        if isempty( ach )
            ach                 = NaN( length( cch ), 1 );
        end
    case 3 % compute raw cch from spike trains (diluted according to periods)
        % get some basic parameters
        nBins                   = halfWidthMS / BinSizeMS;                      % number of bins on each side of the CCH
        BinSize                 = round( BinSizeMS * SpikesFs / 1000 );         % [samples]
        CutOff                  = round( CutOffMS * SpikesFs / 1000 );          % [samples]
        W                       = ceil( 2 * jitWindowMS / ( BinSize / SpikesFs * 1000 ) + 1 );
        t                       = ( -nBins : nBins )' * BinSizeMS;              % [ms]
        t_ROI                   = t >= roiMS( 1 ) & t <= roiMS( 2 );
        % keep only spikes in periods
        if ~isempty( periods )
        end
        % dilute the spike trains
        isis1                   = diff( st1 ); 
        isis2                   = diff( st2 ); 
        idx1                    = find( isis1 <= CutOff ) + 1; 
        idx2                    = find( isis2 <= CutOff ) + 1; 
        st1( idx1 )             = [];
        st2( idx2 )             = [];
        % compute the CCH
        cch                     = calc_cor( st1, st2, nBins );
        N1                      = length( st1 );
        % compute the ACH
        ach                     = calc_cor( st1, st1, nBins );
        ach( nBins + 1 )        = 0;
end
if nBins ~= round( nBins )
    error( 'check CCH length (should be odd)' )
end
if ~isequal( size( cch ), size( ach ), size( t ) )
    error( 'check correspondence between cch, ach, and t lengths (should be the same)' )
end

dt                              = diff( t( 1 : 2 ) ) / 1000;                    % [s]
n                               = 2 * nBins + 1;

% 3. get the convolution predictor (or compute it)
switch cmode
    case 1 % get the predictor from s2s
        alfa                    = s2s.alpha;
        pred                    = s2s.pred( :, n1, n2 );
    case { 2, 3 } % compute based on cch
        alfa                    = alpha0;
        cch0                    = cch;
        [ ~, pred ]             = cch_conv( cch0, W, convType );
end
        
% 4. compute difference: dcch = cch - pred
dcch                            = cch - pred;

% % 2. correct for ACH (if existing)
% ach1                            = ach;
% ach1( nBins + 1 )               = N1;
% %ach1( nBins + 1 + [ -5 5 ] )    = 5000;
% dcch( nBins + 1 )               = mean( dcch( nBins + [ 0 2 ] ) ); % temporary 
% h1                              = deconvfft( dcch, ach1, nBins ); 
% gcch1                           = h1( : ) / dt;

% 5. compute rate: gcch = dcch / ( dt * n1 )
gcch                            = dcch / ( N1 * dt );

%-------------------------------------------------------------------------
% 7. determine if any bin in the ROI is significant (ROI -5 : 5 ms; signficance, 0.001)
nBonf                           = sum( t_ROI );
t_ROI                           = t_ROI & t > 0; % only causal 
gbUpper                         = poissinv( 1 - alfa / nBonf, max( pred( t_ROI), [], 1 ) );
gbLower                         = poissinv( alfa / nBonf, min( pred( t_ROI ), [], 1 ) );
act                             = false( 1 );
sil                             = false( 1 );
if any( cch( t_ROI ) >= gbUpper )
    act                         = 1;
end
if any( cch( t_ROI ) <= gbLower )
    sil                         = 1;
end

%-------------------------------------------------------------------------
% 9. find the first zero crossings before and after the peak, then 
% compute the intergral (can be negative or positive): sum divided by span
% do these for the extrema in the ROI, regardless of the significance 

% support of positive peak
[ ~, maxidx ]                   = max( cch( t_ROI ) );
ft_ROI                          = find( t_ROI );
pidx                            = ft_ROI( maxidx );
si                              = 1;                                            % first preceding negative value
for i                           = pidx : -1 : 1
    if gcch( i ) < 0
        si                      = i + 1;
        break
    end
end
ei                              = n;                                            % first proceeding negative value
for i                           = pidx : n
    if gcch( i ) < 0
        ei                      = i - 1;
        break
    end
end
sidx                            = si : ei;
g1base                          = sidx;
% compute the integral
a1                              = nansum( gcch( sidx ) );
b1                              = sum( ~isnan( gcch( sidx ) ) ) * dt;
g1                              = a1 * b1;

% support of negative peak
[ ~, minidx ]                   = min( cch( t_ROI ) );
ft_ROI                          = find( t_ROI );
pidx                            = ft_ROI( minidx );
si                              = 1;                                            % first preceding positive value
for i                           = pidx : -1 : 1
    if gcch( i ) > 0
        si                      = i + 1;
        break
    end
end
ei                              = n;                                            % first proceeding positive value
for i                           = pidx : n
    if gcch( i ) > 0
        ei                      = i - 1;
        break
    end
end
sidx                            = si : ei;
g2base                          = sidx;
% compute the integral
a2                              = nansum( gcch( sidx ) );
b2                              = sum( ~isnan( gcch( sidx ) ) ) * dt;
g2                              = a2 * b2;

%-------------------------------------------------------------------------
% organize results in a structure
s.cch                           = cch;
s.cchbins                       = t;
s.ach                           = ach;
s.pred                          = pred;
s.dcch                          = dcch;
s.gcch                          = gcch;
s.g1base                        = g1base;
s.g2base                        = g2base;

%-------------------------------------------------------------------------
% plot if requested
if ~graphics
    fig = NaN;
    return
end
fig = figure;
for spn = 1 : 4
    subplot( 2, 2, spn )
    switch spn
        case 1
            bar( s.cchbins, s.cch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            ylabel( 'Counts/bin' )
        case 2
            bar( s.cchbins, s.cch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            line( s.cchbins, s.pred, 'linewidth', 2, 'color', [ 1 0 0 ] )
            ylabel( 'Counts/bin' )
        case 3
            bar( s.cchbins, s.dcch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            ylabel( 'Counts/bin' )
        case 4
            bar( s.cchbins, s.gcch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
            ylabel( 'Conditional rate [spikes/s]' )
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            alines( s.cchbins( minmax( s.g1base ) )' + [ -0.5 0.5 ], 'x', 'linestyle', '--', 'color', [ 1 0 0 ] );
            try
                [ ~, basename ] = fileparts( s2s.filebase );
                str = sprintf( '%s: %d.%d x %d.%d; ASG=%0.3g' ...
                    , basename, n12( 1, 1 ), n12( 1, 2 ) ...
                    , n12( 2, 1 ), n12( 2, 2 ), g1 );
                title( replacetok( str, '\_', '_' ) )
            end
    end
    set( gca, 'box', 'off', 'tickdir', 'out' )
end

return

% EOF

% example 4 (punit)
shankclu1 = [ 1 6 ];
shankclu2 = [ 2 11 ];
filebase        = filebaseLookup( 'mDL5', -16 );

% example 3:
filebase        = filebaseLookup( 'mA234', -15 );
% same shank PYR-INT (mono)
shankclu1 = [ 3 14 ];
shankclu2 = [ 3 16 ];

% example 2:
filebase = filebaseLookup( 'mF84', -3 );
%filebase =  '/Volumes/slab1/mF84/dat/mF84_03/mF84_03';
shankclu1 = [ 1 4 ]; 
shankclu2 = [ 1 3 ];

% example 1 (same shank PYR-INT (mono))
filebase        = filebaseLookup( 'mC41', -33 );
shankclu1 = [ 3 9 ];
shankclu2 = [ 3 7 ];

% load
load( [ filebase '.sst' ], '-mat', 'sst' );
load( [ filebase '.s2s' ], '-mat', 's2s' )

% prepare for call
n12                         = [ shankclu1; shankclu2 ];

% compute the ASG and plot
[ g1, g2, act, sil, s ] = calc_asg( s2s, n12, 'graphics', 1 );


%------------------------------------------------------------
% now run for a full dataset

% example 3:
filebase = filebaseLookup( 'mF108', -2 );

% load
load( [ filebase '.sst' ], '-mat', 'sst' );
load( [ filebase '.s2s' ], '-mat', 's2s' );
% keep only class B and better
shankclu                            = determine_units( filebase );
kidx1                               = ismember( sst.shankclu, shankclu( :, 1 : 2 ), 'rows' );
kidx2                               = ismember( s2s.shankclu, shankclu( :, 1 : 2 ), 'rows' );
sst                                 = struct_select( sst, kidx1 );
s2s                                 = struct_select( s2s, kidx2 );
s2s.ccg                             = s2s.ccg( :, kidx2, kidx2 );
s2s.pred                            = s2s.pred( :, kidx2, kidx2 );
s2s.pvalsUpper                      = s2s.pvalsUpper( :, kidx2, kidx2 );
s2s.pvalsLower                      = s2s.pvalsLower( :, kidx2, kidx2 );
s2s.hiBins                          = s2s.hiBins( :, kidx2, kidx2 );
s2s.loBins                          = s2s.loBins( :, kidx2, kidx2 );
s2s.gbUpper                         = s2s.gbUpper( kidx2, kidx2 );
s2s.gbLower                         = s2s.gbLower( kidx2, kidx2 );

% compute ASGs
nunits                              = size( s2s.shankclu, 1 );
g1mat                               = NaN( nunits, nunits );
g2mat                               = NaN( nunits, nunits );
amat                                = NaN( nunits, nunits );
smat                                = NaN( nunits, nunits );
for n1                              = 1 : nunits
    for n2                          = 1 : nunits
        if n1 == n2
            continue
        end
        n12                         = s2s.shankclu( [ n1 n2 ], : );
        [ g1, g2, act, sil, s ]     = calc_asg( s2s, n12 );
        amat( n1, n2 )              = act;
        smat( n1, n2 )              = sil;
        g1mat( n1, n2 )             = g1;
        g2mat( n1, n2 )             = g2;
    end
end



%-------------------------------------------------------------------------
% to summarize:
% act, sil                % significant effect (activation/silencing)
% g1, g2                  % effect size (gain, in units of spikes/spike)


g1all                       = g1mat( ~isnan( g1mat ) & g1mat ~= 0 ); % ignore zeros, NaNs
g1act                       = g1mat( amat == 1  & g1mat ~= 0 );
g1nact                      = g1mat( amat == 0  & ~isnan( g1mat ) & g1mat ~= 0 );

[ h1, b1 ]                  = hist( log10( g1all ), 100 );
[ h1nact]                   = hist( log10( g1nact ), b1 );
[ h1act ]                   = hist( log10( g1act ), b1 );

g2all                       = g2mat( ~isnan( g2mat ) & g2mat ~= 0 );
g2sil                       = g2mat( smat == 1 & g2mat ~= 0 );
g2nsil                      = g2mat( smat == 0  & ~isnan( g2mat ) & g2mat ~= 0);

[ h2, b2 ]                  = hist( log10( abs( g2all ) ), 100 );
[ h2nsil]                   = hist( log10( abs( g2nsil ) ), b2 );
[ h2sil ]                   = hist( log10( abs( g2sil ) ), b2 );


figure, 
subplot( 1, 2, 1 )
bh = bar( b1', [ h1nact' h1act' ], 1, 'stacked' );
set( bh( 1 ), 'FaceColor', [ 1 1 1 ] * 0.7, 'EdgeColor', [ 1 1 1 ] * 0.7 )
set( bh( 2 ), 'FaceColor', [ 1 0 0 ] * 1, 'EdgeColor', [ 1 0 0 ] * 1 )
title( sprintf( 'Activation (%d/%d pairs); %0.3g', sum( amat( : ) == 1 ), numel( amat ) - size( amat, 1 ), median( g1act ) ) )
xlabel( 'Gain [spikes/spike]' )
lin2log( 'x', 10 )
set( gca, 'tickdir', 'out', 'box', 'off' )
lh1 = alines( log10( median( g1all ) ), 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
lh2 = alines( log10( median( g1act ) ), 'x', 'linestyle', '--', 'color', [ 1 0 0 ] );

subplot( 1, 2, 2 )
bh = bar( b2', [ h2nsil' h2sil' ], 1, 'stacked' );
set( bh( 1 ), 'FaceColor', [ 1 1 1 ] * 0.7, 'EdgeColor', [ 1 1 1 ] * 0.7 )
set( bh( 2 ), 'FaceColor', [ 0 0 1 ] * 0.7, 'EdgeColor', [ 0 0 1 ] * 0.7 )
title( sprintf( 'Inhibition gain (%d/%d pairs); %0.3g', sum( smat( : ) == 1 ), numel( smat ) - size( smat, 1 ), median( g2sil ) ) )
xlabel( '-Gain [spikes/spike]' )
lin2log( 'x', 10 )
set( gca, 'tickdir', 'out', 'box', 'off' )
lh1 = alines( log10( abs( median( g2all ) ) ), 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
lh2 = alines( log10( abs( median( g2sil ) ) ), 'x', 'linestyle', '--', 'color', [ 0 0 0.7 ] );


