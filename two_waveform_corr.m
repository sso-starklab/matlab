% two_waveform_corr         extracts waveforms for two sets of channels during one clu
%
% call                      [ s, dspks1, dspks2 ] = two_waveform_corr( filebase, shankclu1, shankclu2, slct1, slct2 )
% 
% calls
% plot_ss, load_spikes, LoadStims, uniteranges, inranges
% LoadXml, mydetrend, my_xcorr, CCG, ftest
% plot_spk_waveforms, alines, bounds, myjet

% 19-jun-20 AL & ES

% revisions
% 21-jun-20 added xcorr analyses
% 23-jun-20 added plots, significance, and additional metrics
% 24-jun-20 (1) added testing of cc for sync vs. non-sync snippets
%           (2) removed spikes during stims

% to do: 
% change AnatGrps to SpkGrps for chans2
% generalized max
% automatic detection of select
% add number of apikes to plots.
% compute variance at extrimum

function [ s, dspks1, dspks2 ] = two_waveform_corr( filebase, shankclu1, shankclu2, slct1, slct2, varargin )% look at 1.2 and 4.11 in mDS1_08

% win                     = ( -15 : 16 ); % snippet in spkFs
% binSizeSEC              = 0.001;        % CCH bin size [s]
% nBins                   = 50;           % number of CCH bins
nargs                       = nargin;
if nargs < 5 || isempty( filebase ) || isempty( shankclu1 ) || isempty( shankclu2 ) || isempty( slct1 ) || isempty( slct2 ) 
    error( 'Five arguments required' )
end

[ win, binSizeSEC, nBins ...
    , graphics ...
    ]                       = ParseArgPairs(...
    { 'win', 'binSizeSEC', 'nBins' ...
    , 'graphics' ...
    }...
    , { ( -15 : 16 ), 0.001, 50 ...
    , 1 ...
    }...
    , varargin{ : } );


% load all spikes
s                       = load_spikes( filebase );
% remove spikes during stimulation
try
Vals                    = LoadStims( filebase );
uvals                   = uniteranges( Vals( :, 1 : 2 ) ); % combine all the segments
if ~isempty( uvals )
    ridx                = inranges( s.res, uvals ); % remove any spike that is in any segment
    s.res( ridx )       = [];
    s.clu( ridx )       = [];
end
catch
end

% get the spike times of the trigger
u1                      = ismember( s.shankclu( :, 1 : 2 ), shankclu1, 'rows' );
u2                      = ismember( s.shankclu( :, 1 : 2 ), shankclu2, 'rows' );
c1                      = s.map( u1, 1 );
c2                      = s.map( u2, 1 );
st1                     = s.res( s.clu == c1 );
st2                     = s.res( s.clu == c2 );
nspks1                  = length( st1 );
nspks2                  = length( st2 );

% determine which channels were used for 1.2 (the 'triggered' unit)
par                     = LoadXml( filebase );
nchannels               = par.nChannels;
chans2                  = sort( par.SpkGrps( shankclu2( 1 ) ).Channels + 1 );
% note - in mDS1_08 in SpkGrps 4, there is chan 69 instead of 60 (1-based 70 instead of 61)
chans1                  = sort( par.AnatGrps( shankclu1( 1 ) ).Channels + 1 );
nchans1                 = length( chans1 );
nchans2                 = length( chans2 );

scalef                  = 1 / 2^par.nBits * par.VoltageRange / par.Amplification * 1e6; % multiply by this to convert A2D units to microVolts

% map the dat file
datfname                = [ filebase '.dat' ];
a                       = memmapfile( datfname, 'Format', 'int16' );

% get snippets of the channels of triggered (shankclu2) during trigger (shankclu1) spikes
nsamples                = length( win );
spks1                   = zeros( nsamples, nchans1, nspks1, 'single' );
spks2                   = zeros( nsamples, nchans2, nspks1, 'single' );
for i                   = 1 : nspks1
    if ~mod( i, 1000 )
        fprintf( 1, '%d ', i )
    end
    st                  = st1( i );
    idx01               = ( nchannels * ( st - 1 + win ) )' * ones( 1, length( chans1 ) );
    idx02               = ones( length( win ), 1 ) * chans1;
    idx1                = idx01 + idx02;
    idx01               = ( nchannels * ( st - 1 + win ) )' * ones( 1, length( chans2 ) );
    idx02               = ones( length( win ), 1 ) * chans2;
    idx2                = idx01 + idx02;
    spks1( :, :, i )    = a.Data( idx1 );
    spks2( :, :, i )    = a.Data( idx2 );
end
clear a
fprintf( 1, 'done loading\n' )

% % to look at a single spike
% figure
% subplot( 1, 2, 1 ), imagesc( spks1( :, :, i )' ), axis xy
% subplot( 1, 2, 2 ), imagesc( spks2( :, :, i )' ), axis xy
% colormap( myjet )

% now detrend each channel
dspks1                  = zeros( size( spks1 ), 'single' );
for i                   = 1 : nchans1
    mat1                = permute( spks1( :, i, : ), [ 1 3 2 ] ); 
    dmat1               = mydetrend( mat1 );
    dspks1( :, i, : )   = permute( dmat1, [ 1 3 2 ] );
end
dspks2                  = zeros( size( spks2 ), 'single' );
for i                   = 1 : nchans2
    mat2                = permute( spks2( :, i, : ), [ 1 3 2 ] ); 
    dmat2               = mydetrend( mat2 );
    dspks2( :, i, : )   = permute( dmat2, [ 1 3 2 ] );
end

% permute (to conform to plot_spk_waveforms.m conventions) and scale
dspks1                  = permute( dspks1, [ 2 1 3 ] ) * scalef;
dspks2                  = permute( dspks2, [ 2 1 3 ] ) * scalef;

mspks1                  = mean( dspks1, 3 );
sspks1                  = std( dspks1, [], 3 );
mspks2                  = mean( dspks2, 3 );
sspks2                  = std( dspks2, [], 3 );

% normalize all waveform to unit length
mat1                    = permute( dspks1( slct1, :, : ), [ 2 3 1 ] ); 
mat2                    = permute( dspks2( slct2, :, : ), [ 2 3 1 ] ); 
mat1n                   = mat1 ./ ( ones( size( mat1, 1 ), 1 ) * sqrt( sum( mat1.^2 ) ) );
mat2n                   = mat2 ./ ( ones( size( mat2, 1 ), 1 ) * sqrt( sum( mat2.^2 ) ) );
v1n                     = var( mat1n( win == 0, : ) );
v2n                     = var( mat2n( win == 0, : ) );

% compute the xcorr between the mean waveforms at the selected channels

[ cc2vec, lags ]        = my_xcorr( mean( mat2, 2 ), mean( mat1, 2 ), [], -1 );
[ maxval, maxidx ]      = max( cc2vec ); % derive the optimal lag
xls                     = lags / par.SampleRate * 1000; % ms
% compute the xcorr between every pair of simultaneously-recorded snippets
cc2                     = my_xcorr( mat2, mat1, [], -1 ); 
cc                      = cc2( maxidx, : );
% also compute the ACH for each snippet (the cross-correlation at zero lag
% between every snippet and the mean of all snippets) - estimates the variance
ac1vec                  = my_xcorr( mean( mat1, 2 ) * ones( 1, size( mat1, 2 ) ), mat1, [], -1 );
ac2vec                  = my_xcorr( mean( mat2, 2 ) * ones( 1, size( mat2, 2 ) ), mat2, [], -1 );
ac1                     = ac1vec( lags == 0, : );
ac2                     = ac2vec( lags == 0, : );

%----------------------------------------------------------------
% determine, for each spike in st1, whether a spike was detected in st2 or not
res                     = [ st1; st2 ];
clu                     = [ ones( nspks1, 1 ); 2 * ones( nspks2, 1 ) ];
[ ~, sidx ]             = sort( res );
clu                     = clu( sidx );
res                     = res( sidx );
SpikesFs                = par.SampleRate;
BinSize                 = binSizeSEC * SpikesFs;
uClu                    = unique( clu );
[ ccg, t, pairs ]       = CCG( res, clu, BinSize, nBins, SpikesFs, uClu, 'count' );
cch                     = ccg( :, 1, 2 );
ach1                    = ccg( :, 1, 1 );
ach2                    = ccg( :, 2, 2 );
clags                   = t * BinSize / SpikesFs * 1000; % [ms]

% reconstruct the histograms from the pairs:
edges                   = -( nBins + 0.5 ) * BinSize : BinSize : ( nBins + 0.5 ) * BinSize;
that                    = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2 / BinSize;
cmat                    = clu( pairs );
% [ sum( ismember( cmat, [ 2 1 ], 'rows' ) ) sum( ismember( cmat, [ 1 2 ], 'rows' ) ), sum( cch ) ]
% [ sum( ismember( cmat, [ 1 1 ], 'rows' ) ), sum( ach1 ) ]
% [ sum( ismember( cmat, [ 2 2 ], 'rows' ) ), sum( ach2 ) ]

pidx                    = ismember( cmat, [ 1 1 ], 'rows' );
difs                    = diff( res( pairs( pidx, : ) ), [], 2 );
h                       = histc( difs, edges );
h( end - 1 )            = h( end - 1 ) + h( end ); 
h( end )                = [];
ach1hat                 = h;
% isequal( ach1, ach1hat )

pidx                    = ismember( cmat, [ 2 2 ], 'rows' );
difs                    = diff( res( pairs( pidx, : ) ), [], 2 );
h                       = histc( difs, edges );
h( end - 1 )            = h( end - 1 ) + h( end ); 
h( end )                = [];
ach2hat                 = h;
% isequal( ach2, ach2hat )

pidx                    = ismember( cmat, [ 2 1 ], 'rows' );
difs                    = diff( res( pairs( pidx, : ) ), [], 2 );
h                       = histc( difs, edges );
h( end - 1 )            = h( end - 1 ) + h( end ); 
h( end )                = [];
cchhat                  = flipud( h );
% isequal( cch, cchhat )

% now that the association with the pairs is clear, use the pairs to determine which spikes are near-simultaneous
sidx                    = abs( difs ) <= ( BinSize / 2 );
atimes                  = res( pairs( pidx, 2 ) ); % column 2 is for clu=1
sydx                    = ismember( st1, atimes( sidx ) );

%----------------------------------------------------------------
% statistics and summary structure
%[ ~, pval1 ]            = ftest( mat1n( win == 0, : ), mat2n( win == 0, : ) );
pval1 = NaN;
pval2                   = utest( cc( sydx ), cc( ~sydx ) );
pval3                   = utest( ac1, ac2 );
pvals                   = [ pval1 pval2 pval3 ];

clear s
s.filebase              = filebase;
s.shankclu              = [ shankclu1; shankclu2 ];
s.nspks                 = [ nspks1 nspks2 ];
s.s1center              = mat1n( win == 0, : )';
s.s2center              = mat2n( win == 0, : )';
s.ac1                   = ac1( : );
s.ac2                   = ac2( : );
s.cc                    = cc( : );
s.sydx                  = sydx;
s.pvals                 = pvals;

%----------------------------------------------------------------
% graphics
%----------------------------------------------------------------
if ~graphics
    return
end

% plot the units
fig( 1 ) = figure;
plot_ss( filebase, shankclu1 );

fig( 2 ) = figure;
plot_ss( filebase, shankclu2 );


%----------------------------------------------------------------
% plot the mean,SD waveforms for all channels, and every instance of the selected channel
fig( 3 )                = figure;

% trigger unit
subplot( 2, 3, 1 )
%plot_spk_waveforms( dspks1, 0, [], [], 1, [ 0 1 1 ] );
plot_spk_waveforms( cat( 3, mspks1, sspks1 ), 0, [], [], 1, [ 0 1 1 ] );
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 's1: %d.%d (%d spikes)', shankclu1( 1 ), shankclu1( 2 ), nspks1 ) )
xlabel( 'Time [ms]' )
set( gca, 'YTick', [] )

subplot( 2, 3, 2 )
plot( win, mean( mat1, 2 ), 'b' ...
    , win, mspks1( slct1, : )' + sspks1( slct1, : )', '--b' ...
    , win, mspks1( slct1, : )' - sspks1( slct1, : )', '--b'  )
xlim( minmax( win ) )
alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' ); 
alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' ); 
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 's1: %d.%d (chan%d; v=%0.3g)', shankclu1( 1 ), shankclu1( 2 ), slct1, v1n ) )

subplot( 2, 3, 3 )
clim1 = bounds( mat1( : ), 0.01 );
imagesc( win, 1 : nspks1, mat1', clim1( : ).' )
axis xy
title( sprintf( '%0.2g to %0.2g \\muV', clim1( 1 ), clim1( 2 ) ) )
set( gca, 'tickdir', 'out', 'box', 'off' );

% triggered unit
subplot( 2, 3, 4 )
plot_spk_waveforms( cat( 3, mspks2, sspks2 ), 0, [], [], 1, [ 0 1 1 ] );
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 's2: shank%d (%d snippets)', shankclu2( 1 ), nspks1 ) )
xlabel( 'Time [ms]' )
set( gca, 'YTick', [] )

subplot( 2, 3, 5 )
%plot( win, mean( mat2, 2 ) )
plot( win, mean( mat2, 2 ), 'b' ...
    , win, mspks2( slct2, : )' + sspks2( slct2, : )', '--b' ...
    , win, mspks2( slct2, : )' - sspks2( slct2, : )', '--b'  )
xlim( minmax( win ) )
alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' ); 
alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' ); 
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 's2: shank%d (chan%d; v=%0.3g)', shankclu2( 1 ), slct2, v2n ) )
set( gca, 'tickdir', 'out', 'box', 'off' );

subplot( 2, 3, 6 )
clim2                   = bounds( mat2( : ), 0.01 );
imagesc( win, 1 : nspks1, mat2', clim2( : ).' )
axis xy
title( sprintf( '%0.2g to %0.2g \\muV', clim2( 1 ), clim2( 2 ) ) )
set( gca, 'tickdir', 'out', 'box', 'off' );

colormap( myjet )


%----------------------------------------------------------------
% plot CCH and CC statisstics

fig( 4 )                = figure;

subplot( 3, 3, 1 )
bar( clags, cch, 1, 'FaceColor', [ 0 0 0 ], 'EdgeColor', [ 0 0 0 ] )
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Time lag [ms]' )
title( sprintf( 'CCH (%d.%d x %d.%d)', shankclu1( 1 ), shankclu1( 2 ), shankclu2( 1 ), shankclu2( 2 ) ) )

subplot( 3, 3, 2 )
bar( clags, ach1, 1, 'FaceColor', [ 0 0 0 ], 'EdgeColor', [ 0 0 0 ] )
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Time lag [ms]' )
title( sprintf( 'ACH (s1: %d.%d)', shankclu1( 1 ), shankclu1( 2 ) ) )

subplot( 3, 3, 3 )
bar( clags, ach2, 1, 'FaceColor', [ 0 0 0 ], 'EdgeColor', [ 0 0 0 ] )
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Time lag [ms]' )
title( sprintf( 'ACH (s2: %d.%d)', shankclu2( 1 ), shankclu2( 2 ) ) )

subplot( 3, 3, 4 )
plot( xls, cc2vec, '-b' )
xlim( minmax( xls ) )
alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Time lag [ms]' )
ylabel( 'Wavefrom xcorr (Pearson)' )
title( 'mean( s1 ) x mean( s2 )' )
hold on
plot( xls( maxidx ), maxval, '.b' );

subplot( 3, 3, 5 )
bedges                  = -1.02 : 0.04 : 1.02;
binC                    = ( bedges( 1 : end - 1 ) + bedges( 2 : end ) ) / 2;
h1                      = histc( cc( sydx ), bedges )';
h1( end - 1 )           = h1( end - 1 ) + h1( end ); 
h1( end )               = [];
h2                      = histc( cc( ~sydx ), bedges )';
h2( end - 1 )           = h2( end - 1 ) + h2( end ); 
h2( end )               = [];
bar( binC, [ h1 h2 ], 1, 'stacked' )
%hist( cc, 50 )
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'CC' )
ylabel( 'Count' )
title( sprintf( 's1 x s2 (sync: %0.2g; others: %0.2g)', median( cc( sydx ) ), median( cc( ~sydx ) ) ) )


subplot( 3, 3, 7 )
hist( ac1, 50 )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 's1 x mean( s1 ); med=%0.3g', median( ac1 ) ) )

subplot( 3, 3, 8 )
hist( ac2, 50 )
title( sprintf( 's2 x mean( s2 ); med=%0.3g', median( ac2 ) ) )
set( gca, 'tickdir', 'out', 'box', 'off' );

return

% EOF

filebase = filebase_lookup( 'mDS2_07' );
shankclu1 = [2 6];
shankclu2 = [3 2];
slct1 = 11;
slct2 = 16;

[ spks1, spks2, st1, st2 ] = dual_sided_prelim( filebase, shankclu1, shankclu2, slct1, slct2 );
% [ s, dspks1, dspks2 ] = two_waveform_corr( filebase, shankclu2, shankclu1, slct2, slct1 );


