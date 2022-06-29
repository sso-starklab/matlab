% shirly_eps_analysis       EPS analysis by Shirly
%
% call                      [ res, sst ] = shirly_eps_analysis_23dec19( datadir, ilevel )
%
% gets                      datadir         directory with multiple *sst files
%                           ilevel          {B} isolation level threshold
%
% returns                   res             structure with fields...
%                           sst             entire database    
% 
% calls                     minmax, replacetok, uhist                                                      (general)
%                           struct_cat, struct_select                                                      (structs)
%                           check_cluster_quality, classify_waveform,calc_spatial_waveform_features        (spikes)
%                           utest, optclust                                                                (stats)

% 12-sep-19 ES & SSo

% revisions
% 12-sep-19 added anlysis of more features 
% 13-sep-19 (1) some filesep etc to support multiple OS
%           (2) clustering of data
%           (3) added output
% 14-sep-19 (1) identities of EPS added to output structure
% 15-sep-19 (1) added calculation of temporal properties of the ACH
%           (2) expanded the EPS clustering to 4 features
% 04-dec-19 added a zero before the filename if it has less than 8 charecters [to support mice indeces above 99, i.e. mP101)
% 08-dec-19 built routine for spatial properties, presently at EOF
% 23-dec-19 updated the function with the routine for spatial properties 
% 27-dec-19 (1) DONE; include results from calc_spatial_waveform_features in sst
%           (2) DONE; modify the filebase field to a cell array 
%           (3) include the session_fano_factor results in sst
%           (4) DONE; add histograms for all new features
%           (5) DONE; make histograms for signed amplitude
%           (6) make pie charts for neg/pos; pos clusters
%           (7) DONE; improve the tagging at the end
% 29-dec-19 (1) iso-potential lines for classifiers
%           (2) DONE; byProb added for histograms
%           (3) DONE; added KW test for POS clusters (relevant for cases in which there are >=3 positive subgroups)

function [ res, sst ] = shirly_eps_analysis( datadir, ilevel, byprob )

nargs = nargin;
if nargs < 1 || isempty( datadir )
    if ispc
        datadir     = 'G:\dat\EPS\sst';
        datadirFF   = 'G:\dat\EPS\fanofactor';

    elseif ismac
        datadir     = '/Users/eranstark/Documents/da/pos';
        datadirFF   = '/Users/eranstark/Documents/da/pos/fanofactor';
    end
end
if nargs < 2 || isempty( ilevel )
    ilevel          = 'B';
end
if nargs < 3 || isempty( byprob )
    byprob          = 0;
end

%------------------------------------------------------------------
% step 1: load the data
%------------------------------------------------------------------
% load the sst data
filenames           = dir( [ datadir filesep '*sst' ] );
nfiles              = length( filenames );

for i               = 1 : nfiles
    % load a specific dataset
    fname           = sprintf( '%s/%s', filenames( i ).folder, filenames( i ).name );
    L               = load( fname, '-mat' );
    % prepare for accumulation
    nunits          = size( L.sst.shankclu, 1 );
    afilebase       = replacetok( L.sst.filebase, '/', '\' );
    [ ~, filebase ] = fileparts( afilebase );
    L.sst.filebase  = repmat( { filebase }, [ nunits 1 ] );
    L.sst.Tsec      = repmat( L.sst.Tsec, [ nunits 1 ] );
    % accumulate datasets in a combined structure
    if i == 1
        sst         = L.sst;
    else
        [ sst, rc ] = struct_cat( sst, L.sst );
        if sum( rc <= 0 )
            fildnames = fieldnames( sst ); 
            fprintf( '%s: some mismatching fields\n', L.sst.filebase{ 1 } )
            disp( fildnames( rc <= 0 ) )
        end
        if isnan( rc )
            fprintf( '%s: some missing fields\n', L.sst.filebase{ 1 } )
        end
    end
    
end

% load the fanofactor data
estr                = '_fanofactors';
filenames           = dir( [ datadirFF filesep '*' estr '.mat' ] );
nfiles              = length( filenames );

for i               = 1 : nfiles
    % load a specific dataset
    fname           = sprintf( '%s/%s', filenames( i ).folder, filenames( i ).name );
    L               = load( fname, '-mat' );
    
    % prepare for accumulation
    nunits          = size( L.s.shankclu, 1 );
    [ ~, filebase ] = fileparts( fname );
    filebase        = filebase( 1 : end - length( estr ) );
    L.s.filebase    = repmat( { filebase }, [ nunits 1 ] );
    L.s.nwindows    = repmat( L.s.nwindows, [ nunits 1 ] );
    L.s.winsizes    = repmat( L.s.winsizes, [ nunits 1 ] );
    
    % accumulate datasets in a combined structure
    if i == 1
        s           = L.s;
    else
        [ s, rc ] = struct_cat( s, L.s );
        if sum( rc <= 0 )
            fildnames = fieldnames( s ); 
            fprintf( '%s: some mismatching fields\n', L.s.filebase{ 1 } )
            disp( fildnames( rc <= 0 ) )
        end
        if isnan( rc )
            fprintf( '%s: some missing fields\n', L.s.filebase{ 1 } )
        end
    end
    
end

% merge the fano factor data into the sst structure
% 


%------------------------------------------------------------------
% step 2: select cells and extract interesting features
%------------------------------------------------------------------
% select cells according to ilevel
gidx            = check_cluster_quality( sst, ilevel );
sst             = struct_select( sst, gidx );

% get some statistics
ispos           = sst.extremum > 0;
nunits          = size( ispos, 1 );
nbins           = ceil( nunits / 5 );

% compute ACH metrics
% we have at least 5 possible metrics:
% (1) Burst index (Royer et al., 2012 NN)
brst            = ( sum( sst.ach( 51 : 58, : ), 1 ) ./ sum( sst.ach( 82 : 101, : ), 1 ) )'; 
brst( brst == 0 ) = NaN;
brst            = log10( brst );
sst.brst        = brst;
% (2) ACH COM
% sst.ach_com; 
% (3-5) "David & Michal" based on the CDF of the ACH

% run calc_spatial_waveform_features function on all the units
% organize vectors needed for calc_spatial_waveform_features function
nunits = size( sst.shankclu, 1 );
ivals = NaN * ones( nunits, 1 );
vAs = NaN * ones( nunits, 11 );
vTs = NaN * ones( nunits, 11 );
for i = 1 : nunits
    w = sst.mean{ i };
    [ vA, vT, bpi ] = calc_spatial_waveform_features( w );
    %[ vA, vT, bpi ] = calc_spatial_waveform_features( w, [], [], 1  );
    ivals( i ) = bpi;
    vAs( i, 1 : length( vA ) ) = vA';
    vTs( i, 1 : length( vA ) ) = vT';
end
tV = nanstd( vTs, 0, 2 ); % estimate the probability of multi-compartmental recording by the variance in the peak lag
% shirly - think about addition metrics

sst.vAs = vAs;
sst.vTs = vTs;
sst.bpi = ivals;
sst.tV = tV;

% acquire data for positive/negative spikes separately
amp0            = sst.maxp2p( ~ispos );
amp1            = sst.maxp2p(  ispos );
t2p0            = sst.tp( ~ispos );
t2p1            = sst.tp(  ispos );
wid0            = 1 ./ sst.fmax( ~ispos ) * 1000;
wid1            = 1 ./ sst.fmax(  ispos ) * 1000;
gsd0            = sst.geo_sd( ~ispos );
gsd1            = sst.geo_sd(  ispos );
fwh0            = sst.geo_fwhm( ~ispos );
fwh1            = sst.geo_fwhm(  ispos );
asy0            = sst.asy( ~ispos );
asy1            = sst.asy(  ispos );
pyr0            = sst.pyr( ~ispos );
pyr1            = sst.pyr(  ispos );
bst0            = sst.brst( ~ispos );
bst1            = sst.brst(  ispos );
acm0            = sst.ach_com( ~ispos );
acm1            = sst.ach_com(  ispos );

%------------------------------------------------------------------
% step 3: compare each feature between positive/negative spikes
%------------------------------------------------------------------

nbins           = 60;
colors          = [ 0 0.7 0; 1 0 1 ];
figure


subplot( 3, 3, 1 );
[ amp_myus, amp_sds, amp_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.maxp2p, ispos, nbins, colors, 'Amps [mV]', byprob );
subplot( 3, 3, 2 );
[ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.maxp2p .* ( 2 * ispos - 1 ), ispos, nbins, colors, 'Signed Amp [mV]', byprob );
subplot( 3, 3, 3 );
[ t2p_myus, t2p_sds, t2p_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.tp, ispos, nbins, colors, 'T2P [ms]', byprob );
subplot( 3, 3, 4 );
[ wid_myus, wid_sds, wid_pval, bins_x, hx0, hx1 ] = make_one_hist( 1 ./ sst.fmax * 1000, ispos, nbins, colors, 'Width [ms]', byprob );
subplot( 3, 3, 5 );
[ gsd_myus, gsd_sds, gsd_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.geo_sd, ispos, nbins, colors, 'Geo-SD [sites]', byprob );
subplot( 3, 3, 6 );
[ fwh_myus, fwh_sds, fwh_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.geo_fwhm, ispos, nbins, colors, 'Geo-FWHM [sites]', byprob );


subplot( 3, 3, 7 );
[ asy_myus, asy_sds, asy_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.asy, ispos, nbins, colors, 'Assymetry', byprob );
subplot( 3, 3, 8 );
[ bst_myus, bst_sds, bst_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.brst, ispos, nbins, colors, 'Burst', byprob );


subplot( 3, 3, 7 );
[ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.bpi, ispos, nbins, colors, 'BPI', byprob );
subplot( 3, 3, 8 );
[ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.tV, ispos, nbins, colors, 'tV [ms]', byprob );
subplot( 3, 3, 9 );
[ acm_myus, acm_sds, acm_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.ach_com, ispos, nbins, colors, 'ACH-COM [ms]', byprob );

%------------------------------------------------------------------
% step 4: % partition positive spikes into groups
%------------------------------------------------------------------
% amplitude:                    maxp2p
% duration (peak-to-trough)     tp, 1./fmax
% spatial span                  geo_sd, geo_fwhm
% assymetry                     asy

% some exploratory presentation of the data

% 4.1: initially, look at negative and positive spikes by several pairs of
% selected features, tag by PYR/INT
[ ~, ~, fsep ] = classify_waveform( [ t2p0 wid0 ] );
xx = [ 0.1 0.9 ];
yy = [ 0.6 1.4 ];

figure

subplot( 2, 3, 1 )
fh = fimplicit( fsep, [ xx yy ] ); 
set( fh, 'color', [ 0 0 0 ] );
hold on, 
plot( t2p0( pyr0 == 1 ), wid0( pyr0 == 1 ), '.r', t2p0( pyr0 == 0 ), wid0( pyr0 == 0 ), '.b' )
plot( t2p1( pyr1 == 1 ), wid1( pyr1 == 1 ), 'or', t2p1( pyr1 == 0 ), wid1( pyr1 == 0 ), 'ob' )
set( gca, 'tickdir', 'out', 'box', 'off' );
axis square
xlabel( 'T2P [ms]' )
ylabel( 'Width [ms]' )

subplot( 2, 3, 2 )
hold on, 
plot( t2p0( pyr0 == 1 ), fwh0( pyr0 == 1 ), '.r', t2p0( pyr0 == 0 ), fwh0( pyr0 == 0 ), '.b' )
plot( t2p1( pyr1 == 1 ), fwh1( pyr1 == 1 ), 'or', t2p1( pyr1 == 0 ), fwh1( pyr1 == 0 ), 'ob' )
set( gca, 'tickdir', 'out', 'box', 'off' );
axis square
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )

subplot( 2, 3, 4 )
hold on, 
plot( t2p0( pyr0 == 1 ), amp0( pyr0 == 1 ), '.r', t2p0( pyr0 == 0 ), amp0( pyr0 == 0 ), '.b' )
plot( t2p1( pyr1 == 1 ), amp1( pyr1 == 1 ), 'or', t2p1( pyr1 == 0 ), amp1( pyr1 == 0 ), 'ob' )
set( gca, 'tickdir', 'out', 'box', 'off' );
axis square
xlabel( 'T2P [ms]' )
ylabel( 'Amp [mV]' )

subplot( 2, 3, 5 )
hold on, 
plot( t2p0( pyr0 == 1 ), fwh0( pyr0 == 1 ), '.k', t2p0( pyr0 == 0 ), fwh0( pyr0 == 0 ), '.k' )
plot( t2p1( pyr1 == 1 ), fwh1( pyr1 == 1 ), '.k', t2p1( pyr1 == 0 ), fwh1( pyr1 == 0 ), '.k' )
set( gca, 'tickdir', 'out', 'box', 'off' );
axis square
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )

subplot( 2, 3, 3 )
hold on, 
plot( t2p0( pyr0 == 1 ), bst0( pyr0 == 1 ), '.r', t2p0( pyr0 == 0 ), bst0( pyr0 == 0 ), '.b' )
plot( t2p1( pyr1 == 1 ), bst1( pyr1 == 1 ), 'or', t2p1( pyr1 == 0 ), bst1( pyr1 == 0 ), 'ob' )
set( gca, 'tickdir', 'out', 'box', 'off' );
axis square
xlabel( 'T2P [ms]' )
ylabel( 'Burst index' )

subplot( 2, 3, 6 )
hold on, 
plot( t2p0( pyr0 == 1 ), acm0( pyr0 == 1 ), '.r', t2p0( pyr0 == 0 ), acm0( pyr0 == 0 ), '.b' )
plot( t2p1( pyr1 == 1 ), acm1( pyr1 == 1 ), 'or', t2p1( pyr1 == 0 ), acm1( pyr1 == 0 ), 'ob' )
set( gca, 'tickdir', 'out', 'box', 'off' );
axis square
xlabel( 'T2P [ms]' )
ylabel( 'ACH COM [ms]' )

% 4.2: for the full dataset, cluster into unique groups
% use optimized GMM clustering for the T2P vs. FWHM data

colors                      = 'rgbmkcy';
figure

% clustering arguments
didx                        = 1 : 2;
clustMode                   = 'sic';
feature_names               = { 'T2P', 'FWHM', 'AMP', 'WID', 'GSD', 'ASY', 'BURST', 'ACH_COM' };

% prepare for clustering
x                           = [ [ t2p0; t2p1 ] [ fwh0; fwh1 ] [ amp0; amp1 ] [ wid0; wid1 ] [ gsd0; gsd1 ] [ asy0; asy1 ] [ bst0; bst1 ] [ acm0; acm1 ] ];
isposX                      = false( size( x, 1 ), 1 );
isposX( size( t2p0, 1 ) + 1 : end ) = 1;

%------------------------------------
% all data

xhat = x( :, didx );
fprintf( 'Running optclust for %d samples, %d features...', size( xhat, 1 ), size( xhat, 2 ) )
clu = optclust( xhat, clustMode );
nclu = length( unique( clu ) );
fprintf( 'found %d clusters!\n', nclu )

plot_all_datapoints = 0;
subplot( 2, 2, 1 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 2 ), [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, fwh0, 'oc' )
    plot( t2p1, fwh1, 'ok' )
end
title( sprintf( 'Full dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )

%------------------------------------
% negative spikes
xhat = x( isposX == 0, didx );
fprintf( 'Running optclust for %d samples, %d features...', size( xhat, 1 ), size( xhat, 2 ) )
clu = optclust( xhat, clustMode );
nclu = length( unique( clu ) );
fprintf( 'found %d clusters!\n', nclu )

subplot( 2, 2, 2 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 2 ), [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, fwh0, 'oc' )
    plot( t2p1, fwh1, 'ok' )
end
title( sprintf( 'Neg dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )

%------------------------------------
% positive spikes
didx = [ 1 : 2 8 3 ];
xhat = x( isposX == 1, didx );
fprintf( 'Running optclust for %d samples, %d features...', size( xhat, 1 ), size( xhat, 2 ) )
clu = optclust( xhat, clustMode );
nclu = length( unique( clu ) );
fprintf( 'found %d clusters!\n', nclu )

subplot( 2, 3, 4 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 2 ), [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, fwh0, 'oc' )
    plot( t2p1, fwh1, 'ok' )
end
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )
axis square

subplot( 2, 3, 5 )
hold on,
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 3 ), [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, acm0, 'oc' )
    plot( t2p1, acm1, 'ok' )
end
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'ACH_COM [ms]' )
axis square

subplot( 2, 3, 6 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 4 ), [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, amp0, 'oc' )
    plot( t2p1, amp1, 'ok' )
end
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'AMP [mV]' )
axis square

% 4.3: quantify each of the three populations of the EPS:

% before quantifying the EPS, re-run with only three clusters
xhat                    = x( isposX == 1, didx );
clu                     = optclust( xhat, clustMode, [], [ 1 2 3 ] );
% quantify all features of the EPS according to the two clusters defined by the two selected features
xhat                    = x( isposX == 1, : ); 
[ nclu1, uclu ]         = uhist( clu );
nclu1                   = nclu1';
nclu                    = length( uclu );
nfeatures               = length( feature_names );
features1_myu           = NaN * ones( nclu, nfeatures );
features1_sds           = NaN * ones( nclu, nfeatures );
features1_pval          = NaN * ones( 1, nfeatures );

for i = 1 : nclu
    idx                 = clu == uclu( i );
    features1_myu( i, : ) = mean( xhat( idx, : ) );
    features1_sds( i, : ) = std( xhat( idx, : ) );
end

% conduct u-test between the two "major" populations
if nclu >= 2
    [ ~, sidx ]         = sort( nclu1, 'descend' );
    slct                = uclu( sidx( 1 : 2 ) );
    idx1                = clu == slct( 1 );
    idx2                = clu == slct( 2 );
    for j               = 1 : nfeatures
        features1_pval( j ) = utest( xhat( idx1, j ), xhat( idx2, j ) );
    end
end

% conduct KW test between all populations
if nclu >= 2
    for j               = 1 : nfeatures
        features1_pval_KW( j ) = kruskalwallis( xhat( :, j ), clu, 'off' );
    end
end

% summary - negative vs. positive spikes
res.feature_names       = feature_names;
res.nclu                = [ sum( ispos == 0 ); sum( ispos == 1 ) ];
res.features01_myu      = [ t2p_myus' fwh_myus' amp_myus' wid_myus' gsd_myus' asy_myus' bst_myus' acm_myus' ];
res.features01_sds      = [ t2p_sds' fwh_sds' amp_sds' wid_sds' gsd_sds' asy_sds' bst_sds' acm_sds' ];
res.features01_pval     = [ t2p_pval' fwh_pval' amp_pval' wid_pval' gsd_pval' asy_pval' bst_pval' acm_pval' ];

% summary - internal clustering of positive spikes
res.nclu1               = nclu1;
res.features1_myu       = features1_myu;
res.features1_sds       = features1_sds;
res.features1_pval      = features1_pval;
res.features1_pval_KW   = features1_pval_KW;

% add a list of the identities of the positive spikes
res.filebase            = sst.filebase( ispos, : );
res.shankclu            = [ sst.shankclu( ispos, : ) sst.pyr( ispos, : ) ];
res.clu                 = clu;

% %------------------------------------------------------------------
% % step 5: % % develop metrics for small generator size assumed
% %------------------------------------------------------------------
% % include: 
% % (1) bi/polarity over sites
% % (2) time lag of extremum over sites
% 
% % organize vectors needed for calc_spatial_waveform_features function
% nunits = size( sst.shankclu, 1 );
% ivals = NaN * ones( nunits, 1 );
% vAs = NaN * ones( nunits, 11 );
% vTs = NaN * ones( nunits, 11 );
% 
% % run calc_spatial_waveform_features function on all the units
% for i = 1 : nunits
%     w = sst.mean{ i };
%     [ vA, vT, bpi ] = calc_spatial_waveform_features( w );
%     %[ vA, vT, bpi ] = calc_spatial_waveform_features( w, [], [], 1  );
%     ivals( i ) = bpi;
%     vAs( i, 1 : length( vA ) ) = vA';
%     vTs( i, 1 : length( vA ) ) = vT';
% end
% 
% % summary:
% figure
% % (1) bpis in a histogram:
% binsize = 0.1;
% edges = ( -1 - binsize/2 ) : binsize : ( 1 + binsize / 2 );
% 
% vals = ivals;                       % use all
% vals = ivals( abs( ivals ) ~= 1 );  % ignore the purely pos/neg
% %vals = ivals( abs( ivals ) < 0.8 ); 
% 
% h = histc( vals, edges );
% h( end - 1 ) = h( end - 1 ) + h( end ); h( end ) = [];
% bincenters = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
% subplot( 2, 2, 1 )
% bh = bar( bincenters, h, 1 ); 
% set( bh, 'FaceColor', [ 0 0 1 ], 'EdgeColor', [ 0 0 1 ] );
% set( gca, 'tickdir', 'out', 'box', 'off' );
% 
% % check the tendency of positive deflections to peak after the negative:
% subplot( 2, 2, 2 ), plot( vAs( : ), vTs( : ), '.' )
% set( gca, 'tickdir', 'out', 'box', 'off' );
% lsline
% [ hh, pp ] = calc_spearman( vAs( : ), vTs( : ), 1000 );


return


%----------------------------------------------------
% internal function
%----------------------------------------------------

function [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( x, tidx, nbins, colors, str, byprob )

% acquire data
x0 = x( tidx == 0 );
x1 = x( tidx == 1 );

% compute mean, SDs, and pval
myus = [ nanmean( x0 ) nanmean( x1 ) ];
sds = [ nanstd( x0 ) nanstd( x1 ) ];
pval        = utest( x0, x1);

% generate histogram
bords_x           = minmax( [ x0; x1 ] );
edges_x          = linspace( bords_x( 1 ), bords_x( 2 ), nbins + 1 );
hx0              = histc( x0, edges_x );
hx1              = histc( x1, edges_x );
hx0( end - 1 )   = hx0( end - 1 ) + hx0( end );
hx0( end )       = [];
hx1( end - 1 )   = hx1( end - 1 ) + hx1( end );
hx1( end )       = [];
bins_x           = ( edges_x( 1 : end - 1 ) + edges_x( 2 : end ) ) / 2;

if byprob
    hx0         = hx0 / sum( hx0 );
    hx1         = hx1 / sum( hx1 );
end
       
% plot histogram
newplot
bhx0             = stairs( bins_x, hx0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhx1             = stairs( bins_x, hx1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bhx0 bhx1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g; p=%0.2g', myus( 1 ), myus( 2 ), pval ) );
xlabel( str )

return

%----------------------------------------------------

% EOF

% 08-dec-19
% develop metrics for 
% (1) bi/polarity over sites
% (2) time lag of extremum over sites

%dbstop in shirly_eps_analysis at 138
% get the database
[ res, sst ] = shirly_eps_analysis_23dec19;

% use example unit 3.13 from mP23_16:
uidx = ismember( sst.filebase, '0mP23_16', 'rows' ) & ismember( sst.shankclu, [ 3 13 ], 'rows' );
w0 = sst.mean{ uidx };

[ vA, vT, bpi ] = calc_spatial_waveform_features( w0, [], [], 1 );

% run this over all units:
nunits = size( sst.shankclu, 1 );
ivals = NaN * ones( nunits, 1 );
vAs = NaN * ones( nunits, 11 );
vTs = NaN * ones( nunits, 11 );
for i = 1 : nunits
    w = sst.mean{ i };
    [ vA, vT, bpi ] = calc_spatial_waveform_features( w );
    ivals( i ) = bpi;
    vAs( i, 1 : length( vA ) ) = vA';
    vTs( i, 1 : length( vA ) ) = vT';
end

% summary:
figure
% (1) bpis in a histogram:
binsize = 0.1;
edges = ( -1 - binsize/2 ) : binsize : ( 1 + binsize / 2 );

vals = ivals;                       % use all
vals = ivals( abs( ivals ) ~= 1 );  % ignore the purely pos/neg
%vals = ivals( abs( ivals ) < 0.8 ); 

h = histc( vals, edges );
h( end - 1 ) = h( end - 1 ) + h( end ); h( end ) = [];
bincenters = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
subplot( 2, 2, 1 )
bh = bar( bincenters, h, 1 ); 
set( bh, 'FaceColor', [ 0 0 1 ], 'EdgeColor', [ 0 0 1 ] );
set( gca, 'tickdir', 'out', 'box', 'off' );

% check the tendency of positive deflections to peak after the negative:
subplot( 2, 2, 2 ), plot( vAs( : ), vTs( : ), '.' )
set( gca, 'tickdir', 'out', 'box', 'off' );
lsline
[ hh, pp ] = calc_spearman( vAs( : ), vTs( : ), 1000 );
