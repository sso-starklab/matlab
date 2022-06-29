% shirly_eps_analysis       EPS analysis by Shirly
%
% call                      [ res, sst ] = shirly_eps_analysis( datadir, ilevel )
%
% gets                      datadir         directory with multiple *sst files
%                           ilevel          {B} isolation level threshold
%
% returns                   res             structure with fields...
%                           sst             entire database    
% 
% calls                     minmax, replacetok, uhist                       (general)
%                           struct_cat, struct_select                       (structs)
%                           check_cluster_quality, classify_waveform        (spikes)
%                           utest, optclust                                 (stats)

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

function [ res, sst ] = shirly_eps_analysis( datadir, ilevel )

nargs = nargin;
if nargs < 1 || isempty( datadir )
    if ispc
        datadir     = 'G:\dat\EPS';
    elseif ismac
        datadir     = '/Users/eranstark/Documents/da/pos';
    end
end
if nargs < 2 || isempty( ilevel )
    ilevel          = 'B';
end

%------------------------------------------------------------------
% step 1: load the data
%------------------------------------------------------------------
filenames           = dir( [ datadir filesep '*sst' ] );
nfiles              = length( filenames );

for i               = 1 : nfiles
    % load a specific dataset
    fname           = sprintf( '%s/%s', filenames( i ).folder, filenames( i ).name );
    L               = load( fname, '-mat' );
    nunits          = size( L.sst.shankclu, 1 );
    afilebase       = replacetok( L.sst.filebase, '/', '\' );
    if ~isequal( afilebase, L.sst.filebase )
        %fprintf( '%s was generated in a WIN OS!!!\n', L.sst.filebase )
    end
    % expand/rearrange fields to have same number of elements
    [ pathname, filename, extname ] = fileparts( afilebase );
    
    
    % adding a zero before the filename if it has less than 8 charecters
    if length (filename)<8
      str = [ char( ones(1, 1 ) * '0' ) filename] ; 
      filename = str;
     end
    
    if ~ismember( length( filename ), [ 7 8 ] )
        fprintf( '%s filename does not conform to standards\n', afilebase )
        continue
    end
    
    L.sst.filebase  = repmat( filename, [ nunits 1 ] );
    L.sst.Tsec      = repmat( L.sst.Tsec, [ nunits 1 ] );
    if i == 1
        sst         = L.sst;
    else
        [ sst, rc ] = struct_cat( sst, L.sst );
        if sum( rc <= 0 )
            fildnames = fieldnames( sst ); 
            fprintf( '%s: some mismatching fields\n', L.sst.filebase( 1, : ) )
            disp( fildnames( rc <= 0 ) )
        end
        if isnan( rc )
            fprintf( '%s: some missing fields\n', L.sst.filebase( 1, : ) )
        end
    end
    
end

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
% compare positive and negative spikes for each of the features
amp_myus        = [ mean( amp0 ) mean( amp1 ) ];
amp_sds         = [ std( amp0 ) std( amp1 ) ];
amp_pval        = utest( amp0, amp1 );
t2p_myus        = [ mean( t2p0 ) mean( t2p1 ) ];
t2p_sds         = [ std( t2p0 ) std( t2p1 ) ];
t2p_pval        = utest( t2p0, t2p1 );
wid_myus        = [ mean( wid0 ) mean( wid1 ) ];
wid_sds         = [ std( wid0 ) std( wid1 ) ];
wid_pval        = utest( wid0, wid1 );
gsd_myus        = [ mean( gsd0 ) mean( gsd1 ) ];
gsd_sds         = [ std( gsd0 ) std( gsd1 ) ];
gsd_pval        = utest( gsd0, gsd1 );
fwh_myus        = [ mean( fwh0 ) mean( fwh1 ) ];
fwh_sds         = [ std( fwh0 ) std( fwh1 ) ];
fwh_pval        = utest( fwh0, fwh1 );
asy_myus        = [ mean( asy0 ) mean( asy1 ) ];
asy_sds         = [ std( asy0 ) std( asy1 ) ];
asy_pval        = utest( asy0, asy1 );
bst_myus        = [ nanmean( bst0 ) nanmean( bst1 ) ];     % nanmean
bst_sds         = [ nanstd( bst0 ) nanstd( bst1 ) ];       % nanstd
bst_pval        = utest( bst0, bst1 );
acm_myus        = [ mean( acm0 ) mean( acm1 ) ];
acm_sds         = [ std( acm0 ) std( acm1 ) ];
acm_pval        = utest( acm0, acm1 );

% generate a histogram for each 
bords_amp           = minmax( [ amp0; amp1 ] );
edges_amp          = linspace( bords_amp( 1 ), bords_amp( 2 ), nbins + 1 );
hamp0              = histc( amp0, edges_amp );
hamp1              = histc( amp1, edges_amp );
hamp0( end - 1 )   = hamp0( end - 1 ) + hamp0( end );
hamp0( end )       = [];
hamp1( end - 1 )   = hamp1( end - 1 ) + hamp1( end );
hamp1( end )       = [];
bins_amp           = ( edges_amp( 1 : end - 1 ) + edges_amp( 2 : end ) ) / 2;

bords_t2p          = minmax( [ t2p0; t2p1 ] );
edges_t2p          = linspace( bords_t2p( 1 ), bords_t2p( 2 ), nbins + 1 );
ht2p0              = histc( t2p0, edges_t2p );
ht2p1              = histc( t2p1, edges_t2p );
ht2p0( end - 1 )   = ht2p0( end - 1 ) + ht2p0( end );
ht2p0( end )       = [];
ht2p1( end - 1 )   = ht2p1( end - 1 ) + ht2p1( end );
ht2p1( end )       = [];
bins_t2p           = ( edges_t2p( 1 : end - 1 ) + edges_t2p( 2 : end ) ) / 2;

bords_wid          = minmax( [ wid0; wid1 ] );
edges_wid          = linspace( bords_wid( 1 ), bords_wid( 2 ), nbins + 1 );
hwid0              = histc( wid0, edges_wid );
hwid1              = histc( wid1, edges_wid );
hwid0( end - 1 )   = hwid0( end - 1 ) + hwid0( end );
hwid0( end )       = [];
hwid1( end - 1 )   = hwid1( end - 1 ) + hwid1( end );
hwid1( end )       = [];
bins_wid           = ( edges_wid( 1 : end - 1 ) + edges_wid( 2 : end ) ) / 2;

bords_gsd          = minmax( [ gsd0; gsd1 ] );
edges_gsd          = linspace( bords_gsd( 1 ), bords_gsd( 2 ), nbins + 1 );
hgsd0              = histc( gsd0, edges_gsd );
hgsd1              = histc( gsd1, edges_gsd );
hgsd0( end - 1 )   = hgsd0( end - 1 ) + hgsd0( end );
hgsd0( end )       = [];
hgsd1( end - 1 )   = hgsd1( end - 1 ) + hgsd1( end );
hgsd1( end )       = [];
bins_gsd           = ( edges_gsd( 1 : end - 1 ) + edges_gsd( 2 : end ) ) / 2;

bords_fwh          = minmax( [ fwh0; fwh1 ] );
edges_fwh          = linspace( bords_fwh( 1 ), bords_fwh( 2 ), nbins + 1 );
hfwh0              = histc( fwh0, edges_fwh );
hfwh1              = histc( fwh1, edges_fwh );
hfwh0( end - 1 )   = hfwh0( end - 1 ) + hfwh0( end );
hfwh0( end )       = [];
hfwh1( end - 1 )   = hfwh1( end - 1 ) + hfwh1( end );
hfwh1( end )       = [];
bins_fwh           = ( edges_fwh( 1 : end - 1 ) + edges_fwh( 2 : end ) ) / 2;

bords_asy          = minmax( [ asy0; asy1 ] );
edges_asy          = linspace( bords_asy( 1 ), bords_asy( 2 ), nbins + 1 );
hasy0              = histc( asy0, edges_asy );
hasy1              = histc( asy1, edges_asy );
hasy0( end - 1 )   = hasy0( end - 1 ) + hasy0( end );
hasy0( end )       = [];
hasy1( end - 1 )   = hasy1( end - 1 ) + hasy1( end );
hasy1( end )       = [];
bins_asy           = ( edges_asy( 1 : end - 1 ) + edges_asy( 2 : end ) ) / 2;

% add histogram computations and visualizations for Burst index, ACH_COM, .
% visualize
colors          = [ 0 0.7 0; 1 0 1 ];
figure
subplot( 2, 3, 1 )
bhamp0             = stairs( bins_amp, hamp0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhamp1             = stairs( bins_amp, hamp1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bhamp0 bhamp1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g; p=%0.2g', amp_myus( 1 ), amp_myus( 2 ), amp_pval ) );
xlabel( 'Amps [mV]' )

subplot( 2, 3, 2 )
bht2p0             = stairs( bins_t2p, ht2p0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bht2p1             = stairs( bins_t2p, ht2p1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bht2p0 bht2p1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g; p=%0.2g', t2p_myus( 1 ), t2p_myus( 2 ), t2p_pval ) );
xlabel( 'T2P [ms]' )

subplot( 2, 3, 3 )
bhwid0             = stairs( bins_wid, hwid0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhwid1             = stairs( bins_wid, hwid1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bhwid0 bhwid1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.3g/%0.3g; p=%0.2g', wid_myus( 1 ), wid_myus( 2 ), wid_pval ) );
xlabel( 'Width [ms]' )

subplot( 2, 3, 4 )
bhgsd0             = stairs( bins_gsd, hgsd0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhgsd1             = stairs( bins_gsd, hgsd1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bhgsd0 bhgsd1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.3g/%0.3g; p=%0.2g', gsd_myus( 1 ), gsd_myus( 2 ), gsd_pval ) );
xlabel( 'Geo-SD [sites]' )

subplot( 2, 3, 5 )
bhfwh0             = stairs( bins_fwh, hfwh0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhfwh1             = stairs( bins_fwh, hfwh1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bhfwh0 bhfwh1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.3g/%0.3g; p=%0.2g', fwh_myus( 1 ), fwh_myus( 2 ), fwh_pval ) );
xlabel( 'Geo-FWHM [sites]' )

subplot( 2, 3, 6 )
bhasy0             = stairs( bins_asy, hasy0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhasy1             = stairs( bins_asy, hasy1, 'color', colors( 2, : ), 'linewidth', 1 );
legend( [ bhasy0 bhasy1 ], { 'Negative', 'Positive' } )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g; p=%0.2g', asy_myus( 1 ), asy_myus( 2 ), asy_pval ) );
xlabel( 'Asymmetry' )

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

subplot( 2, 2, 1 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 2 ), [ '.' colors( i ) ] )
end
plot( t2p0, fwh0, 'oc' )
plot( t2p1, fwh1, 'ok' )
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
plot( t2p0, fwh0, 'oc' )
plot( t2p1, fwh1, 'ok' )
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
plot( t2p0, fwh0, 'oc' )
plot( t2p1, fwh1, 'ok' )
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )
axis square

subplot( 2, 3, 5 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 3 ), [ '.' colors( i ) ] )
end
plot( t2p0, acm0, 'oc' )
plot( t2p1, acm1, 'ok' )
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'ACH_COM [ms]' )
axis square

subplot( 2, 3, 6 )
hold on, 
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 4 ), [ '.' colors( i ) ] )
end
plot( t2p0, amp0, 'oc' )
plot( t2p1, amp1, 'ok' )
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'AMP [mV]' )
axis square

% 4.3: quantify each of the two populations of the EPS:

% before quantifying the EPS, re-run with only two clusters
xhat                    = x( isposX == 1, didx );
clu                     = optclust( xhat, clustMode, [], [ 1 2 ] );
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

% add a list of the identities of the positive spikes
res.filebase            = sst.filebase( ispos, : );
res.shankclu            = [ sst.shankclu( ispos, : ) sst.pyr( ispos, : ) ];
res.clu                 = clu;

return

% EOF

% 08-dec-19
% develop metrics for 
% (1) bi/polarity over sites
% (2) time lag of extremum over sites

%dbstop in shirly_eps_analysis at 136
% get the database
[ res, sst ] = shirly_eps_analysis;

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
