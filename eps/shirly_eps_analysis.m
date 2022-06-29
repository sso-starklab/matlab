% shirly_eps_analysis       EPS analysis by Shirly
%
% call                      [ res, sst ] = shirly_eps_analysis( datadir )
%
% gets                      datadir         directory with multiple *sst files
%
% optional arguments (given as name/value pairs)
%
%                           ilevel          {B} isolation level threshold
%                           byprob          {0}
%                           onlygather      {0}
%                           ff_suffix       {'_fanofactors_SWS'}
%
% returns                   res             structure with fields...
%                           sst             entire database (only for ilevel units)
%                           ff              fanofactor structure (all units)
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
% 30-dec-19 (1) DONE; onlygather added
%           (2) partial work on making separatrices for NEG
%           (2) partial work on computing mean waveforms for subtypes POS
% 31-dec-19 (1) upon loading, checked for geo_sd and if missing, add NaN
% 17-may-20 (1) datadir ismac modified
%           (2) datadir,datadirFF ispc modified
% 07-jun-20 (1) changed input format to datadir, then varargin
%           (2) fanofactor output
% 10-jun-20 (1) added depth to accumulated sst
%           (2) datadir,datadirFF isunix modified
% 25-jun-20 (1) move gather to after computing burst, bpi, and wavefront features
% 30-mar-21 (1) added spatial waveform features to sst by call to spikes_stats_supplement_bpi 
%               that calls the new version of calc_spatial_waveform_features
% 27-apr-21 (1) added correction for 17th sample extremum
%           (2) added isbip index for bipolar units and added a figure for feature information
% 19-sep-21 (1) added internal function make_two_hist
%           (2) updated output figure 1
% 15-feb-22  call calc_spatial_waveform with sd per unit, instead of unitary TH
% 24-mar-22 added sst.shankcluOrig

function [ res, sst, ff ] = shirly_eps_analysis( datadir, varargin )

% arguments
nargs                   = nargin;
if nargs < 1 || isempty( datadir )
    if ispc             % shirly (nanna) Windows
        datadir         = 'G:\mice\EPS';
        datadirFF       = 'G:\mice\EPS\fanofactor';

    elseif ismac        % eran (hoor)
        datadir         = '/Users/eranstark/Documents/da/punits';
        datadirFF       = '/Users/eranstark/Documents/da/punits/fanofactor';
    elseif isunix      % shirly (nanna) Linux
        datadir         = '/media/shirly/C22865A128659567/mice/EPS';
        datadirFF       = '/media/shirly/C22865A128659567/mice/EPS/fanofactor';
    end
end

% for depth
 %        datadir             = '/media/shirly/C22865A128659567/mice/EPS/sst/depth_linear';
%         datadir             = '/media/shirly/C22865A128659567/mice/EPS/sst/depth';


[ ilevel, byprob, onlygather, ff_suffix ] = ParseArgPairs(...
    { 'ilevel', 'byprob', 'onlygather', 'ff_suffix' }...
    , { 'B', 0, 0, '_fanofactors_SWS' }...
    , varargin{ : } );

% initialize output
res                     = [];
sst                     = [];
ff                      = [];

%------------------------------------------------------------------
% step 1: load the data
%------------------------------------------------------------------
% load the sst data
filenames               = dir( [ datadir filesep 'sst' filesep '*sst' ] );
nfiles                  = length( filenames );

for i                   = 1 : nfiles
    % load a specific dataset
    fname               = sprintf( '%s/%s', filenames( i ).folder, filenames( i ).name );
    L                   = load( fname, '-mat' );
    % prepare for accumulation
    nunits              = size( L.sst.shankclu, 1 );
    afilebase           = replacetok( L.sst.filebase, '/', '\' );
    [ ~, filebase ]     = fileparts( afilebase );
    L.sst.filebase      = repmat( { filebase }, [ nunits 1 ] );
    L.sst.Tsec          = repmat( L.sst.Tsec, [ nunits 1 ] );
    L.sst.shankcluOrig  = (1:nunits)';
    if ~isfield( L.sst, 'geo_sd' )
        L.sst.geo_sd    = NaN( nunits, 1, 'single' ); 
    end
    if ~isfield( L.sst, 'depth' )
        L.sst.depth     = NaN( nunits, 1, 'single' ); 
    end
    % accumulate datasets in a combined structure
    if i == 1
        sst             = L.sst;
    else
        [ sst, rc ]     = struct_cat( sst, L.sst );
        if sum( rc <= 0 )
            fildnames   = fieldnames( sst ); 
            fprintf( '%s: some mismatching fields\n', L.sst.filebase{ 1 } )
            disp( fildnames( rc <= 0 ) )
        end
        if isnan( rc )
            fprintf( '%s: some missing fields\n', L.sst.filebase{ 1 } )
        end
    end
    
end

% load the fanofactor data
try

    filenames           = dir( [ datadirFF filesep '*' ff_suffix '.mat' ] );
    nfiles              = length( filenames );
    
    for i                   = 1 : nfiles
        % load a specific dataset
        fname               = sprintf( '%s/%s', filenames( i ).folder, filenames( i ).name );
        L                   = load( fname, '-mat' );
        
        % prepare for accumulation
        nunits              = size( L.s.shankclu, 1 );
        [ ~, filebase ]     = fileparts( fname );
        filebase            = filebase( 1 : end - length( ff_suffix ) );
        L.s.filebase        = repmat( { filebase }, [ nunits 1 ] );
        L.s.nwindows        = repmat( L.s.nwindows, [ nunits 1 ] );
        L.s.winsizes        = repmat( L.s.winsizes, [ nunits 1 ] );
        
        % accumulate datasets in a combined structure
        if i == 1
            ff              = L.s;
        else
            [ ff, rc ]      = struct_cat( ff, L.s );
            if sum( rc <= 0 )
                fildnames   = fieldnames( ff );
                fprintf( '%s: some mismatching fields\n', L.s.filebase{ 1 } )
                disp( fildnames( rc <= 0 ) )
            end
            if isnan( rc )
                fprintf( '%s: some missing fields\n', L.s.filebase{ 1 } )
            end
        end
        
    end
    
catch
    
    fprintf( 1, '%s: failed compiling fanofactor structure\n', upper( mfilename ) ) 
    
end


%------------------------------------------------------------------
% step 2: select cells and extract interesting features
%------------------------------------------------------------------
% select cells according to ilevel
gidx            = check_cluster_quality( sst, ilevel );
sst             = struct_select( sst, gidx );

% if onlygather
%     return
% end

% get some statistics
ispos           = sst.extremum > 0;
nunits          = size( ispos, 1 );
nbins           = ceil( nunits / 5 );

% compute ACH metrics
% we have at least 5 possible metrics:
% (1) Burst index (Royer et al., 2012 NN)
brst                    = ( sum( sst.ach( 51 : 58, : ), 1 ) ./ sum( sst.ach( 82 : 101, : ), 1 ) )'; 
brst( brst == 0 | isinf( brst ) )       = NaN;
brst                    = log10( brst );
sst.brst                = brst;
% (2) ACH COM
% sst.ach_com; 
% (3-5) "David & Michal" based on the CDF of the ACH

% run calc_spatial_waveform_features function on all the units

sst                 = spikes_stats_supplement_bpi( sst, sst.sd );

% organize vectors needed for calc_spatial_waveform_features function
% nunits                  = size( sst.shankclu, 1 );
% ivals                   = NaN * ones( nunits, 1 );
% vAs                     = NaN * ones( nunits, 11 );
% vTs                     = NaN * ones( nunits, 11 );
% for i                   = 1 : nunits
%     w                   = sst.mean{ i };
%     [ vA, vT, bpi ]     = calc_spatial_waveform_features( w );
%     %[ vA, vT, bpi ] = calc_spatial_waveform_features( w, [], [], 1  );
%     ivals( i )          = bpi;
%     vAs( i, 1 : length( vA ) ) = vA';
%     vTs( i, 1 : length( vA ) ) = vT';
% end
% tV                      = nanstd( vTs, 0, 2 ); % estimate the probability of multi-compartmental recording by the variance in the peak lag
% % shirly - think about addition metrics
% 
% sst.vAs                 = vAs;
% sst.vTs                 = vTs;
% sst.bpi                 = ivals;
% sst.tV                  = tV;

if onlygather
    return
end

is17 = abs(sst.max (17,:))> abs(sst.max(16,:));
is17 = is17';
for w = 1 : length (sst.filebase)              % if the peak sample is 17 and not 16
    if is17 (w)
        if sst.max(17,w) > 0                   % Punit
            maxw = max(sst.max(17,w));
            sst.extremum (w) = maxw;
        else                                   % Nunit
            minw = min(sst.max(17,w));
            sst.extremum (w) = minw;
        end
    end
end

[mm1] = max(sst.max);
[mm2] = min(sst.max);

amp2z = zeros(length(mm2),1)';
for i = 1:length(mm1)
    if mm1(i)>abs(mm2(i))
        amp2z(i) = mm1(i);
    else
        amp2z(i) = mm2(i);
    end
end

isbip           = ~isnan(sst.bpi);
ispos           = sst.extremum > 0 & ~isbip;

% acquire data for positive/negative spikes separately
amp0            = sst.maxp2p( ~ispos & ~isbip);
amp1            = sst.maxp2p(  ispos );
amp2            = sst.maxp2p(  isbip );
t2p0            = sst.tp( ~ispos & ~isbip);
t2p1            = sst.tp(  ispos );
t2p2            = sst.tp(  isbip );
wid0            = 1 ./ sst.fmax( ~ispos & ~isbip) * 1000;
wid1            = 1 ./ sst.fmax(  ispos ) * 1000;
wid2            = 1 ./ sst.fmax(  isbip ) * 1000;
gsd0            = sst.geo_sd( ~ispos & ~isbip);
gsd1            = sst.geo_sd(  ispos );
gsd2            = sst.geo_sd(  isbip );
fwh0            = sst.geo_fwhm( ~ispos & ~isbip);
fwh1            = sst.geo_fwhm(  ispos );
fwh2            = sst.geo_fwhm(  isbip );
asy0            = sst.asy( ~ispos & ~isbip);
asy1            = sst.asy(  ispos );
asy2            = sst.asy(  isbip );
pyr0            = sst.pyr( ~ispos & ~isbip);
pyr1            = sst.pyr(  ispos );
pyr2            = sst.pyr(  isbip );
bst0            = sst.brst( ~ispos & ~isbip);
bst1            = sst.brst(  ispos );
bst2            = sst.brst(  isbip );
acm0            = sst.ach_com( ~ispos & ~isbip );
acm1            = sst.ach_com(  ispos );
acm2            = sst.ach_com(  isbip );

%------------------------------------------------------------------
% step 3: compare each feature between positive/negative spikes
%------------------------------------------------------------------

nbins           = 60;
colors          = [ 138 86 163; 0 163 126; 42 101 176 ] / 255; % Nunits, Punits, Bipolars
title0          = {'Nunits', 'Punits', 'Bipolars'};
figure  % for Punits units

subplot( 3, 3, 1 );
[ amp_myus, amp_sds, amp_pval, bins_x, hx0, hx1 ] = make_two_hist( sst.maxp2p, ispos, isbip, nbins, colors, 'P2P [mV]', byprob, title0 );
subplot( 3, 3, 2 );
[ myus, sds, pval, bins_x, hx0, hx1 ] = make_two_hist( amp2z', ispos, isbip, nbins, colors, 'Amp [mV]', byprob, title0 );
subplot( 3, 3, 3 );
[ t2p_myus, t2p_sds, t2p_pval, bins_x, hx0, hx1 ] = make_two_hist( sst.tp, ispos, isbip, nbins, colors, 'T2P [ms]', byprob, title0 );
subplot( 3, 3, 4 );
[ wid_myus, wid_sds, wid_pval, bins_x, hx0, hx1 ] = make_two_hist( 1 ./ sst.fmax * 1000, ispos, isbip, nbins, colors, 'Width [ms]', byprob, title0 );
subplot( 3, 3, 5 );
[ gsd_myus, gsd_sds, gsd_pval, bins_x, hx0, hx1 ] = make_two_hist( sst.geo_sd, ispos, isbip, nbins, colors, 'Geo-SD [sites]', byprob, title0 );
subplot( 3, 3, 6 );
[ fwh_myus, fwh_sds, fwh_pval, bins_x, hx0, hx1 ] = make_two_hist( sst.geo_fwhm, ispos, isbip, nbins, colors, 'Geo-FWHM [sites]', byprob, title0 );
subplot( 3, 3, 7 );
[ myus, sds, pval, bins_x, hx0, hx1 ] = make_two_hist( sst.tV, ispos, isbip, nbins, colors, 'SD of Time Lag [ms]', byprob, title0 );
subplot( 3, 3, 8 );
[ bst_myus, bst_sds, bst_pval, bins_x, hx0, hx1 ] = make_two_hist( sst.brst, ispos, isbip, nbins, colors, 'Burst', byprob, title0 );
subplot( 3, 3, 9 );
[ acm_myus, acm_sds, acm_pval, bins_x, hx0, hx1 ] = make_two_hist( sst.ach_com, ispos, isbip, nbins, colors, 'ACH-COM [ms]', byprob, title0 );

% subplot( 3, 3, 1 );
% [ amp_myus, amp_sds, amp_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.maxp2p, ispos, nbins, colors, 'Amps [mV]', byprob, 'Positive' );
% subplot( 3, 3, 2 );
% [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.maxp2p .* ( 2 * ispos - 1 ), ispos, nbins, colors, 'Signed Amp [mV]', byprob, 'Positive' );
% subplot( 3, 3, 3 );
% [ t2p_myus, t2p_sds, t2p_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.tp, ispos, nbins, colors, 'T2P [ms]', byprob, 'Positive' );
% subplot( 3, 3, 4 );
% [ wid_myus, wid_sds, wid_pval, bins_x, hx0, hx1 ] = make_one_hist( 1 ./ sst.fmax * 1000, ispos, nbins, colors, 'Width [ms]', byprob, 'Positive' );
% subplot( 3, 3, 5 );
% [ gsd_myus, gsd_sds, gsd_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.geo_sd, ispos, nbins, colors, 'Geo-SD [sites]', byprob, 'Positive' );
% subplot( 3, 3, 6 );
% [ fwh_myus, fwh_sds, fwh_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.geo_fwhm, ispos, nbins, colors, 'Geo-FWHM [sites]', byprob, 'Positive' );
% subplot( 3, 3, 7 );
% [ asy_myus, asy_sds, asy_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.asy, ispos, nbins, colors, 'Assymetry', byprob, 'Positive' );
% subplot( 3, 3, 8 );
% [ bst_myus, bst_sds, bst_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.brst, ispos, nbins, colors, 'Burst', byprob, 'Positive' );
% 
% 
% % subplot( 3, 3, 7 );
% % [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.bpi, ispos, nbins, colors, 'BPI', byprob );
% subplot( 3, 3, 7 );
% [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.tV, ispos, nbins, colors, 'SD of Time Lag [ms]', byprob, 'Positive' );
% subplot( 3, 3, 9 );
% [ acm_myus, acm_sds, acm_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.ach_com, ispos, nbins, colors, 'ACH-COM [ms]', byprob, 'Positive' );


% figure % for bipolar units
% 
% subplot( 3, 3, 1 );
% [ amp_myus, amp_sds, amp_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.maxp2p, isbip, nbins, colors, 'Amps [mV]', 0, 'Bipolar' );
% subplot( 3, 3, 2 );
% [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.maxp2p .* ( 2 * ispos - 1 ), isbip, nbins, colors, 'Signed Amp [mV]', 0, 'Bipolar' );
% subplot( 3, 3, 3 );
% [ t2p_myus, t2p_sds, t2p_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.tp, isbip, nbins, colors, 'T2P [ms]', 0, 'Bipolar' );
% subplot( 3, 3, 4 );
% [ wid_myus, wid_sds, wid_pval, bins_x, hx0, hx1 ] = make_one_hist( 1 ./ sst.fmax * 1000, isbip, nbins, colors, 'Width [ms]', 0, 'Bipolar' );
% subplot( 3, 3, 5 );
% [ gsd_myus, gsd_sds, gsd_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.geo_sd, isbip, nbins, colors, 'Geo-SD [sites]', 0, 'Bipolar' );
% subplot( 3, 3, 6 );
% [ fwh_myus, fwh_sds, fwh_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.geo_fwhm, isbip, nbins, colors, 'Geo-FWHM [sites]', 0, 'Bipolar' );
% 
% 
% subplot( 3, 3, 7 );
% [ asy_myus, asy_sds, asy_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.asy, isbip, nbins, colors, 'Assymetry', byprob, 'Bipolar' );
% subplot( 3, 3, 8 );
% [ bst_myus, bst_sds, bst_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.brst, isbip, nbins, colors, 'Burst', byprob, 'Bipolar' );
% 
% 
% subplot( 3, 3, 7 );
% [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.bpi, isbip, nbins, colors, 'BPI', byprob, 'Bipolar' );
% subplot( 3, 3, 8 );
% [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( sst.tV, isbip, nbins, colors, 'SD of Time Lag [ms]', byprob, 'Bipolar' );
% subplot( 3, 3, 9 );
% [ acm_myus, acm_sds, acm_pval, bins_x, hx0, hx1 ] = make_one_hist( sst.ach_com, isbip, nbins, colors, 'ACH-COM [ms]', byprob, 'Bipolar' );
% 

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

colors                      = 'rgbmkcyrgbmkcy';
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
clu = optclust( xhat, clustMode );%, [], 10 );
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

% prepare for adding iso-probability ellipses
sstN = struct_select( sst, ~ispos );
labels = sstN.pyr;
param1 = sstN.tp; 

subplot( 2, 3, 4 )
hold on, 
sf = 1/10;
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 2 ) * sf, [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, fwh0, 'oc' )
    plot( t2p1, fwh1, 'ok' )
end
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'FWHM [sites]' )
axis square
add_iso_ellipses( param1, sstN.geo_fwhm * sf, labels );            % add ellipses
% scale back the ticks of the y-axis
lbls = get( gca, 'YTickLabel' );
yticks = zeros( 1, length( lbls ) );
for i = 1 : length( lbls )
    yticks( i ) = str2num( lbls{ i } ) / sf;
end
ytick = get( gca, 'ytick' ); 
set( gca, 'ytick', ytick )
set( gca, 'YTickLabel', yticks );

subplot( 2, 3, 5 )
hold on,
sf = 1/50;
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 3 ) * sf, [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, acm0, 'oc' )
    plot( t2p1, acm1, 'ok' )
end
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'ACH_COM [ms]' )
axis square
add_iso_ellipses( param1, sstN.ach_com * sf, labels );            % add ellipses
% scale back the ticks of the y-axis
lbls = get( gca, 'YTickLabel' );
yticks = zeros( 1, length( lbls ) );
for i = 1 : length( lbls )
    yticks( i ) = str2num( lbls{ i } ) / sf;
end
ytick = get( gca, 'ytick' ); 
set( gca, 'ytick', ytick )
set( gca, 'YTickLabel', yticks );

subplot( 2, 3, 6 )
hold on, 
sf = 1;
for i = 1 : nclu
    plot( xhat( clu == i, 1 ), xhat( clu == i, 4 ) * sf, [ '.' colors( i ) ] )
end
if plot_all_datapoints
    plot( t2p0, amp0, 'oc' )
    plot( t2p1, amp1, 'ok' )
end
title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
xlabel( 'T2P [ms]' )
ylabel( 'AMP [mV]' )
axis square
add_iso_ellipses( param1, sstN.maxp2p * sf, labels );            % add ellipses

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

%--------
% 30-dec-19, first part of computing mean waveform per cluster:
%--------
% sstP = struct_select( sst, ispos );
% figure, for i = 1 : 3, subplot( 2, 2, i ), plot( ( -15 : 16 ) / 20, sstP.max( :, clu == i ), 'b' ), end
% figure, for i = 1 : 3, subplot( 2, 2, i ), imagesc( scale( sstP.max( :, clu == i ) ), end


return


%----------------------------------------------------
% internal functions
%----------------------------------------------------

%----------------------------------------------------
% make_one_hist

function [ myus, sds, pval, bins_x, hx0, hx1 ] = make_one_hist( x, tidx, nbins, colors, str, byprob, title0 )

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
legend( [ bhx0 bhx1 ], { 'Negative', sprintf('%s', title0)} )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g; p=%0.2g', myus( 1 ), myus( 2 ), pval ) );
xlabel( str )

return

% make_two_hist

function [ myus, sds, pval, bins_x, hx0, hx1 ] = make_two_hist( x, tidx1, tidx2, nbins, colors, str, byprob, title0 )

% acquire data
x0 = x (tidx1 == 0 & tidx2 == 0);
x1 = x( tidx1 == 1 );
x2 = x( tidx2 == 1 );

% compute mean, SDs, and pval
myus = [ nanmean( x0 ) nanmean( x1 ) nanmean( x2 ) ];
sds = [ nanstd( x0 ) nanstd( x1 ) nanstd( x2 ) ];
pval1        = utest( x0, x1);
pval2        = utest( x0, x2);
pval3        = utest( x1, x2);
pval         = [pval1; pval2; pval3];
% generate histogram
 bords_x           = minmax( [ x0; x1; x2 ] );
edges_x          = linspace( bords_x( 1 ), bords_x( 2 ), nbins + 1 );
hx0              = histcounts( x0, edges_x );
hx1              = histcounts( x1, edges_x );
hx2              = histcounts( x2, edges_x );
% hx0( end - 1 )   = hx0( end - 1 ) + hx0( end );
% hx0( end )       = [];
% hx1( end - 1 )   = hx1( end - 1 ) + hx1( end );
% hx1( end )       = [];
% hx2( end - 1 )   = hx2( end - 1 ) + hx2( end );
% hx2( end )       = [];
bins_x           = ( edges_x( 1 : end - 1 ) + edges_x( 2 : end ) ) / 2;
%  [N,EDGES,BIN] = HISTCOUNTS

if byprob
    hx0         = hx0 / sum( hx0 );
    hx1         = hx1 / sum( hx1 );
    hx2         = hx2 / sum( hx2 );
end
       
% plot histogram
newplot
bhx0             = stairs( bins_x, hx0, 'color', colors( 1, : ), 'linewidth', 1.5 );
hold on
bhx1             = stairs( bins_x, hx1, 'color', colors( 2, : ), 'linewidth', 1.5 );
hold on
bhx2             = stairs( bins_x, hx2, 'color', colors( 3, : ), 'linewidth', 1.5 );
% legend( [ bhx0 bhx1 bhx2 ], { sprintf('%s, %s, %s', title0{1},title0{2} , title0{3})} )
legend( [ bhx0 bhx1 bhx2 ],  title0{1},title0{2} , title0{3})

set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g/%0.2g; p01=%0.2g; p02=%0.2g', myus( 1 ), myus( 2 ),myus( 3 ), pval1, pval2 ) );
xlabel( str )

return

%----------------------------------------------------
% add_iso_ellipses

function h = add_iso_ellipses( param1, param2, labels )

% constants
colors              = [ 0 0 0.7; 1 0 0 ];
ngroups             = 2;
nSDs                = [ 2 4 8 16 32 ];

% compute
nm                  = length( nSDs );
tmp                 = cell( 1, ngroups );
for i               = 1 : ngroups
    ct              = i - 1;
    idx             = labels == ct;
    gmfit           = fitgmdist( [ param1( idx ) param2( idx ) ], 1 );
    tmp{ i }        = gmfit;
end
mixp                =  [   1 - mean( idx ) mean( idx ) ];
mu                  =  [ tmp{ 1 }.mu; tmp{ 2 }.mu ];
Sigma(:,:,1)        = tmp{ 1 }.Sigma;
Sigma(:,:,2)        = tmp{ 2 }.Sigma;
gm                  = gmdistribution( mu, Sigma, mixp );

gm1.mu              = gm.mu( 1, : ); 
gm1.Sigma           = gm.Sigma( :, :, 1 );
gm2.mu              = gm.mu( 2, : ); 
gm2.Sigma           = gm.Sigma( :, :, 2 );

% plot
hold on
xlims               = xlim;
ylims               = ylim;
h                   = NaN( nm, ngroups );
for i               = 1 : ngroups
    for j           = 1 : nm 
        m           = nSDs( j );
        if i == 1
            Cxy     = m * gm1.Sigma;
            Mxy     = gm1.mu;
        else
            Cxy     = m * gm2.Sigma;
            Mxy     = gm2.mu;
        end
        Rho         = Cxy( 1, 2 ) / sqrt( prod( diag( Cxy ) ) );
        if Rho > 0
            Ang     = Rho * pi / 4;
        else
            Ang     = pi + Rho * pi / 4;
        end
        h( j, i )   = ellipse( Mxy( 1 ), Mxy( 2 ), Cxy( 1, 1 ), Cxy( 2, 2 ), Ang, colors( i, : ) );
    end
end
set( gca, 'xlim', xlims, 'ylim', ylims )

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
sd0 = sst.sd{ uidx };

[ vA, vT, bpi ] = calc_spatial_waveform_features( w0, sd0, [], 1 );

% run this over all units:
nunits = size( sst.shankclu, 1 );
ivals = NaN * ones( nunits, 1 );
vAs = NaN * ones( nunits, 11 );
vTs = NaN * ones( nunits, 11 );
for i = 1 : nunits
    w = sst.mean{ i };
    sd = sst.sd{ uidx };
    [ vA, vT, bpi ] = calc_spatial_waveform_features( w, sd);
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
