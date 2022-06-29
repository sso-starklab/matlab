% make_punits_figures           make figures for positive units paper

% 17-may-20 SSo + ES

% revisions
% 31-may-20 (1) completed check_mono derivation of connectivity maps
%           (2) updated figure 3, part 8
% 01-june-20(1) updated datadir and outdir for ispc
% 02-jun-20 (2) added extraction of ASG-e and ASG-i from s2s
%           (3) added plotting of ASG
%           (4) added significance testing for ASG, nPRE, nPOST
% 04-jun-20 (1) added initial place analyses (Fig. 5)
%           (2) added ubiquity measures (Fig. 2)
% 07-jun-20 (1) added ACH_COM, baseline rate, and fanofactor analyses (Fig. 4)
%           (2) externalized make_punits_figures_one_hist
% 10-jun-20 (1) added depth analyses and plots
%           (2) updated datadir and outdir for isunix
%           (3) added pie chart plot (Fig. 2)
%           (4) added Lratio and ISI-index plots (Fig. 4) - change xlim?
% 18-jun-20 (1) added optical responses analysis (preliminary)
% 21-jun-20 (1) finished data accumulation for optical responses
% 24-jun-20 (1) extended data accumulation for optical responses to datadir
%           (2) updated datadir for example in figure 3 for isunix
% 25-jun-20 (1) fixed column names for bar graphs in figure 3
%           (2) added lab convention colors 
%           (3) added figure 8 - punit sub-clusters
% 28-jun-20 (1) added tsne for figure 8
%           (2) added multiple 2D projections to figure 8 (all units, punits)
%           (3) added ubiquity by spikes and by regions
%           (4) added mA234 and mP20 session regions
% 07-jul-20 (1) added regions to all sessions: nCx,CA1,hpx,WM,other (removed reg_CA1_DG, reg_mixed), and updated cases in fig2 (4)
% 12-jul-20 (1) organized figure 7 (optical responses) according to
%               multiple opsin types/targeted cells, and cell types (manipulated cells)
% 13-jul-20 (1) moved to sub-cluster on the local version (not on cluster solution loaded), and changed clunames, 
%               colors and tsne parameters to support over 3 clusters
% 16-jul-20 (1) added cuts through the FF lines and plotted distributions
%               of FF for each unit type, at each window size
%           (2) computed distributions for number of POST peers
%               to do - plot these histograms
% 19-jul-20 (1) started writing algorithm for automatic selection of number
%               of clusters (fig. 8)
% 23-jul-20 (1) finished writing algorithm for automatic selection of number
%               of clusters (fig. 8); added argument clustRunMode
%           (2) added punit cluster selection for all figures
%           (3) corrected punit/nunit tags for fig. 4.1
% 27-jul-20 (1) updated dir for depth of linear probes
%           (2) updated mDL5 regions (fig2), added reg_sub
% 30-jul-20 (1) modified fig5 to run on Hadas's function 'make_PUNIT_figures'
%           (2) added option for log scale axis for fig2.3 subplot( 2, 2, 3 )
% 05-aug-20 (1) updated regions (fig2) - moved from "others" 
% 09-sep-20 (1) updated regions (fig2) - new sessions (mDS1, mDS2)
%           (2) added precentage report fot fig2.4
% 13-sep-20 (1) added flipping and approximation of subplots in fig. 6
% 14-sep-20 (1) added subplots of depth vs width, t2p, firing rate and burstiness in fig. 6
%           (2) removed linear probes from regions in fig2 (because they can include more than one region)
% 11-oct-11 (1) fig. 2: computed SEM for fig2.3
%                       and computed means_all for fig2.4
%           (2) fig. 6: added path to depth directories (presently only for ismac)
%                       scaled depth estimates according to probe resolution (15 vs. 20 um)
%                       added g-test comparing punit and nunit depth distributions
%           (3) fig. 3: plotted fraction of excitatory/inhibitory units as two separate panels
%                       added binomial p-value comparing number E/I units to chance
%                       added g-test comparing fractions of E-punits to E-nunits (and same for I-)
% 29-dec-20 updated fname_hadas for isunix
% 07-jan-21 (1) fig. 4 - fixed bug of stepping over figure 3
%           (2) added cluster num in figs. 2, 3, 4, 5, 6
%           (3) modified call to fig. 5 (make_PUNIT_figures) to dilute by sst
% 13-jan-21 (1) removed region cell arrays from fig2
%           (2) added region field to sst structure after 'onlygather'
%           (3) fig. 2 - added stats on number of punit per region
%           (4) fig. 2 - devided pie charts to NInt and NPyr in case of clunum
% 14-jan-21 added output structure for fig. 2
% 19-jan-21 fig. 2 - added plot of waveforms for subclusters
% 20-jan-21 (1) added ridx_pos to remove 20 nunits from ispos structure that were mistakenly classified as punits
%           (2) fig. 2 - added mean wf and sem for each subcluster
%           (3) added colors for subclusters
%           (4) fig. 2,3 - added meta
% 08-feb-21 (1) added checkSync variable (Fig. 3)
%           (2) added mono.pairsSync analysis (Fig. 3)
% 10-feb-21 added probe type field to sst structure after 'onlygather' (addprobe)
% 26-apr-21 (1) modified fig. 2 to include bipolar units
%           (2) modified colors_NP to include bipolar color (colors_NPB)
% 27-apr-21 (1) added hack for realign bug
%           (2) modified fig. 3 to include bipolar units
% 13-may-21 (1) modified fig. 6 to include bipolar units
% 09-mar-22 (1) modified fig. 7 to include bipolar units
% 14-mar-22 (1) modified call to make_PUNITS_figures in fig. 7 
%           (2) started working on adding gcch and t to fig. 3

% shirly to do:
% fig. 8 - apply smart x/y lims
% fig. 4 - apply correction for sst2 in histograms

% to do:
% figure 8 - make scatter diagrams of all other pairs
% (5) consider wideband spike waveforms (without feature extraction)


% figure 1: 
% figure 2: ubiquity
% figure 3: connectivity - includes loading of data (check_mono)
% figure 4: firing rate statistics (fano factor, ACH_com, Lratio, ISI-index)
% figure 5: place analyses
% figure 6: depth analyses
% figure 7: optical responses
% figure 8: subclusters

% to sync:
% rsync -avh --progress /Volumes/odin/Shirly/eps/xml/*.prm.xml /Users/eranstark/Documents/da/punits/xml/
% % 185
% rsync -avh --progress /Volumes/odin/Shirly/eps/sst/*.sst /Users/eranstark/Documents/da/punits/sst/
% % 187
% rsync -avh --progress /Volumes/odin/Shirly/eps/s2s/*.s2s /Users/eranstark/Documents/da/punits/s2s/
% % 189
% rsync -avh --progress /Volumes/odin/Shirly/eps/sps/*.sps /Users/eranstark/Documents/da/punits/sps/
% % 165
% rsync -avh --progress /Volumes/odin/Shirly/eps/fanofactor/*fanofactors* /Users/eranstark/Documents/da/punits/fanofactor/
% % 150
% rsync -avh --progress /Volumes/odin/Shirly/eps/dc/*celltypeClassification /Users/eranstark/Documents/da/punits/dc/
% % 66

%function make_punits_figures( fignums, savef, savetype, outdir, varargin )
function rez = make_punits_figures( fignums, varargin )


% general constants
Nfigs                   = 10;
Nsubfigs                = 14;

colors_PI               = [ 46 125 50; 106 27 154 ] / 255; % INT, PYR
colors_NP              = [ 75 54 145; 173 20 87] / 255; % Nunit, Punit
colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Bipolar
colors_PI_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR
colors_NP_light         = [ 94 146 243; 227 81 131 ] / 255; % Nunit, Punit
colors_sub              = [ 0 163 126; 145 47 114; 42 101 176; 2 83 65 ] / 255; % subclusters
colors_sub_light        = [ 90 219 179; 241 136 203; 106 151 232; 0 169 127 ] / 255; % subclusters

% argument handling
nargs = nargin;
if nargs < 1 || isempty( fignums )
    fignums = 0;
end
if ~isempty( setdiff( fignums, [ 0 1 ] ) ) || length( fignums ) < Nfigs
    idx = round( fignums );
    fignums = zeros( Nfigs, 1 );
    fignums( idx ) = 1;
end

[ savef, savetype, pstr, outdir ...
    , clustRunMode, clunum, dolinear ...
    , metaplot, rez ... 
    , checkSync] = ParseArgPairs(...
    { 'savef', 'savetype', 'pstr', 'outdir' ...
    , 'clustRunMode', 'clunum', 'dolinear' ...
    , 'metaplot', 'rez' ...
    , 'checkSync'...
    }...
    , { 1, 'pdf', '-dpdf', '' ...
    , 'loadRun', NaN, 0 ...
    , 0, [] ...
    ,0 ...
    }...
    , varargin{ : } );
if isequal( pstr, '-dpdf' )
    resize                      = '-bestfit';
else
    resize                      = '';
end

if ~ismember( clustRunMode, { 'loadRun', 'saveRun', 'localRun' } )
    error( 'mismatch' )
end
if isnan( clunum )
    clunumstr               = '';
else
    clunumstr               = sprintf( '_subclu%d', clunum );
end
if isempty( rez ) || ~metaplot
    rez                     = cell( Nfigs, Nsubfigs );
end

% get some inline functions
[ ~, bino_ci_norm, bino_se_norm ] = binomial_inlines;
binp                                            = inline_stats;

% set plotting parameters:
sea_datadir                 = [];
if isempty( outdir ) || ~exist( outdir, 'dir' )
    if ismac        % eran (hoor)
        outdir          = '/Users/eranstark/Documents/graphics/punits/data_figures';
        datadir         = '/Users/eranstark/Documents/da/punits/';
        % assumes that data are synchronized with odin, if not, write:
        % rsync -avh --progress /Volumes/odin/Shirly/eps/fanofactor/* ~/Documents/da/punits/fanofactor/
        % rsync -avh --progress /Volumes/odin/Shirly/eps/sps/* ~/Documents/da/punits/sps/ 
        % rsync -avh --progress /Volumes/odin/Shirly/eps/sst/* ~/Documents/da/punits/sst/
        % rsync -avh --progress /Volumes/odin/Shirly/eps/s2s/*  ~/Documents/da/punits/s2s/ 
        %datadir_hadas       = '/Users/eranstark/Documents/da/hifi/';
        %fname_hadas         = 'all_animals_data_21May20';
        datadir_hadas       = '/Users/eranstark/Documents/da/punits/pf/';
        fname_hadas         = 'TPR_23-Aug-2021_Z6_extremum_bpi';
        if fignums( 6 )
            if dolinear
                sea_datadir        = '/Users/eranstark/Documents/da/punits/sst/depth_linear/'; % for linear probes
            else
                sea_datadir        = '/Users/eranstark/Documents/da/punits/sst/depth/'; % for dense probes
            end
        end
    elseif isunix   % shirly (nanna) Linux
        outdir              = '/shirly2/_Shirly/Lab/eps/data_figures';
        datadir             = '/probox1/mice/EPS/';
        datadir_hadas       = '/shirly2/_Shirly/Lab/eps/';
        fname_hadas         = 'TPR_23-Aug-2021_Z6_extremum_bpi';
        if fignums( 6 )
            if dolinear
                sea_datadir        = '/shirly2/mice/EPS/sst/depth_linear/'; % for linear probes
            else
                sea_datadir        = '/shirly2/mice/EPS/sst/depth/'; % for dense probes
            end
        end
    elseif ispc     % shirly (nanna) Windows
        outdir              = 'D:\_Shirly\Lab\eps\data_figures';
        datadir             = 'G:\mice\EPS\';
        datadir_hadas       = 'D:\_Shirly\Lab\eps\';
        fname_hadas         = 'TPR_23-Aug-2021_Z6_extremum_bpi';
    end
end
% loadname                = [ datadir '_ridx_pos.mat' ];
% L0 = load(loadname,'-mat');
% ridxpos = L0.idx_new;
% ridxpos = [];
% for depth:
%        datadir             = '/media/shirly/C22865A128659567/mice/EPS/sst/depth_linear';

% get the database
[ res, sst, ff ]            = shirly_eps_analysis( sea_datadir, 'onlygather', 1, 'ff_suffix', '_fanofactors_SWS' );
sst                         = addregion (sst); % add session region to sst structure
sst                         = addprobe (sst); % add probe type to sst structure
% sst.extremum(ridxpos)       = 0; % in order to remedy non-correct extremum detection (due to alignment issues) 
% use only data from a single subcluster
if ~isnan( clunum )
    fignums( 8 )            = 0;
    L                       = load( [ datadir 'punits_clusters.mat' ] );
    ispos                   = sst.extremum > 0;
    if ~isequal( size( L.x, 1 ), size( sst.filebase, 1 ) ) || ~isequal( sum( ispos ), length( L.clu ) )
        error( 'mismatch - apparently you added new data, re-run fignums 8 with saveRun' )
    end
    % add option of 'allclus'
    pnums                   = find( ispos );
    cidx                    = L.clu == clunum;
    if sum( cidx ) == 0
        fprintf( 1, 'No punit cluster named %d!\n', clunum )
        return
    else
        fprintf( 1, 'Focusing on %d punits from cluster %d and %d nunits\n', sum( cidx ), clunum, sum( ~ispos ) )
    end
    rmv                     = pnums( ~cidx );
    kidx                    = true( size( sst.filebase, 1 ), 1 );
    kidx( rmv )             = 0;
    sst                     = struct_select( sst, kidx );
end

% hack for realign bug (exremum on 17th sample instead of 16th)
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
    isbip           = ~isnan(sst.bpi) & ~isinf( sst.bpi );
    ispos           = sst.extremum > 0 & ~isbip;
    
% ----------------------------------------------------------------------
% Figure 1
% 
% ----------------------------------------------------------------------

if fignums( 1 )
    
    % process the data
    
    % plot the figures
    fig1( 1 ) = figure;
    
    fig1( 2 ) = figure;

    %-----
    % save the figures
    fig = fig1;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG1_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 2
% Ubiquity
% ----------------------------------------------------------------------

if fignums( 2 ) && ~metaplot
      
    % cell array of mouse lines
    lines_VIP         = {'mV99'};
    lines_CamK        = {'mF84', 'mF79', 'mF93', 'mF108'};
    lines_PV          = {'mP23', 'mP101', 'mDL5'};
    lines_FVB_CamK    = {'mC41', 'mA234', 'mS234'};
    lines_CCK         = {'mK01'};
    lines_FVB_c57B    = {'mF105'};
    lines_FVB_PV      = {'mB142'};


    % plot the figures
    
    % histograms of signed amplitudes:
    fig2( 1 ) = figure;
    
    % copied from shirly_eps_analysis:
%     ispos           = sst.extremum > 0;

    %nunits          = size( ispos, 1 );
    %nbins           = ceil( nunits / 5 );

    nbins           = 60;
    
    byprob          = 0;
    
    [ myus, sds, pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.maxp2p .* ( 2 * ispos - 1 ) * 1000, ispos, nbins, colors_NP, 'Signed Amp [\muV]', byprob );
    xlim ([-1500 1500])
    if ~isnan( clunum )
        textf( 0.5, 0.975, sprintf( 'cluster %d', clunum ) );
    end
    rez{ 2, 1 }.hx0                 = hx0;
    rez{ 2, 1 }.hx1                 = hx1;
    rez{ 2, 1 }.bins_x              = bins_x;
    rez{ 2, 1 }.myus                = myus;
    rez{ 2, 1 }.sds                 = sds; 
    rez{ 2, 1 }.pval                = pval;
    
    %----------------------------------------------------------------
    % pie charts
    fig2( 2 )                   = figure;
    slice_names                 = {'Nunits', 'Punits', 'Bipolar'};
    if isnan( clunum ) || isempty (clunum)
        for spi                     = 1 : 2
            subplot( 1, 2, spi )
            switch spi
                case 1
                    % pie chart for the fraction of the punits out of the total units
                    sums            = [ sum( ~ispos & ~isbip ) sum( ispos ) sum( isbip )];
                case 2
                    % pie chart for the fraction of the punits out of the total units
                    sums            = [ sum( sst.nspks( ~ispos & ~isbip ) ) sum( sst.nspks( ispos ) ) sum( sst.nspks( isbip ) )];
            end
            ph                      = pie( sums, slice_names );
            rez{ 2, 2 }.slice_names( spi, : )   = slice_names;
            rez{ 2, 2 }.sums( spi, : )          = sums;
            for i                   = 1 : 3
                set( ph( 2 * i - 1 ), 'FaceColor', colors_NPB( i, : ), 'EdgeColor', colors_NPB( i, : ) )
                set( ph( 2 * i ), 'String', sprintf( '%s (%d)', slice_names{ i }, sums( i ) ) )
            end
            title( sprintf( '%0.2g%%', round( sums( 2 ) / sum( sums ) * 1000 ) / 10 ) )
        end
        
    else
        slice_names                 = {'NInt','NPyr', 'Punits'};
        colors_clu = [colors_PI; colors_NPB(2, :)];
        for spi                     = 1 : 2
            subplot( 1, 2, spi )
            switch spi
                case 1
                    % pie chart for the fraction of the punits out of the total units
                    sums            = [ sum( ~ispos&~sst.pyr ) sum( ~ispos&sst.pyr ) sum( ispos ) ];
                case 2
                    % pie chart for the fraction of the punits out of the total units
                    sums            = [ sum( sst.nspks( ~ispos&~sst.pyr ) ) sum( sst.nspks( ~ispos&sst.pyr ) ) sum( sst.nspks( ispos ) ) ];
            end
            ph                      = pie( sums, slice_names );
            rez{ 2, 2 }.slice_names( spi, : )   = slice_names;
            rez{ 2, 2 }.sums( spi, : )          = sums;
            for i                   = 1 : 3
                set( ph( 2 * i - 1 ), 'FaceColor', colors_clu(i,:), 'EdgeColor', colors_clu(i,:) )
                set( ph( 2 * i ), 'String', sprintf( '%s (%d)', slice_names{ i }, sums( i ) ) )
            end
            title( sprintf( '%0.2g%%', round( sums( 3 ) / sum( sums ) * 1000 ) / 10 ) )
        end
    end
    
    %----------------------------------------------------------------
    % histograms/scatter of number of punits/nunits/bipolar per session
    usess           = unique( sst.filebase );
    nsess           = length( usess );
    nums            = NaN( nsess, 3 );
    for i           = 1 : nsess
        asess       = usess{ i };
        idx         = ismember( sst.filebase, asess );
        nums( i, : ) = [ sum( idx & ~ispos & ~isbip ) sum( idx & ispos & ~isbip) sum( idx & ~ispos & isbip) ];
    end
    frcts           = nums ./ [ sum( nums, 2 ) * ones( 1, 3 ) ];
    if ~isempty( clunum ) || ~isnan( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) );
    end
    
    fig2( 3 )       = figure;
    
    subplot( 2, 3, 1 )
    sessnums        = 1 : nsess;
    sessPunitNums   = sort( nums( :, 2 ) );
    %bh              = bar( 1 : nsess, sort( nums( :, 2 ) ), 1 );
    bh              = bar( sessnums, sessPunitNums, 1 );
    set( bh, 'FaceColor', colors_NPB( 2, : ), 'EdgeColor', colors_NPB( 2, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number of Punits' )
    xlabel( 'Session' )
    mean_nums       = mean(nums( :, 2 ) );
    sem_nums       = calc_sem( nums( :, 2 ) );
    n_nums         = size( nums, 1 );
    tstr            = sprintf( '%0.2g +- %0.2g, n=%d', mean_nums, sem_nums, n_nums );
    title( tstr )
    alines( mean_nums, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    axis square
    rez{ 2, 3 }.sessnums            = sessnums;
    rez{ 2, 3 }.sessPunitNums       = sessPunitNums;

    subplot( 2, 3, 4 )
    sessnums        = 1 : nsess;
    sessBipNums   = sort( nums( :, 3 ) );
    %bh              = bar( 1 : nsess, sort( nums( :, 2 ) ), 1 );
    bh              = bar( sessnums, sessBipNums, 1 );
    set( bh, 'FaceColor', colors_NPB( 3, : ), 'EdgeColor', colors_NPB( 3, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number of Punits' )
    xlabel( 'Session' )
    mean_nums       = mean(nums( :, 3 ) );
    sem_nums       = calc_sem( nums( :, 3 ) );
    n_nums         = size( nums, 1 );
    tstr            = sprintf( '%0.2g +- %0.2g, n=%d', mean_nums, sem_nums, n_nums );
    title( tstr )
    alines( mean_nums, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    axis square
    rez{ 2, 3 }.sessBipNums       = sessBipNums;

    subplot( 2, 3, 2 )
    %bh              = bar( 1 : nsess, sort( frcts( :, 2 ) ), 1 );
    sessPunitFracts                 = sort( frcts( :, 2 ) );
    bh                              = bar( sessnums, sessPunitFracts, 1 );
    set( bh, 'FaceColor', colors_NPB( 2, : ), 'EdgeColor', colors_NPB( 2, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of Punits' )
    xlabel( 'Session' )
    mean_frcts                      = mean( frcts( :, 2 ) );
    sem_frcts                       = calc_sem( frcts( :, 2 ) );
    n_frcts                         = size( frcts, 1 );
    tstr                            = sprintf( '%0.2g +- %0.2g, n=%d', mean_frcts, sem_frcts, n_frcts );
    title( tstr )
    alines( mean_frcts, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    %alines( mean_frcts + [ -1 1 ] * sem_frcts, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    axis square
    rez{ 2, 3 }.sessPunitFracts     = sessPunitFracts;
    
    subplot( 2, 3, 3 )
    jit             = ( rand( nsess, 3 ) - 0.5 ) / 2;
    nums_jit        = nums + jit;
    nums_jit( nums == 0 ) = 0;
    scatter ( nums_jit( :, 1 ), nums_jit( :, 2 ), 15,'MarkerFaceColor','b')
%    scatter ( log2(nums_jit( :, 1 )), log2(nums_jit( :, 2 )), 15,'MarkerFaceColor','b')
    alpha(.5);
    [ cc, pp ] = calc_spearman( nums( :, 1 ), nums( :, 2 ), 1000 );
    title( sprintf( 'Number of units/session (CC=%0.2g, pval=%0.3g)', cc, pp ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Nunits' )
    ylabel( 'Punits' )
    rez{ 2, 3 }.nums_jit            = nums_jit;
 
    subplot( 2, 3, 6 )
    jit             = ( rand( nsess, 3 ) - 0.5 ) / 2;
    nums_jit        = nums + jit;
    nums_jit( nums == 0 ) = 0;
    scatter ( nums_jit( :, 1 ), nums_jit( :, 3 ), 15,'MarkerFaceColor','b')
%    scatter ( log2(nums_jit( :, 1 )), log2(nums_jit( :, 3 )), 15,'MarkerFaceColor','b')
    alpha(.5);
    [ cc, pp ] = calc_spearman( nums( :, 1 ), nums( :, 3 ), 1000 );
    title( sprintf( 'Number of units/session (CC=%0.2g, pval=%0.3g)', cc, pp ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Nunits' )
    ylabel( 'Bipolars' )
    
    subplot( 2, 3, 5 )
    slice_names = { 'Only Nunits', 'Only Punits', 'Punits+Nunits', 'Only Bipolars', 'Bipolars+Nunits', 'Punits+Bipolars+Nunits' };
    sums = [ sum( nums( :, 1 ) > 0 & nums( :, 2 ) == 0 & nums( :, 3 ) == 0 )
        sum( nums( :, 1 ) == 0 & nums( :, 2 ) > 0 & nums( :, 3 ) == 0 )
        sum( nums( :, 1 ) > 0 & nums( :, 2 ) > 0 & nums( :, 3 ) == 0) 
        sum( nums( :, 1 ) == 0 & nums( :, 2 ) == 0 & nums( :, 3 ) > 0)
        sum( nums( :, 1 ) > 0 & nums( :, 2 ) == 0 & nums( :, 3 ) > 0)
        sum( nums( :, 1 ) > 0 & nums( :, 2 ) > 0 & nums( :, 3 ) > 0)];
    for j = 1 : length (slice_names)
        slice_names{j} = sprintf('%s(%.2f%%)',slice_names{j}, (sums(j)/sum(sums)*100) );
    end
    ph = pie( sums, slice_names );
%        title( sprintf( '%0.2g%%', sums)
    try
        set( ph( 1 ), 'EdgeColor', colors_NPB( 1, : ), 'FaceColor', colors_NPB( 1, : ) )
        set( ph( 3 ), 'EdgeColor', colors_NPB( 2, : ), 'FaceColor', colors_NPB( 2, : ) )
        set( ph( 5 ), 'EdgeColor', colors_NPB( 3, : ), 'FaceColor', colors_NPB( 3, : ) )
    catch
        fprintf( 'Failed in something\n' )
    end
    if ~isempty( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) );
    end
    rez{ 2, 3 }.sums            = sums;
    rez{ 2, 3 }.sums            = slice_names;
    
    %----------------------------------------------------------------
    % summary for regions and lines
    fig2( 4 )           = figure;
    %nregs               = 6;
    nregs               = length(unique(sst.reg));

    means               = NaN( nregs, 2 );
    sems                = NaN( nregs, 2 );
    counts              = NaN( nregs, 1 );
    means_all               = NaN( nregs, 2 );
    sems_all                = NaN( nregs, 2 );
    counts_all              = NaN( nregs, 1 );
    
    sums1               = NaN( nregs, 1 );
    means1               = NaN( nregs, 1 );
    sems1               = NaN( nregs, 1 );
    
    reg_names           = { 'nCX', 'CA1', 'WM','OTHER' };

    for spi             = 1 : nregs
        
        switch spi
            case 1
                list1   = 0; % nCX
            case 2
                list1   = 1; % CA1
            case 3
                list1   = 2; % WM
            case 4
                list1   = 3; % OTHER          
        end
        name1           = reg_names{ spi };
        
        idx1            = ismember( sst.reg, list1 );
        sst1            = struct_select( sst, idx1 );
        ispos1          = ispos( idx1 );
        
        usess           = unique( sst1.filebase );
        nsess           = length( usess );
        nums1           = NaN( nsess, 2 );
        for i           = 1 : nsess
            asess       = usess{ i };
            idx         = ismember( sst1.filebase, asess );
            nums1( i, : ) = [ sum(idx& ~ispos1 )  sum(idx& ispos1 ) ];
        end
        frcts1          = nums1 ./ ( sum( nums1, 2 ) * ones( 1, 2 ) );
        forpv(spi).nums = nums1;
        sums1( spi )     = sum(nums1(:,2));
        means1( spi )    = mean (nums1(:,2));
        sems1( spi )     = calc_sem (nums1(:,2));

        means( spi, : ) = mean( frcts1, 1 );
        sems( spi, : )  = calc_sem( frcts1, 1 );
        counts( spi, : ) = size( frcts1, 1 );
        
        means_all( spi, : ) = [ sum( ~ispos1 ) sum( ispos1 ) ] / length( ispos1 );
        sems_all( spi, : ) = bino_se_norm( [ sum( ~ispos1 ) sum( ispos1 ) ], length( ispos1 ) );
        counts_all( spi, : ) = length( ispos1 );
        

        
        % plot
        subplot( 2, 4, spi )
        %bh              = bar( 1 : nsess, sort( frcts1( :, 2 ) ), 1 );
        sessnums                = 1 : nsess;
        sessPunitFracts         = sort( frcts1( :, 2 ) );
        bh                      = bar( sessnums, sessPunitFracts, 1 );
        set( bh, 'FaceColor', colors_NP( 2, : ), 'EdgeColor', colors_NP( 2, : ) )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        ylabel( 'Fraction of Punits' )
        xlabel( 'Session' )
        mean_frcts      = mean( frcts1( :, 2 ) );
        alines( mean_frcts, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
        title( sprintf('%s, total %d', name1, sums1(spi)) )
        axis square
        
        rez{ 2, 4 }.sessnums{ spi }             = sessnums;
        rez{ 2, 4 }.sessPunitFracts{ spi }      = sessPunitFracts;
        
    end

    subplot( 2, 4, spi + 1 )
    %bar( 1 : nregs, means( :, 2 ) )
%     [ bh, eh ] = barwerror( 1 : nregs, means( :, 2 ), sems( :, 2 ) ...
%         , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
    reg_nums            = 1 : nregs;
    reg_means           = means_all( :, 2 );
    reg_sems            = sems_all( :, 2 );
    [ bh, eh ] = barwerror( reg_nums, reg_means, reg_sems ...
        , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
%     [ bh, eh ] = barwerror( 1 : nregs, means_all( :, 2 ), sems_all( :, 2 ) ...
%         , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
    set( gca, 'XTickLabel', reg_names );
    ylabel( 'Fraction of Punits' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    axis square
    rez{ 2, 4 }.reg_nums            = reg_nums;
    rez{ 2, 4 }.reg_means           = reg_means;
    rez{ 2, 4 }.reg_sems            = reg_sems;
    
    subplot( 2, 4, spi + 2 )
    %bar( 1 : nregs, means( :, 2 ) )
%     [ bh, eh ] = barwerror( 1 : nregs, means( :, 2 ), sems( :, 2 ) ...
%         , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
%     [ bh, eh ] = barwerror( 1 : nregs, means1, sems1 ...
%         , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
    [ bh, eh ] = barwerror( reg_nums, means1, sems1 ...
        , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
    set( gca, 'XTickLabel', reg_names );
    ylabel( 'Number of Punits' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    axis square
    rez{ 2, 4 }.means1           = means1;
    rez{ 2, 4 }.sems1            = sems1;
    
    pv                  = NaN( 6, 1 );
    for n = 1 : 6
        mat5 = [1 2; 1 3; 1 4; 2 3 ;2 4;3 4];
        pv (n) = utest(forpv(mat5(n,1)).nums(:,2),forpv(mat5(n,2)).nums(:,2));
        fprintf('pv_%d_%d = %0.2g \n', mat5(n,:),pv(n))
    end
    
    if ~isnan( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) );
    end
    
    if ~isnan( clunum )
        sstp = struct_select (sst, ispos);
        nunits = length(sstp.filebase);
        meanwf_clu = NaN(nunits,32);
        fig2(5) = figure;
        subplot (1, 2, 1)
        for p = 1 : nunits
            wf = sstp.mean{p,1}(:,16);
            elec = find(ismember(wf,max(wf)));    % because we are only doing it for Punits
            meanwf_clu(p,:) = sstp.mean{p,1}(elec,:);
            plot(meanwf_clu(p,:))
            hold on
        end
        axis square
        set( gca, 'tickdir', 'out', 'box', 'off' )
        subplot (1, 2, 2)
        meanwf = mean(meanwf_clu);
        semwf = calc_sem(meanwf_clu,1);
        patch_band( 1:32,meanwf, semwf, colors_sub(clunum,:), colors_sub_light(clunum,:) );
        axis square
        set( gca, 'tickdir', 'out', 'box', 'off' )
        
        rez{ 2, 5 }.wf_means           = meanwf_clu;
        rez{ 2, 5 }.wf_mean            = meanwf;
        rez{ 2, 5 }.wf_sem             = semwf;
%         if clunum ==4 %clean data
%             figure;
%             l=0;
%             for p = 1 : nunits
%                 wf = sstp.mean{p,1}(:,7);
%                 elec = find(ismember(wf,min(wf)));    % because we are only doing it for Punits
%                 if min(wf)<-200
%                     l=l+1;
%                     meanwf_clu(p,:) = sstp.mean{p,1}(elec,:);
%                     plot(meanwf_clu(p,:))
%                     fileremove{l} =  sstp.filebase{p};
%                     shankremove(l,:) = sstp.shankclu(p,:);
%                     idxp(l) = p;
%                 hold on
%                 end
%             end   
%             figure;
% %             meanwf_clu(idxp,:)=[];
% %             nunits1 = length(meanwf_clu);
%             for p = 1 : nunits
%                 if ismember(p,idxp)
%                     continue
%                 end
%                 wf = sstp.mean{p,1}(:,16);
%                 elec = find(ismember(wf,max(wf)));
%                 meanwf_clu(p,:) = sstp.mean{p,1}(elec,:);
%                 plot(meanwf_clu(p,:))
%                 hold on
%             end
%         end
        fig_title( sprintf( 'cluster %d \n total: %d Punits', clunum, sum(ispos) ) );
    end
    
% idx_new = false(size(sst.shankclu,1),1);
% 
% for i = 1:length(shankremove)
%     idxf = find(ismember(sst.filebase,fileremove{i}));
%     for j = 1:length(idxf)
%         if ismember(sst.shankclu(idxf(j),:),shankremove(i,:),'rows')
%             idx_new(idxf(j)) = true;
%             idxr(i) = idxf(j);
%             shankremove2(i,:) = sst.shankclu(idxf(j),:);
%             fileremove2{i} =  sst.filebase{idxf(j)};
% 
%             break
%         end
%     end
% end
% 
% savename                = [ datadir '_ridx_pos.mat' ];
% save( savename, 'idx_new' )



    %-----
    % save the figures
    fig = fig2;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG2_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
end    

if fignums( 2 ) && metaplot
    
        
    % meta-plot of clusters of fig.2:
    clu_names = { 'Clu1', 'Clu2', 'Clu3' };
    reg_names           = { 'nCX', 'CA1', 'WM','OTHER' };
    nregs = length( reg_names );
    nclus = length( clu_names );
    
    fig2m( 1 )           = figure;
    
    byprob          = 1;
    if byprob
       rez{1}{ 2, 1 }.hx0                     = rez{1}{ 2, 1 }.hx0 / sum( rez{1}{ 2, 1 }.hx0 );
    end
    bhx0                        = stairs( rez{1}{ 2, 1 }.bins_x, rez{1}{ 2, 1 }.hx0, 'color', colors_NP( 1, : ), 'linewidth', 1.5 );
    for clunum = 1 : nclus % loop over clusters
        hold on
        if byprob
            rez{clunum}{ 2, 1 }.hx1                 = rez{clunum}{ 2, 1 }.hx1 / sum( rez{clunum}{ 2, 1 }.hx1 );
        end
        bhx1(clunum)            = stairs( rez{clunum}{ 2, 1 }.bins_x, rez{clunum}{ 2, 1 }.hx1, 'color', colors_sub( clunum, : ), 'linewidth', 1.5  );
        set( gca, 'tickdir', 'out', 'box', 'off' );
    end
    lh                      = legend( [ bhx0 bhx1(1:end) ], { 'Negative', 'P1', 'P2', 'P3', 'P4' } );
    xlim ([-1500 1500])

%     fig2m( 2 )           = figure;
%     for clunum = 1 : nclus % loop over clusters
%         for spi = 1:2
%             subplot( 1, 2, spi )
%             slice_names = rez{clunum}{ 2, 2 }.slice_names( spi, : );
%             sums = rez{clunum}{ 2, 2 }.sums( spi, : );            
%             ph                      = pie( sums, slice_names );
%             for i                   = 1 : 2
%                 set( ph( 2 * i - 1 ), 'FaceColor', colors_PI( i, : ), 'EdgeColor', colors_PI( i, : ) )
%                 set( ph( 2 * i ), 'String', sprintf( '%s (%d)', slice_names{ i }, sums( i ) ) )
%             end
%             set( ph( 2 * clunum +4 ), 'String', sprintf( '%0.2g%% P%d (%d)', round( sums( 2+clunum ) / sum( sums ) * 1000 ) / 10, clunum, sums( 3 ) ) )
%             hold on
%         end
%          set( ph( 2 * clunum + 4 ), 'FaceColor', colors_sub( clunum, : ), 'EdgeColor', colors_sub( i, : ) )
% %        title( sprintf( '%0.2g%%', round( sums( 2+clunum ) / sum( sums ) * 1000 ) / 10 ) )
% %         set( ph( 2 * clunum +2 ), 'String', sprintf( 'P%d (%d)', clunum{ i }, sums( 3 ) ) )
%         hold on
%     end
    
    fig2m( 2 )           = figure;
    for spi = 1 : nregs
        means = NaN( nclus, 1 );
        sems = NaN( nclus, 1 );
        for clunum = 1 : nclus % loop over clusters
            ridx = rez{ clunum }{ 2, 4 }.reg_nums == spi;
            means( clunum ) = rez{ clunum }{ 2, 4 }.reg_means( ridx );
            sems( clunum ) = rez{ clunum }{ 2, 4 }.reg_sems( ridx );
        end
        
        subplot( 2, 2, spi )
        barwerror( 1 : nclus, means, sems ...
            , colors_NP( 2, : ), 0.8, [ 0 0 0 ], 2 );
        set( gca, 'XTickLabel', clu_names );
        ylabel( 'Fraction of Punits' )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        axis square
        title( reg_names{ spi } )
    end

    fig2m( 3 )           = figure;
    for clunum = 1 : nclus       
        m1                      = min( rez{ clunum }{2,5}.wf_mean );            % define as the global minimum 
        m2                      = max( rez{ clunum }{2,5}.wf_mean );            % define as the global max
        m1                      = min( m1, 0 );              % force m1 to be non-positive
        m2                      = max( m2, 0 );              % force m2 to be non-negative
        m1                      = abs( m1 );                 % take the absolute value of the negative extremum
        m2                      = abs( m2 );                 %                   " positive "
        bpi (clunum)            = ( m1 - m2 ) ./ ( m1 + m2 );
        pb = patch_band( 1:32,rez{ clunum }{2,5}.wf_mean, rez{ clunum }{2,5}.wf_sem, colors_sub(clunum,:), colors_sub_light(clunum,:) );
        hold on
    end
        axis square
        set( gca, 'tickdir', 'out', 'box', 'off' )
        title( sprintf('Mean waveforms\n bpi: %0.4f, %0.4f, %0.4f', bpi))
               
    fig = fig2m;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG2m_part%d%s', outdir, i, clunumstr );
%             fig_out( fig( i ), 1, figname, savetype ),
            figure( fig( i ) );
            figi = gcf;
            figi.Renderer = 'painters';
            pause( 0.2 )
            print( fig( i ), pstr, figname, resize )
            
        end
    end
    
end

% ----------------------------------------------------------------------
% Figure 3
% Connectivity
% ----------------------------------------------------------------------

if fignums( 3 ) && ~metaplot
    
    % process the data
%     filebase        = filebaseLookup( 'mDL5', -16 );
%     mono            = check_mono( filebase );
%     
%     ilevel          = 'B';
%     L0              = load( [ filebase '.sst' ], '-mat' );
%     L1              = load( [ filebase '.s2s' ], '-mat' );
%     gidx            = check_cluster_quality( L0.sst, ilevel );    
%     mono1           = check_mono( L1.s2s, gidx );
    
    % go over global sst stucture and run check_mono for each session
    usess           = unique( sst.filebase );
    nsess           = length( usess );
    asg1all         = [];
    asg2all         = [];
    ind             = 0;
    for i = 1 : nsess
        % the indices in the global sst
        asess       = usess{ i };
        uidx        = ismember( sst.filebase, asess );
        % initialize an output structure for the session
        ssti            = struct_select( sst, uidx );
        nunits          = size( ssti.filebase, 1 );
        si.filebase     = repmat( ssti.filebase( 1 ), [ nunits 1 ] );
        si.shankclu     = NaN( nunits, 2 );
        lv              = false( nunits, 1 );
        nv              = zeros( nunits, 1 );
        si.exc          = lv;                           % boolean flag - is an excitatory unit
        si.inh          = lv;                           % boolean flag - is an inhibitory unit
        si.nexcPost     = nv;                           % number of excited post-synaptic peers
        si.ninhPost     = nv;                           % number of inhibited post-synaptic peers
        si.nexcPre      = nv;                           % number of exciting pre-synaptic peers 
        si.ninhPre      = nv;                           % number of inhibiting pre-synaptic peers 
        si.asg1         = nv;                           % mean ASG to excited post-synaptic peers
        si.asg2         = nv;                           % mean ASG to inhibited post-synaptic peers
        si.pairsSync    = nv;                           % number of units in sync 

        % load the sst and s2s for the relevant session
        monoOK          = 1;
        try
            L0          = load( [ datadir 'sst/' ssti.filebase{ 1 } '.sst' ], '-mat' );
            L1          = load( [ datadir 's2s/' ssti.filebase{ 1 } '.s2s' ], '-mat' );
        catch
            fprintf( 1, '%d: Missing data (either sst, s2s, or both) for %s\n', i, usess{ i } )
            monoOK      = 0;
        end
        if ~isequal( L1.s2s.shankclu, L0.sst.shankclu )
            fprintf( 1, '%d: Data mismatch for %s\n', i, usess{ i } )
            monoOK      = 0;
        end
        if monoOK
            % call check_mono locally with the indices of only the relevant units
            gidx        = ismember( L0.sst.shankclu, ssti.shankclu, 'rows' );
            mono        = check_mono( L1.s2s, gidx );
            % populate the structure with the relevant details
            si.exc( mono.exc ) = 1;
            si.inh( mono.inh ) = 1;
            si.pairsSync       = size(mono.pairsSync,1);

            % loop over all units, and counts number of peers for each
            for j = 1 : nunits
                si.nexcPost( j ) = sum( mono.pairsExc( :, 1 ) == j );
                si.ninhPost( j ) = sum( mono.pairsInh( :, 1 ) == j );
                si.nexcPre( j ) = sum( mono.pairsExc( :, 2 ) == j );
                si.ninhPre( j ) = sum( mono.pairsInh( :, 2 ) == j );
                asg1            = mono.pairsExc( mono.pairsExc( :, 1 ) == j, 5 );
                asg2            = mono.pairsInh( mono.pairsInh( :, 1 ) == j, 5 );
                si.asg1( j )    = nangeomean( asg1 );
                si.asg2( j )    = -nangeomean( -asg2 );
                asg1all         = [ asg1all; asg1 ];
                asg2all         = [ asg2all; asg2 ];
            end
            
            if checkSync
                if isunix
                    filebase_mono           = filebase_lookup(ssti.filebase{1});
                end 
                         s_mono1                  = load_spikes( filebase_mono );

                % go over all sync pairs
                for pairs               = 1:size(mono.pairsSync,1)
                    idx_shankclu1       = mono.pairsSync (pairs,1);
                    idx_shankclu2       = mono.pairsSync (pairs,2);
                    pidxi               = ssti.extremum > 0;

                    if pidxi(idx_shankclu1) || pidxi(idx_shankclu2) % if at least one of the synced units are Punits
%                         s_mono                  = load_spikes( filebase_mono );
                        s_mono = s_mono1;
                        par_mono                = LoadXml( filebase_mono );
                        ind                 = ind+1;

                        % remove spikes during stimuli
                        vals                    = LoadStims( filebase_mono );
                        if ~isempty (vals)
                            uvals                   = uniteranges( vals( :, 1 : 2 ) );
                            ridx                    = inranges( s_mono.res, uvals );
                            s_mono.clu( ridx )      = [];
                            s_mono.res( ridx )      = [];
                        end
                        
                        shankclu1               = mono.shankclu(idx_shankclu1,:);
                        shankclu2               = mono.shankclu(idx_shankclu2,:);
                        % get the spike times of the trigger
                        clunum1                 = s_mono.map( ismember( s_mono.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
                        clunum2                 = s_mono.map( ismember( s_mono.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
                        st1                     = s_mono.res( s_mono.clu == clunum1 );
                        st2                     = s_mono.res( s_mono.clu == clunum2 );

                        [ g1, g2, act, sil]     = calc_asg( st1, st2, 'cmode', 3, 'graphics', 0 );
                        [ g1, g2, act, sil,stc]     = calc_asg( st1, st2, 'cmode', 3, 'graphics', 0, 'BinSizeMS', 0.5 );
                        rez{ 3, 1 }.syncfilebase{ind,1}   = filebase_mono;
                        rez{ 3, 1 }.syncshankclu1(ind,1:2)  = shankclu1;
                        if pidxi(idx_shankclu1)
                            rez{ 3, 1 }.syncshankclu1(ind,3) = 2; % identity of the punit
                        else
                            rez{ 3, 1 }.syncshankclu1(ind,3) = 0;
                        end
                        rez{ 3, 1 }.syncshankclu2(ind,1:2)  = shankclu2;
                        if pidxi(idx_shankclu2)
                            rez{ 3, 1 }.syncshankclu2(ind,3) = 2;
                        else
                            rez{ 3, 1 }.syncshankclu2(ind,3) = 0;
                        end
                        rez{ 3, 1 }.syncASG(ind)        = g1;
                        
%                         %temp
%                         maxlag = (size(stc.cch,1)-1)/2;
%                         tidx = stc.g1base(:)-maxlag;
%                         kval                     = stc.gcch( stc.g1base( : ) );
%                         kval(isnan(kval))=0;
% 
%                         z                         = zeros( tidx( 1 ) - 1, 1 );
%                         STC                       = [ z; kval ];                              % [spks/s]
%                         kP                          = STC * stc.dt;                               % [spks/bin]
%                         ktP                         = ( tidx - 1 ) * stc.dt;      % [ms]
% 
%                         nk                       = length( kP );
%                         z0                          = zeros( nk - 1, 1 );
%                         k                       =  [ z0; kP ];                              % zeropad the kernel to be symmetric (but causal)
%                         t0                       = -( nk - 1 : -1  : 1 )' * stc.dt;
%                         kt                       = [ t0; ktP ];
% 
%                         % keep
%                         STCs{ j }                   = k;
%                         STC_time{ j }               = kt;
                        
                        
                    end
                end
            end
        end
        % accumulate over sessions
        if i == 1
            s       = si;
        else
            s       = struct_cat( s, si );
        end
            
    end
    
    
    % check if equal
    if ~isequal( s.filebase, sst.filebase  )
        error( 'fix all issues first' )
    end
    for i = 1 : size( s.filebase, 1 )
        eq( i ) = isequal( s.filebase( i ), sst.filebase( i ) ); 
    end
    neq = find( ~eq );
    % fix mK01_10 - mK1_10; mS###; mP23_06; mDL5_28; mP23_04_3
    
    
    % 31-may-20: note problem in check sum:
    % sum( [ double( s.exc ) - double( s.nexcPost > 0 ) ] )
    
    % plot the figures
%     pidx        = sst.extremum > 0;
    pidx        = ispos;
    bidx        = isbip;
    
    eidx        = s.exc == 1;                   % excitatory
    iidx        = s.inh == 1;                   % inhibitory
    eidx2       = s.nexcPre ~= 0;               % is post-synaptic to an excitatory (excited)
    iidx2       = s.ninhPre ~= 0;               % inhibited
    cidx        = eidx | iidx | eidx2 | iidx2;  % connected to the simultaneously recorded units
    
    % numbers/fractions of E-cells (punits/nunits/bipolars)
    emat        = NaN( 4, 3 );
    emat( 1 : 2, 1 ) = [ sum( pidx & eidx ), sum( pidx ) ]'; % Punits
    emat( 1 : 2, 3 ) = [ sum( bidx & eidx ), sum( bidx ) ]'; % Bipolars
    emat( 1 : 2, 2 ) = [ sum( ~pidx & ~bidx & eidx ), sum( ~pidx & ~bidx ) ]'; % Nunits
    emat( 3, : ) = emat( 1, : ) ./ emat( 2, : );
    emat( 4, : ) = bino_se_norm( emat( 1, : ), emat( 2, : ) )';
      
    % numbers/fractions of I-cells (punits/nunits/bipolars)
    imat        = NaN( 4, 3 );
    imat( 1 : 2, 1 ) = [ sum( pidx & iidx ), sum( pidx ) ]'; % Punits
    imat( 1 : 2, 3 ) = [ sum( bidx & iidx ), sum( bidx ) ]'; % Bipolars
    imat( 1 : 2, 2 ) = [ sum( ~pidx & ~bidx & iidx ), sum( ~pidx & ~bidx ) ]';
    imat( 3, : ) = imat( 1, : ) ./ imat( 2, : );
    imat( 4, : ) = bino_se_norm( imat( 1, : ), imat( 2, : ) )';
        
    % number/fractions of null-cells (punits/nunits/bipolars)
    nmat        = NaN( 4, 3 );
    nmat( 1 : 2, 1 ) = [ sum( pidx & ~eidx & ~iidx ), sum( pidx ) ]'; % Punits
    nmat( 1 : 2, 3 ) = [ sum( bidx & ~eidx & ~iidx ), sum( bidx ) ]'; % Bipolars
    nmat( 1 : 2, 2 ) = [ sum( ~pidx & ~bidx & ~eidx & ~iidx ), sum( ~pidx & ~bidx ) ]';
    nmat( 3, : ) = nmat( 1, : ) ./ nmat( 2, : );
    nmat( 4, : ) = bino_se_norm( nmat( 1, : ), nmat( 2, : ) )';

    % number/fractions of connected-cells (punits/nunits/bipolars)
    cmat        = NaN( 4, 3 );
    cmat( 1 : 2, 1 ) = [ sum( pidx & cidx ), sum( pidx ) ]';
    cmat( 1 : 2, 3 ) = [ sum( bidx & cidx ), sum( bidx ) ]'; % Bipolars
    cmat( 1 : 2, 2 ) = [ sum( ~pidx & ~bidx & cidx ), sum( ~pidx & ~bidx ) ]';
    cmat( 3, : ) = cmat( 1, : ) ./ cmat( 2, : );
    cmat( 4, : ) = bino_se_norm( cmat( 1, : ), cmat( 2, : ) )';
    
    % numbers/fractions of excited units (punits/nunits/bipolars)
    emat2       = NaN( 4, 3 );
    emat2( 1 : 2, 1 ) = [ sum( pidx & eidx2 ), sum( pidx ) ]';
    emat2( 1 : 2, 3 ) = [ sum( bidx & eidx2 ), sum( bidx ) ]'; % Bipolars
    emat2( 1 : 2, 2 ) = [ sum( ~pidx & ~bidx & eidx2 ), sum( ~pidx & ~bidx ) ]';
    emat2( 3, : ) = emat2( 1, : ) ./ emat2( 2, : );
    emat2( 4, : ) = bino_se_norm( emat2( 1, : ), emat2( 2, : ) )';
    
    % numbers/fractions of inhibited units (punits/nunits/bipolars)
    imat2       = NaN( 4, 3 );
    imat2( 1 : 2, 1 ) = [ sum( pidx & iidx2 ), sum( pidx ) ]';
    imat2( 1 : 2, 3 ) = [ sum( bidx & iidx2 ), sum( bidx ) ]'; % Bipolars
    imat2( 1 : 2, 2 ) = [ sum( ~pidx & ~bidx & iidx2 ), sum( ~pidx & ~bidx ) ]';
    imat2( 3, : ) = imat2( 1, : ) ./ imat2( 2, : );
    imat2( 4, : ) = bino_se_norm( imat2( 1, : ), imat2( 2, : ) )';
    
    rez{ 3, 1 }.emat           = emat;
    rez{ 3, 1 }.imat           = imat;
    rez{ 3, 1 }.emat2          = emat2;
    rez{ 3, 1 }.imat2          = imat2;
    rez{ 3, 1 }.cmat           = cmat;
    rez{ 3, 1 }.nmat           = nmat;

    % number of connected units (excitatory/inhibitory/excited/inhibited)
    
    % number of post-synaptic peers
    npost           = NaN( 3, 6 );
    npost( :, 1 )   = [ mean( s.nexcPost( eidx & ~pidx & ~bidx ) ) calc_sem( s.nexcPost( eidx & ~pidx & ~bidx ) ) sum( eidx & ~pidx & ~bidx ) ]';
    npost( :, 2 )   = [ mean( s.nexcPost( eidx & pidx ) ) calc_sem( s.nexcPost( eidx & pidx ) ) sum( eidx & pidx ) ]';
    npost( :, 3 )   = [ mean( s.ninhPost( iidx & ~pidx & ~bidx ) ) calc_sem( s.ninhPost( iidx & ~pidx & ~bidx ) ) sum( iidx & ~pidx & ~bidx ) ]';
    npost( :, 4 )   = [ mean( s.ninhPost( iidx & pidx ) ) calc_sem( s.ninhPost( iidx & pidx ) ) sum( iidx & pidx ) ]';
    npost( :, 5 )   = [ mean( s.nexcPost( eidx & bidx ) ) calc_sem( s.nexcPost( eidx & bidx ) ) sum( eidx & bidx ) ]';
    npost( :, 6 )   = [ mean( s.nexcPost( eidx & bidx ) ) calc_sem( s.nexcPost( eidx & bidx ) ) sum( eidx & bidx ) ]';

    % E-nunits, E-punits, I-nunits, I-punits
    pval_NPOST      = [ utest( s.nexcPost( eidx & ~pidx & ~bidx ), s.nexcPost( eidx & pidx ) )
        utest( s.ninhPost( iidx & ~pidx & ~bidx ), s.ninhPost( iidx & pidx ) ) ];
    [ ~, pp1 ]      = kstest2( s.nexcPost( eidx & ~pidx & ~bidx ), s.nexcPost( eidx & pidx ) );
    [ ~, pp2 ]      = kstest2( s.ninhPost( iidx & ~pidx & ~bidx ), s.ninhPost( iidx & pidx ) );
    pval_NPOST_KS2  = [ pp1 pp2 ];
    
    rez{ 3, 1 }.npost           = npost;

    ebinC                       = ( 1 : max( s.nexcPost( eidx ) ) )';
    nbins                       = length( ebinC );
    ecounts                     = zeros( nbins, 2 );
    for i                       = 0 : 1
        switch i
            case 0 
                [ cnts, bidxx ]  = uhist( s.nexcPost( eidx & ~pidx & ~bidx ) );
            case 1
                [ cnts, bidxx ]  = uhist( s.nexcPost( eidx & pidx ) );
        end
        cidx                    = ismember( ebinC, bidxx );
        ecounts( cidx, i + 1 )  = cnts;
    end
    
    ibinC                       = ( 1 : max( s.ninhPost( iidx ) ) )';
    nbins                       = length( ibinC );
    icounts                     = zeros( nbins, 2 );
    for i                       = 0 : 1
        switch i
            case 0 
                [ cnts, bidxx ]  = uhist( s.ninhPost( iidx & ~pidx & ~bidx ) );
            case 1
                [ cnts, bidxx ]  = uhist( s.ninhPost( iidx & pidx ) );
        end
        cidx                    = ismember( ibinC, bidxx );
        icounts( cidx, i + 1 )  = cnts;
    end
    % to do - plot these histograms and check whether scaled or different
    %[ ibinC icounts ]
    %[ ebinC ecounts ]
    %[ hh pp ] = kstest2( s.nexcPost( eidx & ~pidx ), s.nexcPost( eidx & pidx ) )
    
    % number of pre-synaptic peers
    npre            = NaN( 3, 6 );
    npre( :, 1 )    = [ mean( s.nexcPre( eidx2 & ~pidx & ~bidx ) ) calc_sem( s.nexcPost( eidx2 & ~pidx & ~bidx ) ) sum( eidx2 & ~pidx & ~bidx ) ]';
    npre( :, 2 )    = [ mean( s.nexcPre( eidx2 & pidx ) ) calc_sem( s.nexcPost( eidx2 & pidx ) ) sum( eidx2 & pidx ) ]';
    npre( :, 3 )    = [ mean( s.ninhPre( iidx2 & ~pidx & ~bidx ) ) calc_sem( s.ninhPre( iidx2 & ~pidx & ~bidx ) ) sum( iidx2 & ~pidx & ~bidx ) ]';
    npre( :, 4 )    = [ mean( s.ninhPre( iidx2 & pidx ) ) calc_sem( s.ninhPre( iidx2 & pidx ) ) sum( iidx2 & pidx ) ]';
    npre( :, 5 )    = [ mean( s.nexcPre( eidx2 & bidx ) ) calc_sem( s.nexcPre( eidx2 & bidx ) ) sum( eidx2 & bidx ) ]';
    npre( :, 6 )    = [ mean( s.ninhPre( iidx2 & bidx ) ) calc_sem( s.ninhPre( iidx2 & bidx ) ) sum( iidx2 & bidx ) ]';
    
    % nunits-E, punits-E, nunits-I, punits-I
    H_JB = [ jbtest( s.nexcPre( eidx2 & ~pidx & ~bidx ) ), jbtest( s.nexcPre( eidx2 & pidx ) ) ...
        jbtest( s.ninhPre( iidx2 & ~pidx & ~bidx ) ), jbtest( s.ninhPre( iidx2 & pidx ) ) ];
    if all( H_JB == 0 )
        [ ~, pp1 ]  = ttest2( s.nexcPre( eidx2 & ~pidx & ~bidx ) , s.nexcPre( eidx2 & pidx ) );
        [ ~, pp2 ]  = ttest2( s.ninhPre( iidx2 & ~pidx & ~bidx ), s.ninhPre( iidx2 & pidx ) );
    else
        pp1         = utest( s.nexcPre( eidx2 & ~pidx & ~bidx ), s.nexcPre( eidx2 & pidx ) );
        pp2         = utest( s.ninhPre( iidx2 & ~pidx & ~bidx ), s.ninhPre( iidx2 & pidx ) );
    end
    pval_NPRE   = [ pp1 pp2 ];
    rez{ 3, 1 }.npre           = npre;
        
    % ASG
    asg             = NaN( 3, 6 );
    asg( :, 1 )    = [ nangeomean( s.asg1( eidx & ~pidx & ~bidx ) ) calc_sem( s.asg1( eidx & ~pidx & ~bidx ) ) sum( eidx & ~pidx & ~bidx ) ]';
    asg( :, 2 )    = [ nangeomean( s.asg1( eidx & pidx ) ) calc_sem( s.asg1( eidx & pidx ) ) sum( eidx & pidx ) ]';
    asg( :, 3 )    = [ -nangeomean( -s.asg2( iidx & ~pidx & ~bidx ) ) calc_sem( s.asg2( iidx & ~pidx & ~bidx ) ) sum( iidx & ~pidx & ~bidx ) ]';
    asg( :, 4 )    = [ -nangeomean( -s.asg2( iidx & pidx ) ) calc_sem( s.asg2( iidx & pidx ) ) sum( iidx & pidx ) ]';
    asg( :, 5 )    = [ nangeomean( s.asg1( eidx & bidx ) ) calc_sem( s.asg1( eidx & bidx ) ) sum( eidx & bidx ) ]';
    asg( :, 6 )    = [ -nangeomean( -s.asg2( iidx & bidx ) ) calc_sem( s.asg2( iidx & bidx ) ) sum( iidx & bidx ) ]';

    % E-nunits, E-punits, I-nunits, I-punits
    pval_ASG        = [ utest( s.asg1( eidx & ~pidx ), s.asg1( eidx & pidx ) )
        utest( s.asg2( iidx & ~pidx ), s.asg2( iidx & pidx ) ) ];
    rez{ 3, 1 }.asg           = asg;

    fig3( 1 ) = figure;
    subplot( 3, 3, 1 )
    [ bh, eh ] = barwerror( 1 : 6, [ emat( 3, [ 2 1 3 ] ) imat( 3, [ 2 1 3 ] ) ], [ emat( 4, [ 2 1 3 ] ) imat( 4, [ 2 1 3 ] )] ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'E-nunits', 'E-punits','E-bipolars', 'I-nunits', 'I-punits','I-bipolars' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Units that are excitatory (E) or inhibitory (I)' )
    
    subplot( 3, 3, 2 )
    [ bh, eh ] = barwerror( 1 : 6, [ emat2( 3, [ 2 1 3 ] ) imat2( 3, [ 2 1 3 ] ) ] ...
        , [ emat2( 4, [ 2 1 3 ] ) imat2( 4, [ 2 1 3 ] ) ] ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Excited-bipolars', 'Inhibited-nunits', 'Inhibited-punits' , 'Inhibited-bipolars' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Units that receive excitation or inhibition' )
    
    subplot( 3, 3, 3 )
    [ bh, eh ] = barwerror( 1 : 3, cmat( 3, [ 2 1 3 ] ), cmat( 4, [ 2 1 3 ] ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Connected-nunits', 'Connected-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Connected units' )
    
    subplot( 3, 3, 4 )
    [ bh, eh ] = barwerror( 1 : 6, npost( 1, : ), npost( 2, : ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'E-nunits', 'E-punits', 'E-bipolars', 'I-nunits', 'I-punits', 'I-bipolars'};
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number of units' )
    title( 'Post-synaptic peers for E and I units' )

    subplot( 3, 3, 5 )
    [ bh, eh ] = barwerror( 1 : 6, npre( 1, : ), npre( 2, : ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Excited-bipolars', 'Inhibited-nunits', 'Inhibited-punits', 'Inhibited-bipolars' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number units' )
    title( 'Pre-synaptic peers for E (excitatory) and I (inhibitory) units' )
    
    subplot( 3, 3, 6 )
    [ bh, eh ] = barwerror( 1 : 6, abs( asg( 1, : ) ), asg( 2, : ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'E-nunits', 'E-punits', 'E-bipolars', 'I-nunits', 'I-punits', 'I-bipolars' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for E and I units' )

    subplot( 3, 3, 7 )
    [ bh, eh ] = barwerror( 1 : 3, asg( 1, [1 : 2 , 5] ), asg( 2, [1 : 2 , 5] ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'E-nunits', 'E-punits', 'E-bipolars' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for E units' )

    subplot( 3, 3, 8 )
    [ bh, eh ] = barwerror( 1 : 3, asg( 1, [3 : 4 , 6] ), asg( 2, [3 : 4 , 6] ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'I-nunits', 'I-punits', 'I-bipolars'};
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for I units' )
    
    if ~isempty( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) );
    end
    
    %-----
    % example of punit and CCH and ACH
    if 0 % showExample
        if ismac
            filebase = filebaseLookup( 'mDL5', -16 );
        elseif ispc
            filebase = 'G:\mice\mDL5_16\mDL5_16';
        elseif isunix
            filebase = '/media/shirly/C22865A128659567/mice/mDL5/dat/mDL5_16/mDL5_16';
        end
        fig3( 2 ) = figure;
        plot_ss( filebase, [ 1 6 ] ); % punit with 5 post-synaptic peers and no pre
        %mono.pairsExc( mono.pairsExc( :, 1 ) == 3, : )
        % punit is on S1 (1.6)
        % show CCH with a post-synaptic INT (S2; 2.11), and a post-synaptic PYR (S3; 3.52)
        % note that [ 3 50 ] is also a punit, with multiple pre-synaptic peers
        fig3( 3 ) = figure;
        plot_ss( filebase, [ 2 11 ] );
        fig3( 4 ) = figure;
        plot_ss( filebase, [ 3 52 ] );
        fig3( 5 ) = figure;
        subplot( 2, 2, 1 )
        plot_s2s( filebase, [ 1 6; 2 11 ], 'plotmode', -4 );
        subplot( 2, 2, 2 )
        plot_s2s( filebase, [ 1 6; 2 12 ], 'plotmode', -4 );
        subplot( 2, 2, 3 )
        plot_s2s( filebase, [ 1 6; 3 52 ], 'plotmode', -4 );
        
        % quantify the strength
        load( [ filebase '.s2s' ], '-mat', 's2s' )
        n12                         = [ 1 6; 2 11 ];
        [ g1, g2, act, sil, s5, fig ] = calc_asg( s2s, n12, 'graphics', 1 );
        fig3( 6 ) = fig;
        
        n12                         = [ 1 6; 2 12 ];
        [ g1, g2, act, sil, s6, fig ] = calc_asg( s2s, n12, 'graphics', 1 );
        fig3( 7 ) = fig;
        
        n12                         = [ 1 6; 3 52 ];
        [ g1, g2, act, sil, s7, fig ] = calc_asg( s2s, n12, 'graphics', 1 );
        fig3( 8 ) = fig;
    end
    
    %--------------------------------------------------------------------
    % focused plots + statistical analyses
    alfa                            = 0.001; % from s2s
    fig3( 2 ) = figure;
    
    subplot( 2, 2, 1 )
    hold on
    nn                              = emat( 1, [ 2 1 3 ] );
    tt                              = emat( 2, [ 2 1 3 ] );
    xx                              = emat( 3, [ 2 1 3 ] );
    ss                              = emat( 4, [ 2 1 3 ] );
    pp                              = binp( emat( 1, [ 2 1 3 ] ), emat( 2, [ 2 1 3 ] ), alfa );
    pval                            = gtest( [ tt; nn ], 'lr', 1, 'ind' ); % test of independence
    rez{ 3, 2 }.nn(1,:)               = nn;
    rez{ 3, 2 }.tt(1,:)               = tt;
    rez{ 3, 2 }.xx(1,:)               = xx;
    rez{ 3, 2 }.ss(1,:)               = ss;
    rez{ 3, 2 }.pp(1,:)               = pp;
    rez{ 3, 2 }.pval(1)             = pval;

    for i = 1 : 3
        barwerror( i, xx( i ), ss( i ), colors_NPB( i, : ), 0.8, [ 0 0 0 ], 2 );
        str = sprintf( '%0.2g%% (%d/%d)\n%0.3g', nn( i ) / tt( i ) * 100, nn( i ), tt( i ), pp( i ) );
        th = text( i, max( xx + ss ) * 1.1, str );
        set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
    end
    ylim( [ 0 max( xx + ss ) * 1.3 ] )
    colnames                        = { 'E-nunits', 'E-punits', 'E-bipolars' };
    set( gca, 'xtick', 1 : 3, 'XTickLabel', colnames )
    xlim( [ 1 3 ] + [ -1 1 ] * 0.8 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( sprintf( 'Units that are excitatory (E), p=%0.3g', pval ) )
    alines( alfa, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 2, 2 )
    hold on
    nn                              = imat( 1, [ 2 1 3 ] );
    tt                              = imat( 2, [ 2 1 3 ] );
    xx                              = imat( 3, [ 2 1 3 ] );
    ss                              = imat( 4, [ 2 1 3 ] );
    pp                              = binp( imat( 1, [ 2 1 3 ] ), imat( 2, [ 2 1 3 ] ), alfa );
    pval                            = gtest( [ tt; nn ], 'lr', 1, 'ind' ); % test of independence
    
    rez{ 3, 2 }.nn(2,:)               = nn;
    rez{ 3, 2 }.tt(2,:)               = tt;
    rez{ 3, 2 }.xx(2,:)               = xx;
    rez{ 3, 2 }.ss(2,:)               = ss;
    rez{ 3, 2 }.pp(2,:)               = pp;
    rez{ 3, 2 }.pval(2)             = pval;

    for i = 1 : 3
        barwerror( i, xx( i ), ss( i ), colors_NPB( i, : ), 0.8, [ 0 0 0 ], 2 );
        str = sprintf( '%0.2g%% (%d/%d)\n%0.3g', nn( i ) / tt( i ) * 100, nn( i ), tt( i ), pp( i ) );
        th = text( i, max( xx + ss ) * 1.1, str );
        set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
    end
    ylim( [ 0 max( xx + ss ) * 1.3 ] )
    colnames = { 'I-nunits', 'I-punits', 'I-bipolars' };
    set( gca, 'xtick', 1 : 3, 'XTickLabel', colnames )
    xlim( [ 1 3 ] + [ -1 1 ] * 0.8 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( sprintf( 'Units that are inhibitory (I), p=%0.3g', pval ) )
    alines( alfa, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    if ~isempty( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) );
    end
    
    %-----
    % save the figures
    fig = fig3;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG3_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
end 
if fignums( 3 ) && metaplot

        fig3m( 1 )           = figure;
                   
                subplot( 2, 2, 1 )
                [ bh, eh ] = barwerror( 1 : 4, [rez{1,1}{3,1}.emat( 3, [ 2 1 ] ) rez{2,1}{3,1}.emat( 3, 1 ) ...
                    rez{3,1}{3,1}.emat( 3, 1 )] , [rez{1,1}{3,1}.emat( 4, [ 2 1 ] ) rez{2,1}{3,1}.emat( 4, 1  )...
                    rez{3,1}{3,1}.emat( 4, 1 )] , colors_NP( 1, : ), 0.8, [0 0 0], 1 );
                colnames = { 'Nunits', 'P1', 'P2', 'P3' };
                set( gca, 'XTickLabel', colnames )
                set( gca, 'tickdir', 'out', 'box', 'off' )
                ylabel( 'Fraction of units' )
                title( 'Excitatory units' )
                str = sprintf( '%0.2g%% (%d/%d)\n%0.3g', rez{1,1}{3,2}.nn( 1,1 ) / rez{1,1}{3,2}.tt( 1,1 ) * 100, rez{1,1}{3,2}.nn( 1,1 ), rez{1,1}{3,2}.tt( 1,1 ), rez{1,1}{3,2}.pp( 1,1 ) );
                th = text( 1, max( rez{1,1}{3,2}.xx(1) + rez{1,1}{3,2}.ss(1) ), str, 'FontSize',6 );
                set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
                for pi = 1:3  
                    str = sprintf( '%0.2g%% (%d/%d)\n%0.3g', rez{pi,1}{3,2}.nn(1, 2 ) / rez{pi,1}{3,2}.tt(1, 2 ) * 100, rez{pi,1}{3,2}.nn( 1,2 ), rez{pi,1}{3,2}.tt( 1,2 ), rez{pi,1}{3,2}.pp( 1,2 ) );
                    th = text( 1+pi, max( rez{pi,1}{3,2}.xx(1,2) + rez{pi,1}{3,2}.ss(1,2) ), str, 'FontSize',6 );
                    set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
                end 
                
                subplot( 2, 2, 2 )
                [ bh, eh ] = barwerror( 1 : 4, [rez{1,1}{3,1}.imat( 3, [ 2 1 ] ) rez{2,1}{3,1}.imat( 3, 1 ) ...
                    rez{3,1}{3,1}.imat( 3, 1 )] , [rez{1,1}{3,1}.imat( 4, [ 2 1 ] ) rez{2,1}{3,1}.imat( 4, 1  )...
                    rez{3,1}{3,1}.imat( 4, 1 )] , colors_NP( 1, : ), 0.8, [0 0 0], 1 );
                set( gca, 'XTickLabel', colnames )
                set( gca, 'tickdir', 'out', 'box', 'off' )
                ylabel( 'Fraction of units' )
                title( 'Inhibitory units' )
                str = sprintf( '%0.2g%% (%d/%d)\n%0.3g', rez{1,1}{3,2}.nn( 2,1 ) / rez{1,1}{3,2}.tt( 2,1 ) * 100, rez{1,1}{3,2}.nn( 2,1 ), rez{1,1}{3,2}.tt( 2,1 ), rez{1,1}{3,2}.pp( 2,1 ) );
                th = text( 1, max( rez{1,1}{3,2}.xx(2,1) + rez{1,1}{3,2}.ss(2,1) ), str, 'FontSize',6 );
                set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
                for pi = 1:3  
                    str = sprintf( '%0.2g%% (%d/%d)\n%0.3g', rez{pi,1}{3,2}.nn(2, 2 ) / rez{pi,1}{3,2}.tt(2, 2 ) * 100, rez{pi,1}{3,2}.nn( 2,2 ), rez{pi,1}{3,2}.tt( 2,2 ), rez{pi,1}{3,2}.pp( 2,2 ) );
                    th = text( 1+pi, max( rez{pi,1}{3,2}.xx(2,2) + rez{pi,1}{3,2}.ss(2,2) ), str, 'FontSize',6 );
                    set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )
                end 
                
                subplot( 2, 2, 3 )
                [ bh, eh ] = barwerror( 1 : 4, [ rez{1,1}{3,1}.emat2( 3, [ 2 1 ] ) rez{2,1}{3,1}.emat2( 3, 1 ) ...
                    rez{3,1}{3,1}.emat2( 3, 1 )] , [rez{1,1}{3,1}.emat2( 4, [ 2 1 ] ) rez{2,1}{3,1}.emat2( 4, 1  )...
                    rez{3,1}{3,1}.emat2( 4, 1 )] , colors_NP( 1, : ), 0.8, [0 0 0], 1 );
                set( gca, 'XTickLabel', colnames )
                set( gca, 'tickdir', 'out', 'box', 'off' )
                ylabel( 'Fraction of units' )
                title( 'Units that receive excitation' )
                                
                subplot( 2, 2, 4 )
                [ bh, eh ] = barwerror( 1 : 4, [  rez{1,1}{3,1}.imat2( 3, [ 2 1 ] ) rez{2,1}{3,1}.imat2( 3, 1 ) ...
                    rez{3,1}{3,1}.imat2( 3, 1 )] , [rez{1,1}{3,1}.imat2( 4, [ 2 1 ] ) rez{2,1}{3,1}.imat2( 4, 1  )...
                    rez{3,1}{3,1}.imat2( 4, 1 )] , colors_NP( 1, : ), 0.8, [0 0 0], 1 );
                set( gca, 'XTickLabel', colnames )
                set( gca, 'tickdir', 'out', 'box', 'off' )
                ylabel( 'Fraction of units' )
                title( 'Units that receive inhibition' )
                
          fig3m( 2 )           = figure;

            [ bh, eh ] = barwerror( 1 : 4, [  rez{1,1}{3,1}.cmat( 3, [ 2 1 ] ) rez{2,1}{3,1}.cmat( 3, 1 ) ...
                rez{3,1}{3,1}.cmat( 3, 1 ) ], [rez{1,1}{3,1}.cmat( 4, [2 1] ) rez{2,1}{3,1}.cmat( 4, 1 ) ...
                rez{3,1}{3,1}.cmat( 4, 1 )] , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 1 );
                set( gca, 'XTickLabel', colnames )
                set( gca, 'tickdir', 'out', 'box', 'off' )
                ylabel( 'Fraction of units' )
                title( 'Connected units' )
                
        fig = fig3m;
        if savef
            for i = 1 : length( fig )
                if fig( i ) == 0 || isempty( fig( i ) )
                    continue
                end
                figname = sprintf( '%s/punits_FIG3m_part%d%s', outdir, i, clunumstr );
                figure( fig( i ) );
                figi = gcf;
                figi.Renderer = 'painters';
                pause( 0.2 )
                print( fig( i ), pstr, figname, resize )

            end
        end
    
end

% ----------------------------------------------------------------------
% Figure 4
% firing rate statistics (fano factor)
% ----------------------------------------------------------------------

if fignums( 4 ) && ~metaplot
    
    % intersect the units in the sst and the ff structures
    u1                      = unique( sst.filebase );
    u2                      = unique( ff.filebase );
    u                       = unique( [ u1; u2 ] );
    [ ~, f1 ]               = ismember( sst.filebase, u );
    [ ~, f2 ]               = ismember( ff.filebase, u );
    s1                      = [ f1 sst.shankclu ];
    s2                      = [ f2 ff.shankclu( :, 1 : 2 ) ];
    % since each item appears only once, we can use intersect
    [ ~, i1, i2 ]           = intersect( s1, s2, 'rows' );
    i1l                     = false( size( sst.shankclu, 1 ), 1 );
    i1l( i1 )               = 1;
    i2l                     = false( size( ff.shankclu, 1 ), 1 );
    i2l( i2 )               = 1;
    sst1                    = struct_select( sst, i1l );
    ff1                     = struct_select( ff, i2l );
    % merge the two structures
    sst1                    = struct_merge( rmfield( sst1, 'shankclu' ), rmfield( ff1, 'filebase' ) );
    
    % remove units with negative FF??
    ridx                    = sum( sst1.fanomat < 0, 2 );
    sst1                    = struct_select( sst1, ~ridx );
    
    rez{ 4, 1 }.sst1                 = sst1;

    % process the data
    
    
    % plot the figures
    %--------------------------------
    % ACH COM
    fig4( 1 ) = figure; 
    
%     ispos                   = sst1.extremum > 0;
    nbins                   = 60;
    isbip           = ~isnan(sst1.bpi);
    ispos           = sst1.extremum > 0 & ~isbip;
    
    rez{ 4, 1 }.ispos                 = ispos;
    rez{ 4, 1 }.nbins                 = nbins;

    byprob                  = 0;
    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ach_com, ispos, nbins, colors_NP, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'Punits'; 
    lh.String{ 2 }          = 'Nunits';

    subplot( 2, 2, 2 )
    sst2                    = struct_select( sst1, ~ispos & ~isbip );
    rez{ 4, 1 }.sstN                 = sst2;

    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst2.ach_com, sst2.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    nbins                   = 30;
    byprob                  = 1;
    subplot( 2, 2, 3 )
    sst1i                   = struct_select( sst1, ispos );
    rez{ 4, 1 }.sstP                 = sst1i;

    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1i.ach_com, true( sum( ispos ), 1 ), nbins, colors_NP, 'ACH-COM [ms]', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst2.ach_com, sst2.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    % statistics - KW and multiple comparisons
    Xkw                     = sst1.ach_com;
    Gkw                     = sst1.shankclu( :, 3 );
    Gkw( ispos )            = 2;
    [ Pkw, ~, Skw ]         = kruskalwallis( Xkw, Gkw,'off' );
    Tmc                     = multcompare( Skw, 'displayopt', 'off' );
    pvals                   = Tmc( :, 6 );
    cnames                  = { 'INT-PYR', 'INT-Punits', 'PYR-Punits' };

    rez{ 4, 1 }.cnames      = cnames;
    rez{ 4, 1 }.pvals       = pvals;
        
    fig4( 2 )               = figure;

    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ach_com, isbip, nbins, colors_NPB ([1 3],:), 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'Bipolars'; 
    lh.String{ 2 }          = 'Nunits';

    subplot( 2, 2, 2 )
    sst2                    = struct_select( sst1, ~ispos & ~isbip);
    rez{ 4, 1 }.sstN                 = sst2;

    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst2.ach_com, sst2.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    nbins                   = 30;
    byprob                  = 1;
    subplot( 2, 2, 3 )
    sst1i                   = struct_select( sst1, ispos );
    rez{ 4, 1 }.sstP                 = sst1i;

    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1i.ach_com, true( sum( isbip ), 1 ), nbins, colors_NPB([1 3],:), 'ACH-COM [ms]', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst2.ach_com, sst2.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';

    title( sprintf( '%s: %0.2g; %s: %0.2g; %s: %0.2g', cnames{ 1 }, pvals( 1 ) ...
        , cnames{ 2 }, pvals( 2 ), cnames{ 3 }, pvals( 3 ) ) )
    if ~isempty( clunum ), fig_title( sprintf( 'cluster %d', clunum ) ); end

    %--------------------------------
    % FF - baseline rates
    fig4( 2 )               = figure;
    ispos                   = sst1.extremum > 0;
    nbins                   = 60;
    byprob                  = 0;
    logsst_BL               = log2( sst1.baseline );

    rez{ 4, 2 }.ispos       = ispos;
    rez{ 4, 2 }.nbins       = nbins;
    rez{ 4, 2 }.logsst_BL   = logsst_BL;

    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( logsst_BL, ispos, nbins, colors_NP, 'Log Firing rate [spk/s]', byprob );
    lh.String{ 1 }          = 'Punits'; 
    lh.String{ 2 }          = 'Nunits';
    lin2log( 'x', 2, 5 );

    subplot( 2, 2, 2 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( logsst_BL, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'Log Firing rate [spk/s]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    lin2log( 'x', 2, 5 );
    
    nbins                   = 30;
    byprob                  = 1;
    subplot( 2, 2, 3 )
    sst1i                   = struct_select( sst1, ispos );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( log2( sst1i.baseline ), true( sum( ispos ), 1 ), nbins, colors_NP, 'Log Firing rate [spk/s]', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( log2( sst1.baseline ), sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'Log Firing rate [spk/s]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    lin2log( 'x', 2, 5 );
    
    % statistics - KW and multiple comparisons
    Xkw                     = sst1.ach_com;
    Gkw                     = sst1.shankclu( :, 3 );
    Gkw( ispos )            = 2;
    [ Pkw, ~, Skw ]         = kruskalwallis( Xkw, Gkw,'off' );
    Tmc                     = multcompare( Skw, 'displayopt', 'off' );
    pvals                   = Tmc( :, 6 );
    cnames                  = { 'INT-PYR', 'INT-Punits', 'PYR-Punits' };
    title( sprintf( '%s: %0.2g; %s: %0.2g; %s: %0.2g', cnames{ 1 }, pvals( 1 ) ...
        , cnames{ 2 }, pvals( 2 ), cnames{ 3 }, pvals( 3 ) ) )
    if ~isempty( clunum ), fig_title( sprintf( 'cluster %d', clunum ) ); end
    
    rez{ 4, 2 }.pvals       = pvals;
    rez{ 4, 2 }.cnames      = cnames;
    
    %--------------------------------
    % FF 
    % to do - look at individual images
    fig4( 3 )                   = figure;
    hold on
    for i                       = 0 : 2
        
        idx                     = Gkw == i;
        y                       = sst1.fanomat( idx, : );
        y                       = outliers( y, 3, [], 1, 2 );       % remove outliers
        x                       = sst1.winsizes( 1, : );
        if i < 2
            c1                  = colors_PI( i + 1, : );
        else
            c1                  = colors_NP( 2, : );
        end
        c2                      = [];
        
        %plot( log2( x ), nangeomean( y, 1 ), 'b' );
        ph                      = patch_band( log2( x / 20 ), nangeomean( y, 1 ), calc_sem( y, 1 ), c1, c2 );
        
    end
    for i                       = 0 : 2
        
        idx                     = Gkw == i;
        y                       = sst1.fanomat( idx, : );
        y                       = outliers( y, 3, [], 1, 2 );       % remove outliers
        x                       = sst1.winsizes( 1, : );
        if i < 2
            c1                  = colors_PI( i + 1, : );
        else
            c1                  = colors_NP( 2, : );
        end
        c2                      = [];
        
        ph                      = plot( log2( x / 20 ), nangeomean( y, 1 ), '.-' );
        set( ph, 'MarkerSize', 20, 'color', c1 )
        
    end
    lin2log( 'x', 2, 1 );
    
    xlim( log2( x( [ 1 end ] ) / 20 ) )
    xlabel( 'Window [ms]' )
    
    ylims = ylim; 
    set( gca, 'yscale', 'log' ), 
    ylim( [ 0.95 ylims( 2 ) ] )
    ylabel( 'Fano factor' )

    alines( 1, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    set( gca, 'tickdir', 'out', 'box', 'off' )

    if ~isempty( clunum ), fig_title( sprintf( 'cluster %d', clunum ) ); end
    
    
    %--------------------------------
    % make a cut through the FF figure and present histograms of the FF at individual windows
    colors                      = [ colors_PI; colors_NP( 2, : ) ];

    nbins                       = 30;
    byprob                      = 1;
    
    winsizes                    = sst1.winsizes( 1, : );
    nwinsizes                   = length( winsizes );
    hh                          = NaN( nbins, 3, nwinsizes );
    binC                        = NaN( nbins, nwinsizes );
    edgesL                      = NaN( nbins, nwinsizes );
    for wi                      = 1 : nwinsizes
        winsize                 = winsizes( wi );
        mm                      = log2( minmax( sst1.fanomat( :, winsizes == winsize ) ) );
        mm                      = [ floor( mm( 1 ) ) ceil( mm( 2 ) ) ];
        binsize                 = diff( mm ) / nbins;
        %edges                   = mm( 1 ) : diff( mm ) / nbins : mm( 2 );
        edges                   = mm( 1 ) - binsize / 2 : binsize : mm( 2 );
        edgesL( :, wi )         = 2 .^ edges( 1 : end - 1 );
        binC( :, wi )           = 2 .^ ( ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2 );
        edges                   = 2 .^ edges;
        for i                   = 0 : 2
            idx                 = Gkw == i;
            y                   = sst1.fanomat( idx, : );
            y                   = outliers( y, 3, [], 1, 2 );       % remove outliers
            y                   = y( :, winsizes == winsize );
            h                   = histc( y, edges );
            h( end - 1 )        = h( end - 1 ) + h( end );
            h( end )            = [];
            hh( :, i + 1, wi )  = h;
        end
    end
    
    
    fig4( 4 ) = figure;
    for wi                      = 1 : nwinsizes
        subplot( 3, 4, wi )
        hold on
        for i = 0 : 2
            hi                  = hh( :, i + 1, wi );
            if byprob
                hi              = hi / sum( hi );
            end
            %bh( i + 1 )         = stairs( log2( binC( :, wi ) ), hi, 'color', colors( i + 1, : ), 'linewidth', 2 );
            %bh( i + 1 )         = stairs( log2( binC( :, wi ) - diff( log2( binC( 1 : 2, wi ) ) ) / 2 ), hi, 'color', colors( i + 1, : ), 'linewidth', 2 );
            bh( i + 1 )         = stairs( log2( edgesL( :, wi ) ), hi, 'color', colors( i + 1, : ), 'linewidth', 1 );
        end
        lin2log( 'x', 2, 1 );
        if wi == nwinsizes
            lh                  = legend( bh, { 'INT', 'PYR', 'Punit' } );
        end
        set( gca, 'tickdir', 'out', 'box', 'off' );
        title( sprintf( '%0.3g ms', winsizes( wi ) / 20 ) )
    end
    if ~isempty( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) ); 
    end
    

    % LR/ID parameter
    
    fig4( 5 ) = figure; 

    ispos                   = sst1.extremum > 0;
    nbins                   = 60;
    byprob                  = 0;
    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.Lratio, ispos, nbins, colors_NP, 'Lratio', byprob );
    lh.String{ 1 }          = 'Punits'; 
    lh.String{ 2 }          = 'Nunits';

    subplot( 2, 2, 2 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.Lratio, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'Lratio', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    nbins                   = 30;
    byprob                  = 1;
    subplot( 2, 2, 3 )
    sst1i                   = struct_select( sst1, ispos );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1i.Lratio, true( sum( ispos ), 1 ), nbins, colors_NP, 'Lratio', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.Lratio, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'Lratio', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    % statistics - KW and multiple comparisons
    Xkw                     = sst1.Lratio;
    Gkw                     = sst1.shankclu( :, 3 );
    Gkw( ispos )            = 2;
    [ Pkw, ~, Skw ]         = kruskalwallis( Xkw, Gkw,'off' );
    Tmc                     = multcompare( Skw, 'displayopt', 'off' );
    pvals                   = Tmc( :, 6 );
    cnames                  = { 'INT-PYR', 'INT-Punits', 'PYR-Punits' };
    title( sprintf( '%s: %0.2g; %s: %0.2g; %s: %0.2g', cnames{ 1 }, pvals( 1 ) ...
        , cnames{ 2 }, pvals( 2 ), cnames{ 3 }, pvals( 3 ) ) )

    if ~isempty( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) ); 
    end
    
    
    % ISI-index parameter

    fig4( 6 ) = figure; 

    ispos                   = sst1.extremum > 0;
    nbins                   = 60;
    byprob                  = 0;
    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ISIindex, ispos, nbins, colors_NP, 'ISIindex', byprob );
    lh.String{ 1 }          = 'Punits'; 
    lh.String{ 2 }          = 'Nunits';

    subplot( 2, 2, 2 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ISIindex, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ISIindex', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    nbins                   = 30;
    byprob                  = 1;
    subplot( 2, 2, 3 )
    sst1i                   = struct_select( sst1, ispos );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1i.ISIindex, true( sum( ispos ), 1 ), nbins, colors_NP, 'ISIindex', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ISIindex, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ISIindex', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    % statistics - KW and multiple comparisons
    Xkw                     = sst1.ISIindex;
    Gkw                     = sst1.shankclu( :, 3 );
    Gkw( ispos )            = 2;
    [ Pkw, ~, Skw ]         = kruskalwallis( Xkw, Gkw,'off' );
    Tmc                     = multcompare( Skw, 'displayopt', 'off' );
    pvals                   = Tmc( :, 6 );
    cnames                  = { 'INT-PYR', 'INT-Punits', 'PYR-Punits' };
    title( sprintf( '%s: %0.2g; %s: %0.2g; %s: %0.2g', cnames{ 1 }, pvals( 1 ) ...
        , cnames{ 2 }, pvals( 2 ), cnames{ 3 }, pvals( 3 ) ) )

    if ~isempty( clunum )
        fig_title( sprintf( 'cluster %d', clunum ) ); 
    end
    
    
    %-----
    % save the figures
    fig = fig4;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG4_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

if fignums( 4 ) && metaplot
    
    fig4m( 1 )           = figure;
    
    a1 = rez{1,1}{4,1}.sst1.ach_com;
    idx1 = rez{1,1}{4,1}.ispos;
    a2 = rez{2,1}{4,1}.sst1.ach_com;
    idx2 = rez{2,1}{4,1}.ispos;
    a3 = rez{3,1}{4,1}.sst1.ach_com;
    idx3 = rez{3,1}{4,1}.ispos;
    pyridx = rez{1,1}{4,1}.sst1.pyr;

    ach_com_punits1 = a1(idx1);
    ach_com_punits2 = a2(idx2);
    ach_com_punits3 = a3(idx3);
    ach_com_nunits = a1(~idx1);
    ach_com_nunitsPYR = a1(~idx1&pyridx);
    ach_com_nunitsINT = a1(~idx1&~pyridx);

    [ y1 x1 ] = cdfcalc( ach_com_punits1);
    x                               = [ x1 ]; 
    y                               = y1(2:end);
    xq                              = x( 1 ) : 0.1 : x( end ); 
    yq                              = interp1( x, y, xq );
    med00i                          = x( find( y >= 0.5, 1, 'first' ) );
    xq00                            = x;
    yq00                            = y;
    ph                              = plot( xq00, yq00, '-','color',colors_sub(1,:));
    alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
    alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
    set( gca, 'tickdir', 'out', 'box', 'off' )
    text(med00i,1,sprintf('%0.2g',med00i));
    axis square
    hold on

    [ y1 x1 ] = cdfcalc( ach_com_punits2);
    x                               = [ x1 ]; 
    y                               = y1(2:end);
    xq                              = x( 1 ) : 0.1 : x( end ); 
    yq                              = interp1( x, y, xq );
    med00i                          = x( find( y >= 0.5, 1, 'first' ) );
    xq00                            = x;
    yq00                            = y;
    ph                              = plot( xq00, yq00, '-','color',colors_sub(2,:));
    alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
    alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
    text(med00i,1,sprintf('%0.2g',med00i));


    [ y1 x1 ] = cdfcalc( ach_com_punits3);
    x                               = [ x1 ]; 
    y                               = y1(2:end);
    xq                              = x( 1 ) : 0.1 : x( end ); 
    yq                              = interp1( x, y, xq );
    med00i                          = x( find( y >= 0.5, 1, 'first' ) );
    xq00                            = x;
    yq00                            = y;
    ph                              = plot( xq00, yq00, '-','color',colors_sub(3,:));
    alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
    alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
    text(med00i,1,sprintf('%0.2g',med00i));

    [ y1 x1 ] = cdfcalc( ach_com_nunitsPYR);
    x                               = [ x1 ]; 
    y                               = y1(2:end);
    xq                              = x( 1 ) : 0.1 : x( end ); 
    yq                              = interp1( x, y, xq );
    med00i                          = x( find( y >= 0.5, 1, 'first' ) );
    xq00                            = x;
    yq00                            = y;
    ph                              = plot( xq00, yq00, '-','color',colors_PI(2,:));
    alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
    alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
    text(med00i,1,sprintf('%0.2g',med00i));

    [ y1 x1 ] = cdfcalc( ach_com_nunitsINT);
    x                               = [ x1 ]; 
    y                               = y1(2:end);
    xq                              = x( 1 ) : 0.1 : x( end ); 
    yq                              = interp1( x, y, xq );
    med00i                          = x( find( y >= 0.5, 1, 'first' ) );
    xq00                            = x;
    yq00                            = y;
    ph                              = plot( xq00, yq00, '-','color',colors_PI(1,:));
    alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
    alines( 0.5, 'y', 'color','k', 'linestyle', '--' );
    text(med00i,1,sprintf('%0.2g',med00i));

    % [ y1 x1 ] = cdfcalc( ach_com_nunits);
    % x                               = [ x1 ]; 
    % y                               = y1(2:end);
    % xq                              = x( 1 ) : 0.1 : x( end ); 
    % yq                              = interp1( x, y, xq );
    % med00i                          = x( find( y >= 0.5, 1, 'first' ) );
    % xq00                            = x;
    % yq00                            = y;
    % ph                              = plot( xq00, yq00, '-');
    % alines( med00i, 'x', 'color','k' , 'linestyle', '--' );
    % alines( 0.5, 'y', 'color','k', 'linestyle', '--' );

    ylim([0 1.1])
    xlim([10 35])
    ylabel('Columative distribution')
    xlabel('ACH - COM [ms]')
    legend

    byprob                  = 0;
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( rez{1,1}{4,1}.sst1.ach_com, rez{1,1}{4,1}.ispos, 100, [colors_NP_light(1,:); colors_sub(1,:)], 'ACH-COM [ms]', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( rez{2,1}{4,1}.sst1.ach_com, true( sum(rez{2,1}{4,1}.ispos),1), 30, [colors_NP(1,:); colors_sub(2,:)], 'ACH-COM [ms]', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( rez{3,1}{4,1}.sst1.ach_com,  true( sum(rez{3,1}{4,1}.ispos),1), 100, [colors_NP(1,:); colors_sub(3,:)], 'ACH-COM [ms]', byprob );
  
    title( sprintf( '%s: %0.2g; %s: %0.2g; %s: %0.2g', cnames{ 1 }, pvals( 1 ) ...
        , cnames{ 2 }, pvals( 2 ), cnames{ 3 }, pvals( 3 ) ) )
    
    fig = fig4m;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG4m_part%d%s', outdir, i, clunumstr );
            figure( fig( i ) );
            figi = gcf;
            figi.Renderer = 'painters';
            pause( 0.2 )
            print( fig( i ), pstr, figname, resize )

        end
    end
end
% ----------------------------------------------------------------------
% Figure 5
% Place analyses
% ----------------------------------------------------------------------

if fignums( 5 )

%     % load Hadas's data
%     L                       = load( [ datadir_hadas fname_hadas '.mat' ], '-mat' );
%     
%     % partition into two separate structures
%     fieldnames              = fields( L );
%     nfields                 = length( fieldnames );
%     n_all                   = size( getfield( L.uD, 'session_tag' ), 1 );
%     n_fields                = size( getfield( L.fD, 'session_tag' ), 1 );
%     for i                   = 1 : nfields
%         
%         nrow                = size( getfield( L, fieldnames{ i } ), 1 );
%         if nrow == n_all
%             eval( sprintf( 's_all.%s = L.%s;', fieldnames{ i }, fieldnames{ i } ) )
%         elseif nrow == n_fields
%             eval( sprintf( 's_field.%s = L.%s;', fieldnames{ i }, fieldnames{ i } ) )
%         end
%            
%     end
%     
%     % add the extremum field from sst to the relevant items in s_all
%     n_all                   = size (L.uD.shankclu, 1);
%     L.uD.extremum           = NaN( n_all, 1 );
%     L.uD.bpi                = NaN( n_all, 1 );
%     L.uD.upol               = NaN( n_all, 1 );
%     u1                      = unique( sst.filebase );
%     u2                      = unique( L.uD.session_tag );
%     u                       = unique( [ u1; u2 ] );
%     [ ~, f1 ]               = ismember( sst.filebase, u );
%     [ ~, f2 ]               = ismember( L.uD.session_tag, u );
%     s1                      = [ f1 sst.shankclu ];
%     s2                      = [ f2 L.uD.shankclu( :, 1 : 2 ) ];
%     % if Hadas's data were unique, then we could write:
%     %     [ ~, i1, i2 ]           = intersect( s1, s2, 'rows' );
%     %     s_all.extremum( i2, : ) = sst.extremum( i1, : );
%     % since it is not unique, we use the following:
%     us2                     = unique( s2, 'rows' );
%     nus2                    = size( us2, 1 );
%     for i                   = 1 : nus2
%         i1                  = ismember( s1, us2( i, : ), 'rows' );
%         i2                  = ismember( s2, us2( i, : ), 'rows' );
%         if sum( i1 ) > 0
%             L.uD.extremum( i2, : ) = repmat( sst.extremum( i1, : ), [ sum( i2 ) 1 ] );
%             L.uD.bpi( i2, : )      = repmat( sst.bpi( i1, : ), [ sum( i2 ) 1 ] );
%             L.uD.upol( i2, : )     = repmat( sst.upol( i1, : ), [ sum( i2 ) 1 ] );
%         end
%     end
%     
%     % do the same for the s_field:
%     n_fields                = size (L.fD.shankclu, 1);
%     L.fD.extremum           = NaN( n_fields, 1 );
%     L.fD.bpi                = NaN( n_fields, 1 );
%     L.fD.upol               = NaN( n_fields, 1 );
%     u2                      = unique( L.fD.session_tag );
%     u                       = unique( [ u1; u2 ] );
%     [ ~, f2 ]               = ismember( L.fD.session_tag, u );
%     s2                      = [ f2 L.fD.shankclu( :, 1 : 2 ) ];
%     us2                     = unique( s2, 'rows' );
%     nus2                    = size( us2, 1 );
%     for i                   = 1 : nus2
%         i1                  = ismember( s1, us2( i, : ), 'rows' );
%         i2                  = ismember( s2, us2( i, : ), 'rows' );
%         if sum( i1 ) > 0
%             L.fD.extremum( i2, : ) = repmat( sst.extremum( i1, : ), [ sum( i2 ) 1 ] );
%             L.fD.bpi( i2, : )      = repmat( sst.bpi( i1, : ), [ sum( i2 ) 1 ] );
%             L.fD.upol( i2, : )     = repmat( sst.upol( i1, : ), [ sum( i2 ) 1 ] );
%         end
%     end
% 
%     n_track                = size (L.tD.shankclu, 1);
%     L.tD.extremum           = NaN( n_track, 1 );
%     L.tD.bpi                = NaN( n_track, 1 );
%     L.tD.upol               = NaN( n_track, 1 );    
%     u1                      = unique( sst.filebase );
%     u2                      = unique( L.tD.session_tag );
%     u                       = unique( [ u1; u2 ] );
%     [ ~, f2 ]               = ismember( L.tD.session_tag, u );
%     s2                      = [ f2 L.tD.shankclu( :, 1 : 2 ) ];
%     us2                     = unique( s2, 'rows' );
%     nus2                    = size( us2, 1 );
%     for i                   = 1 : nus2
%         i1                  = ismember( s1, us2( i, : ), 'rows' );
%         i2                  = ismember( s2, us2( i, : ), 'rows' );
%         if sum( i1 ) > 0
%             L.tD.extremum( i2, : ) = repmat( sst.extremum( i1, : ), [ sum( i2 ) 1 ] );
%             L.tD.bpi( i2, : )      = repmat( sst.bpi( i1, : ), [ sum( i2 ) 1 ] );
%             L.tD.upol( i2, : )     = repmat( sst.upol( i1, : ), [ sum( i2 ) 1 ] );
%         end
%     end
% 
%     % combine back for Hadas
% %    s_combined              = struct_merge( s_all, s_field );
%     savename                = [ datadir_hadas fname_hadas '_extremum_bpi.mat' ];
%     save( savename, 'L' )
% %     
%     % process the data
%     eidx                    = L.uD.shankclu( :, 3 ) == 1; % PYR
%     iidx                    = L.uD.shankclu( :, 3 ) == 0; % INT
%     pidx                    = L.uD.extremum > 0; % PUNIT
%     eidx( pidx )            = 0; % NUNIT PYR 
%     iidx( pidx )            = 0; % NUNIT INT
%     colnames                = { 'bits/s', 'bits/spike' };
%     row_names               = { 'PYR', 'INT', 'Punits' };
%     info_meds               = [ median( L.uD.info( eidx, : ) )
%                                 median( L.uD.info( iidx, : ) )
%                                 median( L.uD.info( pidx, : ) ) ];
%     info_sem                = [ calc_sem( L.uD.info( eidx, : ) )
%                                 calc_sem( L.uD.info( iidx, : ) )
%                                 calc_sem( L.uD.info( pidx, : ) ) ];
                            
%     % plot the figures
%     fig5( 1 )               = figure;
%     for i = 1 : 2
%         subplot( 2, 2, i )
%         [ bh, eh ] = barwerror( 1 : 3, info_meds( :, i ), info_sem( :, i ) ...
%             , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
%         set( gca, 'XTickLabel', row_names )
%         set( gca, 'tickdir', 'out', 'box', 'off' )
%         ylabel( colnames{ i } )
%         title( 'Spatial information' )
%     end
    
    % do the same for many other interesting spatial features at the unit
    % and at the field level...
    %fig5( 1:14 ) = make_PUNIT_figures( 1:14, [], [], 'sst', sst );
    fig5( 1:17 ) = make_PUNIT_figures( 1:17, fname_hadas, datadir_hadas, 'sst', sst );
    %fig5 = make_PUNIT_figures( [ 9 11 13 ], fname_hadas, datadir_hadas, 'sst', sst );
    if ~isempty( clunum ) && ~isnan( clunum )
        for i = 1 : length( fig5 )
            figure( fig5( i ) );
            fig_title( sprintf( 'cluster %d', clunum ) ); 
        end
    end
    
    
 
    %-----
    % save the figures
    fig = fig5;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG5_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 6
% Depth
% ----------------------------------------------------------------------

if fignums( 6 )
    
    % process the data
    didx                = ~isnan( sst.depth );
    sst1                = struct_select( sst, didx );
%     pidx                = sst1.extremum > 0;
    pidx                = ispos;
    bidx                = isbip;
    pidx                = pidx(didx);
    bidx                = bidx(didx);
    f                   = sst1.geo_fwhm;
    w                   = 1 ./ sst1.fmax * 1000; % width
    t                   = sst1.tp;               % t2p
%     d                   = sst1.depth * 20;       % [um]
    amp                 = abs( sst1.extremum );
    r                   = sst1.frate;            % firing rate
    a                   = sst1.ach_com;          % burstiness
    fnames              = sst1.filebase;
    
    % compute depth in um based on probe resolution
    ufnames             = unique( sst1.filebase );
    nufnames            = length( ufnames );
    tmpcell             = strfind( ufnames, 'mC41' ); % presently supporting only one animal
    % to support multiple animals (strings), must expand the code
    tmpvec              = NaN( nufnames, 1 );
    for i               = 1 : nufnames
        if isempty( tmpcell{ i } )
            tmpvec( i ) = 0;
        else
            tmpvec( i ) = tmpcell{ i };
        end
    end
    tmpvec              = logical( tmpvec );
    hridx               = ismember( fnames, ufnames( tmpvec ) );
    d( hridx )          = sst1.depth( hridx ) * 15;         % [um]
    d( ~hridx )         = sst1.depth( ~hridx ) * 20;        % [um]
    
    % determine spatial bins:
    binsize             = 10; 
    minVal              = floor( min( d ) / binsize ) * binsize;
    maxVal              = ceil( max( d ) / binsize ) * binsize;
    rside               = ( 0 + binsize / 2 ) : binsize : ( maxVal + binsize / 2 );
    lside               = ( 0 + binsize / 2 ) : binsize : ( -minVal + binsize / 2 );
    edges               = [ fliplr( -lside ) rside ];
    binC                = ( edges( 1 : end - 1 )  + edges( 2 : end ) )' / 2;
    
    % compute histograms
    hh                  = histc( d( ~pidx & ~bidx ), edges );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h0                  = hh;
    
    hh                  = histc( d( pidx ), edges );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h1                  = hh;
    
    hh                  = histc( d( bidx ), edges );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h2                  = hh;
    
    % compute histograms according to arbitrary parsing (in layer, above, below)
    edges2              = [ min( d ) -25 25 max( d ) ];
    hh                  = histc( d( ~pidx & ~bidx ), edges2 );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h20                  = hh;
    
    hh                  = histc( d( pidx ), edges2 );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h21                  = hh;
    
    hh                  = histc( d( bidx ), edges2 );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h22                  = hh;
    
    mat                 = [ h21( : ) h20( : ) h22( : )]; % punits vs. nunits vs. bipolar
    pval                = gtest( mat, 'lr', 1, 'ind' ); % test of independence
    
    % plot the figures
    fig6( 1 )           = figure;
    
    sp1                 = subplot( 2, 2, 1 );
    sp2                 = subplot( 2, 2, 2 );
    pos1                = get( sp1, 'Position' );
    pos2                = get( sp2, 'Position' );
    pos2( 1 )           = sum( pos1( [ 1 3 ] ) );
    set( sp2, 'Position', pos2 );
    
    subplot( sp1 )
    barh( binC, -h0, 1, 'EdgeColor', colors_NP( 1, : ), 'FaceColor', colors_NP( 1, : ) )
    xlabel( 'Nunit count' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    % to add dashed lines at sub-population median:
%     xx                  = d( ~pidx );
%     sepXX               = 300;
%     med1                = median( xx( xx <= sepXX ) );
%     med2                = median( xx( xx > sepXX ) );
%     alines( med1, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
%     alines( med2, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( sp2 )
    barh( binC, h1, 1, 'EdgeColor', colors_NP( 2, : ), 'FaceColor', colors_NP( 2, : ) )
    xlabel( 'Punit count' )
    %ylabel( 'Depth [\mum]' )
    set( gca, 'YTickLabel', [] )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    % to add dashed lines at sub-population median:
%     xx                  = d( pidx );
%     sepXX               = 300;
%     med3                = median( xx( xx <= sepXX ) );
%     med4                = median( xx( xx > sepXX ) );
%     alines( med3, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
%     alines( med4, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );

    subplot( 2, 2, 3 )
    frcts               = h1 ./ ( h0 + h1 + h2 );
    minCount            = 2;
    frcts( ( h0 + h1 + h2 ) <= minCount ) = NaN;
    barh( binC, frcts, 1, 'EdgeColor', colors_NP( 2, : ), 'FaceColor', colors_NP( 2, : ) )
    xlabel( 'Fraction of Punits' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    title( sprintf( 'pval (g-test) = %0.3g', pval ) )
    alines( sum( h1 ) / sum( h0 + h1 + h2 ), 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    rez{ 6, 1 }.binC                    = binC;
    rez{ 6, 1 }.h1                      = h1;
    rez{ 6, 1 }.h0                      = h0;
 
    fig6( 2 )           = figure;
    
    sp1                 = subplot( 2, 2, 1 );
    sp2                 = subplot( 2, 2, 2 );
    pos1                = get( sp1, 'Position' );
    pos2                = get( sp2, 'Position' );
    pos2( 1 )           = sum( pos1( [ 1 3 ] ) );
    set( sp2, 'Position', pos2 );
    
    subplot( sp1 )
    barh( binC, -h0, 1, 'EdgeColor', colors_NPB( 1, : ), 'FaceColor', colors_NPB( 1, : ) )
    xlabel( 'Nunit count' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( sp2 )
    barh( binC, h2, 1, 'EdgeColor', colors_NPB( 3, : ), 'FaceColor', colors_NPB( 3, : ) )
    xlabel( 'Bipolar count' )
    %ylabel( 'Depth [\mum]' )
    set( gca, 'YTickLabel', [] )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );

    subplot( 2, 2, 3 )
    frcts               = h2 ./ ( h0 + h1 + h2 );
    minCount            = 2;
    frcts( ( h0 + h1 + h2 ) <= minCount ) = NaN;
    barh( binC, frcts, 1, 'EdgeColor', colors_NPB( 3, : ), 'FaceColor', colors_NPB( 3, : ) )
    xlabel( 'Fraction of Bipolar units' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    title( sprintf( 'pval (g-test) = %0.3g', pval ) )
    alines( sum( h2 ) / sum( h0 + h1 + h2 ), 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    % plot some statstic by depth
    fig6( 3 ) = figure;
    subplot( 2, 3, 1 )
    hold on, 
    ph( 1 ) = plot( f( ~pidx & ~bidx ), d( ~pidx & ~bidx ), '.' ); 
    ph( 2 ) = plot( f( pidx ), d( pidx ), '.' ); 
    ph( 3 ) = plot( f( bidx ), d( bidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NPB( 1, : ) )    
    set( ph( 2 ), 'color', colors_NPB( 2, : ) )
    set( ph( 3 ), 'color', colors_NPB( 3, : ) )
    xlabel( 'FWHM' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 3, 2 )
    hold on, 
    ph( 1 ) = plot( amp( ~pidx & ~bidx ), d( ~pidx & ~bidx ), '.' ); 
    ph( 2 ) = plot( amp( pidx ), d( pidx ), '.' ); 
    ph( 3 ) = plot( f( bidx ), d( bidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NPB( 1, : ) )    
    set( ph( 2 ), 'color', colors_NPB( 2, : ) )    
    set( ph( 3 ), 'color', colors_NPB( 3, : ) )
    xlabel( 'Amplitude [\muV]' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
      
    subplot( 2, 3, 3 )
    hold on, 
    ph( 1 ) = plot( w( ~pidx & ~bidx ), d( ~pidx & ~bidx ), '.' ); 
    ph( 2 ) = plot( w( pidx ), d( pidx ), '.' ); 
    ph( 3 ) = plot( w( bidx ), d( bidx ), '.' );  
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NPB( 1, : ) )    
    set( ph( 2 ), 'color', colors_NPB( 2, : ) )    
    set( ph( 3 ), 'color', colors_NPB( 3, : ) )    
    xlabel( 'Width [ms]' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 3, 4 )
    hold on, 
    ph( 1 ) = plot( t( ~pidx & ~bidx ), d( ~pidx & ~bidx ), '.' ); 
    ph( 2 ) = plot( t( pidx ), d( pidx ), '.' ); 
    ph( 3 ) = plot( f( bidx ), d( bidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NPB( 1, : ) )    
    set( ph( 2 ), 'color', colors_NPB( 2, : ) )    
    set( ph( 3 ), 'color', colors_NPB( 3, : ) )
    xlabel( 't2p [ms]' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );    
    
    subplot( 2, 3, 5 )
    hold on, 
    ph( 1 ) = plot( r( ~pidx & ~bidx ), d( ~pidx & ~bidx ), '.' ); 
    ph( 2 ) = plot( r( pidx ), d( pidx ), '.' ); 
    ph( 3 ) = plot( f( bidx ), d( bidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NPB( 1, : ) )    
    set( ph( 2 ), 'color', colors_NPB( 2, : ) )    
    set( ph( 3 ), 'color', colors_NPB( 3, : ) )
    xlabel( 'Firing rate [spk/s]' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );   
    
    subplot( 2, 3, 6 )
    hold on, 
    ph( 1 ) = plot( a( ~pidx & ~bidx ), d( ~pidx & ~bidx ), '.' ); 
    ph( 2 ) = plot( a( pidx ), d( pidx ), '.' ); 
    ph( 3 ) = plot( f( bidx ), d( bidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NPB( 1, : ) )    
    set( ph( 2 ), 'color', colors_NPB( 2, : ) )    
    set( ph( 3 ), 'color', colors_NPB( 3, : ) )
    xlabel( 'Burstiness' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );       
    
    %-----
    % save the figures
    fig = fig6;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG6_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 7
% optical responses
% ----------------------------------------------------------------------

if fignums( 7 )
    
    % process the data
    
    % go over sessions
    % for each session, check if there is a celltypeClassification 
    % if there is, load it, and check the following:
    % -opsinType field is no empty
    % -peth field has second dimension with 201 elements (default of celltypeClassification)
    
    % select, for each unit, the local stimulation channel and extract the 
    % peth, bins, pact, psup, gain
    % also, keep the act
    
    % module #1 - 
    % go over global sst stucture and get celltypeClassification results for each session

    % parameters that depend on the parameters used during celltypeClassification 
    nPethBins       = 201;
    ilevel          = 'B'; 
    
    % where to load the data from - either datadir or the actual filebase
    loadMode        = 'datadir';
    wavRange        = [ 350 500 ];                          % UV LED, blue LD, blue LED
    %wavRange        = [ 550 700 ];
    
    usess           = unique( sst.filebase );
    
    %%%%%%%%%%%%%%
    % remove this line 
    %usess           = usess( [ 76 78 ] );
    %%%%%%%%%%%%%%
    nsess           = length( usess );
    
    ii                              = 0;
    for i                           = 1 : nsess
        % the indices in the global sst
        asess                       = usess{ i };
        uidx                        = ismember( sst.filebase, asess );
        % initialize an output structure for the session
        ssti                        = struct_select( sst, uidx );
        nunits                      = size( ssti.filebase, 1 );
        si.filebase                 = repmat( ssti.filebase( 1 ), [ nunits 1 ] );
        si.shankclu                 = NaN( nunits, 3 );

        lv                          = false( nunits, 1 );
        nv                          = zeros( nunits, 1 );
        cv                          = cell( nunits, 1 );
        pv                          = zeros( nunits, nPethBins );
        
        si.opsinType                = cv;               % cell array of strings, e.g. pv::chr2
        si.act                      = lv;               % boolean flag - is light controlled (activated or silenced)
        si.exc                      = lv;               % boolean flag - is an excitatory unit
        si.inh                      = lv;               % boolean flag - is an inhibitory unit
        
        si.pact                     = nv;               % p-value for activation
        si.psup                     = nv;               % p-value for silencing
        si.gain                     = nv;               % gain
        
        si.peth                     = pv;               % PETH during local stimulation
        si.bins                     = pv;               % PETH bins during local stimulation

        % load the s and the shank/stimulation asssociation for each filebase
        try
            switch loadMode
                case 'filebase'
                    filebase                = filebaseLookup( si.filebase{ 1 } );
                    % get the association between channels, targets, and wavelengths
                    [ chans, targets, ~, ~, wavelength ]      = get_stimchans( filebase );
                    % get the results of the celltype classification 
                    sctc                    = celltypeClassification( filebase, 'Overwrite', -1 );
                case 'datadir'
                    suffix                  = si.filebase{ 1 };
                    % get the association between channels, targets, and wavelengths
                    xmlfname                = [ datadir 'xml/' suffix '.prm.xml' ];
                    if ~exist( xmlfname, 'file' )
                        fprintf( 1, '%d: Missing data (*.prm.xml file) for %s\n', i, usess{ i } )
                        continue
                    end
                    par                     = LoadXml( xmlfname );
                    [ chans, targets, ~, ~, wavelength ]      = get_stimchans( par );
                    % decide which cell type file to use
                    kidx                   = inrange( wavelength, wavRange );
                    if sum( kidx ) == 0
                        fprintf( 1, '%d: Missing data (no stim channels within %d to %d nm ) for %s\n'...
                            , i, wavRange( 1 ), wavRange( 2 ), usess{ i } )
                        continue
                    end
                    chans( ~kidx )          = [];
                    targets( ~kidx )        = [];
                    wavelength( ~kidx )     = [];
                    ctwildcard              = sprintf( '%sdc/%s.*_c%s.celltypeClassification', datadir, suffix, ilevel );
                    fnames                  = dir( ctwildcard );
                    if isempty( fnames )
                        fprintf( 1, '%d: Missing data (no celltypeClassification at ilevel %s) for %s\n'...
                            , i, ilevel, usess{ i } )
                        continue

                    end
                    ctfname                 = '';
                    for j = 1 : length( fnames )
                        dots                = strfind( fnames( j ).name, '.' );
                        corename            = fnames( j ).name( ( dots( 1 ) + 1 ) : ( dots( 2 ) - 1 ) );
                        if contains( corename, num2strs( chans ) )
                            ctfname         = fnames( j ).name;
                            break
                        end
                    end
                    if isempty( ctfname )
                        fprintf( 1, '%d: Missing data (no celltypeClassification within %d to %d nm) for %s\n'...
                            , i, wavRange( 1 ), wavRange( 2 ), usess{ i } )
                        continue
                    end
                    % get the results of the celltype classification 
                    ctfullname              = sprintf( '%sdc/%s', datadir, ctfname );
%                     L                       = load( ctfullname, '-mat' );
%                     s = L.z.s; stats = L.z.stats; save( ctfullname, 's', 'stats' )
                    L                       = load( ctfullname, '-mat' );
                    try
                        sctc                    = L.s;
                    catch
                        sctc                = L.z.s;
                    end
                    % temporary hack - convert nested cell arrays to strings
                    for l                   = 1 : size( si.opsinType, 1 )
                        otype               = si.opsinType{ l };
                        if isa( otype, 'cell' ) && isa( otype{ 1 }, 'char' )
                            si.opsinType{ l } = otype{ 1 };
                        end
                    end
            end
        catch
            fprintf( 1, '%d: Missing data (either celltype, stimchans, or both) for %s\n', i, usess{ i } )
            continue
        end
        if ~isequal( sctc.shankclu( :, 1 :  2 ), ssti.shankclu )
            fprintf( 1, '%d: Data mismatch for %s\n', i, usess{ i } )
            continue
        end
        shanks                      = sctc.shankclu( :, 1 );
        ushanks                     = unique( shanks );
        nshanks                     = length( ushanks );
        for j                       = 1 : nshanks
            ashank                  = ushanks( j );
            achan                   = chans( targets == ashank );
            if isempty( achan )
                continue
            end
            uidx                    = sctc.shankclu( :, 1 ) == ashank;
            imat                    = permute( sctc.pinfo( uidx, 6, : ), [ 1 3 2 ] ) == achan;
            if size( imat, 1 ) ~= max( sum( imat, 1 ) ) && max( sum( imat, 1 ) ) ~= 0
                fprintf( 1, '%d: Data mismatch for %s in shank %d\n', i, usess{ i }, ashank )
                continue
                % thus, all units in the shank are associated with the same trigchan
            end
            [ ~, cidx ]             = find( imat ); % vector of the index of the trigchan
            if isempty( cidx )
                continue
            end
            if length( unique( cidx ) ) ~= 1
                error( 'should never happen' )
            end
            si.shankclu( uidx, : )  = sctc.shankclu( uidx, : );
            si.opsinType( uidx, : ) = sctc.opsinType( uidx, : );
            si.act( uidx, : )       = sctc.act( uidx, : );
            si.exc( uidx, : )       = sctc.exc( uidx, : );
            si.inh( uidx, : )       = sctc.inh( uidx, : );
            si.pact( uidx, : )      = sctc.pact( uidx, cidx( 1 ) );
            si.psup( uidx, : )      = sctc.psup( uidx, cidx( 1 ) );
            si.gain( uidx, : )      = sctc.gain( uidx, cidx( 1 ) );
            si.peth( uidx, : )      = sctc.peth( uidx, :, cidx( 1 ) );
            si.bins( uidx, : )      = sctc.bins( uidx, :, cidx( 1 ) );
        end
        % keep only populated units
        kidx                        = ~isnan( si.shankclu( :, 1 ) );
        si                          = struct_select( si, kidx, 1 );
        % accumulate over sessions
        ii                          = ii + 1;
        if ii == 1
            s                       = si;
        else
            s                       = struct_cat( s, si );
        end
            
    end
    
    
    % check if equal
%     if ~isequal( s.filebase, sst.filebase  )
%         error( 'fix all issues first' )
%     end
%     for i = 1 : size( s.filebase, 1 )
%         eq( i ) = isequal( s.filebase( i ), sst.filebase( i ) ); 
%     end
%     neq = find( ~eq );
    
% module #3 - intersect units

    % intersect the units in the sst and the ff structures
    u1                      = unique( sst.filebase );
    u2                      = unique( s.filebase );
    u                       = unique( [ u1; u2 ] );
    [ ~, f1 ]               = ismember( sst.filebase, u );
    [ ~, f2 ]               = ismember( s.filebase, u );
    s1                      = [ f1 sst.shankclu ];
    s2                      = [ f2 s.shankclu( :, 1 : 2 ) ];
    % check that unique and no NaNs
    if ~isequal( size( unique( s1, 'rows' ) ), size( s1 ) )
        error( 'Repeating items in s1\n' )
    end
    if ~isequal( size( unique( s2, 'rows' ) ), size( s2 ) )
        error( 'Repeating items in s2\n' )
    end
    if sum( isnan( s1( : ) ) ) ~= 0
        error( 'NaNs in s1\n' )
    end
    if sum( isnan( s2( : ) ) ) ~= 0
        error( 'NaNs in s2\n' )
    end
    % since each item appears only once, we can use intersect
    [ ~, i1, i2 ]           = intersect( s1, s2, 'rows' );
    i1l                     = false( size( sst.shankclu, 1 ), 1 );
    i1l( i1 )               = 1;
    i2l                     = false( size( s.shankclu, 1 ), 1 );
    i2l( i2 )               = 1;
    sst1                    = struct_select( sst, i1l );
    ff1                     = struct_select( s, i2l );
    % merge the two structures
    sst1                    = struct_merge( rmfield( sst1, 'shankclu' ), rmfield( ff1, 'filebase' ) );
    % add tag for punits
%     pidx                    = sst1.extremum > 0;
    isbip           = ~isnan(sst1.bpi) & ~isinf( sst1.bpi );
    ispos           = sst1.extremum > 0 & ~isbip;
    sst1.shankclu( ispos, 3 ) = 2;  % step over the INT/PYR for punits
    sst1.shankclu( isbip, 3 ) = 3;  % step over the INT/PYR for biphasic
    
    % temporary hack - convert nested cell arrays to strings
    for i                   = 1 : size( sst1.opsinType, 1 )
        otype               = sst1.opsinType{ i };
        if isa( otype, 'cell' ) && isa( otype{ 1 }, 'char' )
            sst1.opsinType{ i } = otype{ 1 };
            
        end
    end
    
    % temporary hack #2 - remove items for which opsinType is empty
    ridx = false( length( sst1.opsinType ), 1 ); 
    for i = 1 : length( sst1.opsinType )
        if isempty( sst1.opsinType{ i } )
            ridx( i ) = 1;
        end
    end
    if sum( ridx ) > 0
        fprintf( 1, 'Removing %d units without opsinType\n', sum( ridx ) )
        sst1 = struct_select( sst1, ~ridx );
    end
    
    % temporary hack #3 - change camk::chr2 to camkii::chr2
    hidx                    = ismember( sst1.opsinType, 'camk::chr2' ); 
    sst1.opsinType( hidx )  = repmat( { 'camkii::chr2' }, [ sum( hidx ) 1 ] );
    
    % plot the figures
    % constants
    xlims                   = [ -100 150 ]; % [ms]
    ctypes                  = { 'INT', 'PYR', 'Punit', 'Biphasic' };
    opsinTypes              = { 'camkii::chr2', 'vip::chr2', 'pv::chr2', 'pv::jaws', 'cck::chr2' };
    
    % select a subset of units
    aidx                    = sst1.act;
    j                       = 0;
    for onum                = 1 : length( opsinTypes )
        opsinType           = opsinTypes{ onum };
        oidx                = ismember( sst1.opsinType, opsinType );
        if sum( oidx ) == 0
            fprintf( 1, 'No units with %s found\n', opsinType )
        end
        if sum( oidx & aidx ) == 0
            fprintf( 1, 'No units with %s and optical activation/silencing found\n', opsinType )
        end
        
        j                   = j + 1;
        [ ah, fh ]          = tilefig( 1, 4, -2, 0.75, 'center' );
        fig7( j )           = fh;
        xvals                   = median( sst1.bins( oidx, : ), 1 )' * 1000; % time [ms]
        for ct                  = 0 : 3
            cidx                = sst1.shankclu( :, 3 ) == ct;
            idx                 = aidx & oidx & cidx;
            % to do: 
            % keep the following and create bar with error bars
            %sum( idx )
            %sum( oidx & cidx )
            if sum( idx ) == 0
                subplot( ah( ct + 1, 1 ) )
                axis off
                subplot( ah( ct + 1, 2 ) )
                axis off
                continue
            end
            peth                = sst1.peth( idx, : )';
            yvals               = 1 : size( peth, 2 );
            speth               = scale( peth );
            myu                 = mean( peth, 2 );
            sem                 = calc_sem( peth, 2 );
            %myu                 = mean( speth, 2 );
            %sem                 = calc_sem( speth, 2 );
            switch ct
                case { 0, 1 }
                    color1      = colors_PI( ct + 1, : );
                    color2      = colors_PI_light( ct + 1, : );
                case 2
                    color1      = colors_NP( 2, : );
                    color2      = colors_NP_light( 2, : );
                case 3
                    color1      = colors_NPB( 2, : );
                    color2      = colors_NP_light( 2, : );
            end
            
            subplot( ah( ct + 1, 1 ) )
            patch_band( xvals, myu, sem, color1, color2 );
            xlim( xlims )
            set( gca, 'tickdir', 'out', 'box', 'off' );
            ylabel( 'Firing rate [spks/s]' )
            title( sprintf( '%s, %s', opsinType, ctypes{ ct + 1 } ) )
            
            subplot( ah( ct + 1, 2 ) )
            imagesc( xvals, yvals, speth' )
            axis xy,
            colormap( myjet )
            xlim( xlims )
            xlabel( 'Time [ms]' )
            ylabel( 'Unit number' )
            set( gca, 'ytick', yvals )
            %yvals( 1 : ceil( length( yvals ) / 10 ) : length( yvals ) )
            set( gca, 'tickdir', 'out', 'box', 'off' );
            
            % add stim lines and patch
            
        end
    end
    %s_light=find_pv_punits (sst1,pidx);
    %fig7( 2 ) = figure;

    %-----
    % save the figures
    fig = fig7;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG7_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 8
% subclusters
% ----------------------------------------------------------------------

trif fignums( 8 )
    
    % process the data
    ispos                       = sst.extremum > 0;

    sst.shankclu( :, 3 )        = sst.pyr;
    sst.shankclu( ispos, 3 )    = 2;
    
    % plot the figures
    
    %------------------------------------------------------------------
    % (1) histograms of each feature - nunits vs. punits (similar to figure 2 - "by prob")

    nbins1          = 60;
    nbins           = 30;
    byprob          = 1;
    make_legend     = 0;
    
    fig8( 1 ) = figure;
    subplot( 2, 3, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.maxp2p * 1000, ispos, nbins1, colors_NP, 'Amp [\muV]', byprob );
    
    subplot( 2, 3, 2 );
    [ t2p_myus, t2p_sds, t2p_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.tp, ispos, nbins, colors_NP, 'T2P [ms]', byprob, make_legend );

    subplot( 2, 3, 3 );
    [ t2p_myus, t2p_sds, t2p_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.tp, ispos, nbins, colors_NP, 'T2P [ms]', byprob, make_legend );
    ylim( [ 0 0.2 ] )
    
    subplot( 2, 3, 4 );
    [ wid_myus, wid_sds, wid_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( 1 ./ sst.fmax * 1000, ispos, nbins, colors_NP, 'Width [ms]', byprob, make_legend );
    
    subplot( 2, 3, 5 );
    [ gsd_myus, gsd_sds, gsd_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.geo_sd, ispos, nbins, colors_NP, 'Geo-SD [sites]', byprob, make_legend );
    
    subplot( 2, 3, 6 );
    [ fwh_myus, fwh_sds, fwh_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.geo_fwhm, ispos, nbins, colors_NP, 'Geo-FWHM [sites]', byprob, make_legend );

    for i               = 1 : 6
        subplot( 2, 3, i )
        axis square
    end
    
    fig8( 2 ) = figure;

    subplot( 2, 3, 1 );
    [ asy_myus, asy_sds, asy_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.asy, ispos, nbins, colors_NP, 'Assymetry', byprob, make_legend );
    
    subplot( 2, 3, 2 );
    [ myus, sds, pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.bpi, ispos, nbins, colors_NP, 'Bi-Polarity', byprob, make_legend );
    alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 3, 3 );
    [ myus, sds, pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.bpi, ispos, nbins, colors_NP, 'Bi-Polarity', byprob, make_legend );
    ylim( [ 0 0.1 ] )
    alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 3, 4 );
    [ myus, sds, pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.tV, ispos, nbins, colors_NP, 'SD of peak time [ms]', byprob, make_legend );

    subplot( 2, 3, 5 );
    [ bst_myus, bst_sds, bst_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.brst, ispos, nbins, colors_NP, 'Burst', byprob, make_legend );
    alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 3, 6 );
    [ acm_myus, acm_sds, acm_pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.ach_com, ispos, nbins, colors_NP, 'ACH-COM [ms]', byprob, make_legend );
    
    for i               = 1 : 6
        subplot( 2, 3, i )
        axis square
    end
    
    %------------------------------------------------------------------
    % (2) scatter of t2p and width - PYR, INT, punits (scatter with jitter)
    cti                     = sst.shankclu( :, 3 );
    t2p                     = sst.tp;
    wid                     = 1 ./ sst.fmax * 1000;
    fwh                     = sst.geo_fwhm;
    amp                     = sst.maxp2p * 1000;
    acm                     = sst.ach_com;

    gsd                     = sst.geo_sd;
    asy                     = sst.asy;
    bst                     = sst.brst;
    bpi                     = sst.bpi;
    tv                      = sst.tV;
    x                       = [ t2p fwh amp acm bpi wid gsd asy bst tv ];
    
    % prepare the separatrix
    [ ~, ~, fsep ]          = classify_waveform( [ t2p( cti ~= 2 ) wid( cti ~= 2 ) ] );
%     xx                  = [ 0.1 0.9 ];
%     yy                  = [ 0.6 1.4 ];
    xx                      = [ 0.025 0.95 ];
    yy                      = [ 0.35 1.4 ];
%     xx                  = [ 0.03 1.13 ];
%     yy                  = [ 0.25 1.55 ];
    MSs                     = [ 8 8 12 ];
    jitflag                 = 1;
    
    fig8( 3 )                   = figure;
    fig8( 4 )                   = figure;
    for spi = 1 : 4
        switch spi
            case 1
                yval            = wid;
                param2_ylabel   = 'Width [ms]';
                param2_ylim     = yy;
            case 2
                yval            = fwh;
                param2_ylabel   = 'FWHM [sites]';
                param2_ylim     = [ 0.2 9.8 ];
            case 3
                yval            = amp;
                param2_ylabel   = 'AMP [\muV]';
                param2_ylim     = [ 20 1200 ];
            case 4
                yval            = acm;
                param2_ylabel   = 'ACH COM [ms]';
                param2_ylim     = [ 5 40 ];
        end
        
        % jitter points to improve visualization
        if jitflag
            nclu                    = length( t2p );
            dd                      = sort( diff( sort( t2p ) ) );
            dd                      = dd( dd > 0 );
            tt                      = sort( diff( sort( yval ) ) );
            tt                      = tt( tt > 0 );
            if isempty( dd )
                jx                  = 0;
            else
                dd                  = dd( ceil(  length( dd ) * 0.05 ) );
                jx                  = dd / 2 * randn( nclu, 1 );
            end
            if isempty( tt )
                jy                  = 0;
            else
                tt                  = tt( ceil(  length( tt ) * 0.05 ) );
                jy                  = tt * randn( nclu, 1 );
            end
        else
            jx                      = 0;
            jy                      = 0;
        end
        
        if spi == 1
            figure( fig8( 3 ) )
            fh                          = fimplicit( fsep, [ xx yy ] );
            set( fh, 'color', [ 0 0 0 ] );
        else
            figure( fig8( 4 ) )
            subplot( 2, 2, spi - 1 )
        end
        
        hold on
        for ct                      = 0 : 2
            cidx                    = cti == ct;
            ph                      = plot( t2p( cidx ) + jx( cidx ), yval( cidx ) + jy( cidx ), '.' );
            switch ct
                case 0
                    set( ph, 'color', colors_PI( 1, : ), 'MarkerSize', MSs( 1 ) )
                case 1
                    set( ph, 'color', colors_PI( 2, : ), 'MarkerSize', MSs( 2 ) )
                case 2
                    set( ph, 'color', colors_NP_light( 2, : ), 'MarkerSize', MSs( 3 ) )
            end
        end
        axis square
        xlim( xx )
        ylim( param2_ylim )
        set( gca, 'tickdir', 'out', 'box', 'off' );
        axis square
        xlabel( 'T2P [ms]' )
        ylabel( param2_ylabel )
    end

    % (2.1) scatter of other pairs of features
    
    %------------------------------------------------------------------
    % (3) cluster punits

    clustMode               = 'sic';
    feature_names           = { 'T2P', 'FWHM', 'AMP', 'ACM_COM', 'BPI', 'WID', 'GSD', 'ASY', 'BURST', 'SD_tV' };
    
    % cluster all punits
    didx                    = ( 1 : 5 );
    xhat                    = x( ispos, didx );
   
    % making bipolar index binary
    idx_hat = xhat(:,5) > -0.7;
    xhat(idx_hat,5) = 1;
    xhat(~idx_hat,5) = -1;
    
    %clustRunMode                = 'saveRun';
    
    switch clustRunMode
        
        case 'localRun'
            
            fprintf( 'Running optclust for %d samples, %d features...', size( xhat, 1 ), size( xhat, 2 ) )
            clu                     = optclust( xhat, clustMode );
            nclu                    = length( unique( clu ) );
            fprintf( 'found %d clusters!\n', nclu )
            
            clear L
            L.x                     = x; 
            L.feature_names         = feature_names; 
            L.xhat                  = xhat; 
            L.clu                   = clu;
            
        case 'saveRun'
            
            % after dataset is finalized, should
            % (1) run clustering many (e.g. 20) times to determine number of clusters (nclus):
            nreps1                  = 20;
            nclus                   = NaN( nreps1, 1 );
            for i                   = 1 : nreps1
                fprintf( 1, 'repetition #%d: clustering with all possible number of clusters...\n', i )
                rng( round( rand( 1 ) * sum( 100 * clock ) ), 'v4' ) % randomly initialize the RNG used by randi
                clu                 = optclust( xhat, clustMode );
                nclus( i )          = length( unique( clu ) );
            end
            
            % (2) for each possible number of clusters, run clustering many times (e.g. 20) to determine stability of this number:
            npunits                 = sum( ispos );
            nreps2                  = 20;
            
            ncluAll                 = unique( nclus );
            nnclus                  = length( ncluAll );
            ncluActual              = NaN( nnclus, nreps2 );
            for i                   = 1 : nnclus
                nclu                = ncluAll( i );
                clutags                 = NaN( npunits, nreps2 );
                fprintf( 1, 'starting with %d clusters', nclu )
                for j                   = 1 : nreps2
                    fprintf( 1, '.' )
                    rng( round( rand( 1 ) * sum( 100 * clock ) ), 'v4' ) % randomly initialize the RNG used by randi
                    clu                 = optclust( xhat, clustMode, [], [ nclu nclu ] );
                    clutags( :, j )     = clu;
                    ncluActual( i, j )  = length( unique( clutags( :, j ) ) );
                end
                fprintf( 1, '\n' )
            end
            % generate a distance metric that quantifies the fraction of cases in
            % which the actual number of clusters is the same as the requested ("stability")
            dm                      = mean( ncluAll * ones( 1, nreps2 ) - ncluActual, 2 ) ./ ncluAll; % weighted
            %dm                      = mean( ( ncluAll * ones( 1, nreps2 ) - ncluActual ) ~= 0, 2 ); % non-weighted
            minval                  = min( dm );
            ncluSelected            = ncluAll( dm == minval );
            if length( ncluSelected ) > 1
                fprintf( 'More than one number of clusters with same score\n' )
            end
            nclu                    = ncluSelected( 1 ); % nclu=3
            
            % (3) fix the number of clusters, and run clustering many times (e.g. 20) to determine tagging of each individual unit:
            nreps3                      = 20;
            clutags                     = NaN( npunits, nreps3 );
            ncluFinal                   = NaN( 1, nreps3 );
            fprintf( 1, 'Decided on %d clusters', nclu )
            for j                       = 1 : nreps3
                fprintf( 1, '.' )
                rng( round( rand( 1 ) * sum( 100 * clock ) ), 'v4' ) % randomly initialize the RNG used by randi
                % optclust -> gmdistribution -> gmclust -> randsample -> randi
                clu                     = optclust( xhat, clustMode, [], [ nclu nclu ] );
                clutags( :, j )         = clu;
                ncluFinal( j )          = length( unique( clutags( :, j ) ) );
            end
            nidx                        = ncluFinal ~= nclu;
            if sum( nidx )
                fprintf( 1, 'Removed %d solutions\n', sum( nidx ) )
            end
            clutags( :, nidx )          = [];
            nreps                       = sum( ~nidx );
            [ ii, jj ]                  = find( clutags - round( nanmean( clutags, 2 ) ) * ones( 1, nreps ) );
            nUnstable                   = length( ii );
            if nUnstable > 0
                fprintf( 1, '%d/%d unstable solutions (%d data points)\n', nUnstable, numel( clutags ), npunits )
            else
                fprintf( 1, 'All points are stable (in terms of clustering and naming)\n' )
            end
            clu                         = round( nanmean( clutags, 2 ) );
            
            % (4) save this to file and use the fixed clustering results
            save( [ datadir 'punits_clusters' ], 'x', 'feature_names', 'xhat', 'clu' )
            clear L
            L.x                     = x; 
            L.feature_names         = feature_names; 
            L.xhat                  = xhat; 
            L.clu                   = clu;
            
        case 'loadRun'
            
            L                       = load( [ datadir 'punits_clusters.mat' ] );
    end

    clear x feature_names xhat clu
    nclu                    = length( unique( L.clu ) );
    
    % assume nclu of three clusters, and assign colors and names
    clunames                = cell( size( L.clu ) ); 
    uclunames               = { 'SPU (INT)', 'LPU (axons/dendrites)', 'JUX (PYR)' };
    if nclu>3
        for i               =  1 : nclu
        uclunames{i} = sprintf('%d', i);
        end
    end
    for i                   = 1 : nclu
        clunames( L.clu == i )  = uclunames( i );
    end
    colors                  = [ colors_PI( 1, : ); colors_NP( 2, : ); colors_PI( 2, : ) ; [1 0 0]; [0 0 1]; [0 0 0]];
   
    % scatter of t2p and FWHM - punits only, data-dependent clustering (in high-dimensional space)
    
    %sstN                    = struct_select( sst, ~ispos );
    %labels                  = sstN.pyr;
    %param1                  = sstN.tp;
    
    labels                  = sst.pyr( ~ispos );

    SFs                     = [ 1 1/10 1 1/50 ];
    MS_Punits               = 10;
    
    fig8( 5 )               = figure;
    %colors                      = 'rgbmkcyrgbmkcy';
    pairs                   = make_ai( size( L.xhat, 2 ) );
    for spi                 = 1 : 6

        f1n                 = pairs( spi, 1 );
        f2n                 = pairs( spi, 2 );
        param1_label        = L.feature_names{ f1n };
        param2_label        = L.feature_names{ f2n };
        param1n             = L.x( ~ispos, f1n );
        param2n             = L.x( ~ispos, f2n );
        param1p             = L.xhat( :, f1n );
        param2p             = L.xhat( :, f2n );
        param1              = L.x( :, f1n );
        param2              = L.x( :, f2n );
        m1                  = max( abs( param1 ) );
        m2                  = max( abs( param2 ) );
        
        subplot( 2, 3, spi )
        hold on,
        for i                   = 1 : nclu
            %ph                  = plot( param1p( L.clu == i ) / m1, param2p( L.clu == i ) * SFs( f2n ) / m2, '.' );
            ph                  = plot( param1p( L.clu == i ) / m1, param2p( L.clu == i ) / m2, '.' );
            set( ph, 'color', colors( i, : ), 'markersize', MS_Punits )
        end
        title( sprintf( 'Pos dataset (%d clusters)', nclu ) )
        xlabel( param1_label )
        ylabel( param2_label )
        axis square
        %make_punits_figures_add_iso_ellipses( param1n / m1, param2n * SFs( f2n ) / m2, labels, colors_PI_light );            % add ellipses
        make_punits_figures_add_iso_ellipses( param1n / m1, param2n / m2, labels, colors_PI_light );            % add ellipses
        % scale back the ticks of the y-axis
        lbls                    = get( gca, 'YTickLabel' );
        yticks                  = zeros( 1, length( lbls ) );
        for i                   = 1 : length( lbls )
            %yticks( i )         = str2double( lbls{ i } ) / SFs( f2n );
            try
                yticks( i )         = str2double( lbls{ i } ) * m2;
            catch
                yticks( i )         = str2double( lbls( i ) ) * m2;
            end
        end
        ytick                   = get( gca, 'ytick' );
        set( gca, 'ytick', ytick )
        set( gca, 'YTickLabel', round( yticks * 10 ) / 10 );
        
        % scale back the ticks of the x-axis
        lbls                    = get( gca, 'XTickLabel' );
        xticks                  = zeros( 1, length( lbls ) );
        for i                   = 1 : length( lbls )
            %yticks( i )         = str2double( lbls{ i } ) / SFs( f2n );
            try
                xticks( i )         = str2double( lbls{ i } ) * m1;
            catch
                xticks( i )         = str2double( lbls( i ) ) * m1;
            end
        end
        xtick                   = get( gca, 'xtick' );
        set( gca, 'xtick', xtick )
        set( gca, 'XTickLabel', round( xticks * 10 ) / 10 );
        set( gca, 'tickdir', 'out', 'box', 'off' );
    end
    
    % organize x/y limits
    
    % (3.1) tSNE projection of the high-D data
    fig8( 6 )               = figure;
    Y                       = tsne( L.x( ispos, didx ), 'Algorithm', 'exact', 'Standardize', true, 'Perplexity', 20 ); 
    gh                      = gscatter( Y( :, 1 ), Y( :, 2 ), clunames );
    for i                   = 1 : nclu
        acolor              = colors( ismember( uclunames, get( gh( i ), 'DisplayName' ) ), : );
        set( gh( i ), 'color', acolor, 'markersize', MS_Punits );
    end
    axis off
%     axis square
%     set( gca, 'tickdir', 'out', 'box', 'off' );
    
    fig8( 7 ) = figure;
    for i                   = 1: nclu
        sums (i)            = sum (ismember(L.clu, i));   
    end
    ph                      = pie( sums, uclunames );   
    for i                   = 1 : nclu
%         set( ph( 2 * i - 1 ), 'FaceColor', colors(i), 'EdgeColor', colors( i ) )
        set( ph( 2 * i ), 'String', sprintf( '%s (%d)', uclunames{ i }, sums( i ) ) )
    end

    % demonstrate all features
    % including the calc_spatial_waveform_features for a couple of
    % interesting units (nunits/punits/bunits)
    

    %-----
    % save the figures
    fig = fig8;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG8_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
% end

% ----------------------------------------------------------------------
% Figure 9
% 
% ----------------------------------------------------------------------

if fignums( 9 )
    
    % process the data
    
    % plot the figures
    fig9( 1 ) = figure;
    
    fig9( 2 ) = figure;

    %-----
    % save the figures
    fig = fig9;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG9_part%d%s', outdir, i, clunumstr );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

return


%----------------------------------------------------
% internal functions
%----------------------------------------------------



%----------------------------------------------------

% EOF

% run over all clusters
rez = cell( 4, 1 );
for clunum = 1 : 3
    rez{ clunum } = make_punits_figures( 3, 'clustRunMode', [], 'clunum', clunum, 'savef', 0 );
end

make_punits_figures( 3, 'metaplot', 1, 'rez', rez );


