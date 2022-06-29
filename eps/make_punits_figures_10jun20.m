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

% next:
% update figure 3, part 8, with ASG and STC statistics

% figure 1: 
% figure 2: ubiquity
% figure 3: connectivity
% figure 4: firing rate statistics (fano factor)
% figure 5: place analyses
% figure 6: depth analyses
% figure 7:

function make_punits_figures( fignums, savef, savetype, outdir )

% general constants
Nfigs                   = 10;
colors_PI               = [ 0 0 0.7; 1 0 0 ]; % PYR, INT
colors_NP               = [ 0 0.7 0; 1 0 1 ]; % Nunit, Punit

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

if nargs < 2 || isempty( savef )
    savef = 1;
end
if nargs < 3 || isempty( savetype )
    savetype = 'pdf';
end
if nargs < 4 || isempty( outdir )
    outdir = '';
end

% get some inline functions
[ bino_ci_exact, bino_ci_norm, bino_se_norm ] = binomial_inlines;

% set plotting parameters:
if isempty( outdir ) || ~exist( outdir, 'dir' )
    if ismac        % eran (hoor)
        outdir          = '/Users/eranstark/Documents/graphics/punits/data_figures';
        datadir         = '/Users/eranstark/Documents/da/punits/';
        % assumes that data are synchronized with odin, if not, write:
        % rsync -avh --progress /Volumes/odin/Shirly/eps/fanofactor/* ~/Documents/da/punits/fanofactor/
        % rsync -avh --progress /Volumes/odin/Shirly/eps/sps/* ~/Documents/da/punits/sps/ 
        % rsync -avh --progress /Volumes/odin/Shirly/eps/sst/* ~/Documents/da/punits/sst/
        % rsync -avh --progress /Volumes/odin/Shirly/eps/s2s/*  ~/Documents/da/punits/s2s/ 
        datadir_hadas       = '/Users/eranstark/Documents/da/hifi/';
        fname_hadas         = 'all_animals_data_21May20';
    elseif isunix   % shirly (nanna) Linux
    elseif ispc     % shirly (nanna) Windows
        outdir              = 'D:\_Shirly\Lab\eps\data_figures';
        datadir             = 'G:\mice\EPS\';
    end
end

% get the database
[ res, sst, ff ]            = shirly_eps_analysis( [], 'onlygather', 1, 'ff_suffix', '_fanofactors_SWS' );

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
            figname = sprintf( '%s/punits_FIG1_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 2
% Ubiquity
% ----------------------------------------------------------------------

if fignums( 2 )
    
    % process the data
    
    % plot the figures
    
    % histograms of signed amplitudes:
    fig2( 1 ) = figure;
    
    % copied from shirly_eps_analysis:
    ispos           = sst.extremum > 0;
    %nunits          = size( ispos, 1 );
    %nbins           = ceil( nunits / 5 );

    nbins           = 60;
    
    byprob          = 0;
    [ myus, sds, pval, bins_x, hx0, hx1 ] = make_punits_figures_one_hist( sst.maxp2p .* ( 2 * ispos - 1 ) * 1000, ispos, nbins, colors_NP, 'Signed Amp [\muV]', byprob );

    % pie charts
    
    % histograms/scatter of number of punits/nunits per session
    usess           = unique( sst.filebase );
    nsess           = length( usess );
    nums            = NaN( nsess, 2 );
    for i           = 1 : nsess
        asess       = usess{ i };
        idx         = ismember( sst.filebase, asess );
        nums( i, : ) = [ sum( idx & ~ispos ) sum( idx & ispos ) ];
    end
    frcts           = nums ./ [ sum( nums, 2 ) * ones( 1, 2 ) ];

    fig2( 2 ) = figure;
    subplot( 2, 2, 1 )
    bh = bar( 1 : nsess, sort( nums( :, 2 ) ), 1 );
    set( bh, 'FaceColor', colors_NP( 1, : ), 'EdgeColor', colors_NP( 1, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number of Punits' )
    xlabel( 'Session' )
    
    subplot( 2, 2, 2 )
    bh = bar( 1 : nsess, sort( frcts( :, 2 ) ), 1 );
    set( bh, 'FaceColor', colors_NP( 1, : ), 'EdgeColor', colors_NP( 1, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of Punits' )
    xlabel( 'Session' )
    
    subplot( 2, 2, 3 )
    jit             = ( rand( nsess, 2 ) - 0.5 ) / 2;
    nums_jit        = nums + jit;
    nums_jit( nums == 0 ) = 0;
    plot( nums_jit( :, 1 ), nums_jit( :, 2 ), '.b' )
    [ cc, pp ] = calc_spearman( nums( :, 1 ), nums( :, 2 ), 1000 );
    title( sprintf( 'Number of units/session (CC=%0.2g, pval=%0.3g)', cc, pp ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlabel( 'Nunits' )
    ylabel( 'Punits' )
    
    subplot( 2, 2, 4 )
    sums = [ sum( nums( :, 1 ) > 0 & nums( :, 2 ) == 0 )
        sum( nums( :, 1 ) == 0 & nums( :, 2 ) > 0 )
        sum( nums( :, 1 ) > 0 & nums( :, 2 ) > 0 ) ];
    slice_names = { 'Only Nunits', 'Only Punits', 'Both' };
    ph = pie( sums, slice_names );
    set( ph( 1 ), 'EdgeColor', colors_NP( 1, : ), 'FaceColor', colors_NP( 1, : ) )
    set( ph( 3 ), 'EdgeColor', [ 0 0 0 ], 'FaceColor', [ 0 0 0 ] )
    set( ph( 5 ), 'EdgeColor', colors_NP( 2, : ), 'FaceColor', colors_NP( 2, : ) )
    
    %-----
    % save the figures
    fig = fig2;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG2_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 3
% Connectivity
% ----------------------------------------------------------------------

if fignums( 3 )
    
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
        si.exc          = lv;               % boolean flag - is an excitatory unit
        si.inh          = lv;               % boolean flag - is an inhibitory unit
        si.nexcPost     = nv;               % number of excited post-synaptic peers
        si.ninhPost     = nv;               % number of inhibited post-synaptic peers
        si.nexcPre      = nv;               % number of exciting pre-synaptic peers 
        si.ninhPre      = nv;               % number of inhibiting pre-synaptic peers 
        si.asg1         = nv;               % mean ASG to excited post-synaptic peers
        si.asg2         = nv;               % mean ASG to inhibited post-synaptic peers

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
    pidx        = sst.extremum > 0;
    eidx        = s.exc == 1;                   % excitatory
    iidx        = s.inh == 1;                   % inhibitory
    eidx2       = s.nexcPre ~= 0;               % is post-synaptic to an excitatory (excited)
    iidx2       = s.ninhPre ~= 0;               % inhibited
    cidx        = eidx | iidx | eidx2 | iidx2;  % connected to the simultaneously recorded units
    
    % numbers/fractions of E-cells (punits/nunits)
    emat        = NaN( 4, 2 );
    emat( 1 : 2, 1 ) = [ sum( pidx & eidx ), sum( pidx ) ]';
    emat( 1 : 2, 2 ) = [ sum( ~pidx & eidx ), sum( ~pidx ) ]';
    emat( 3, : ) = emat( 1, : ) ./ emat( 2, : );
    emat( 4, : ) = bino_se_norm( emat( 1, : ), emat( 2, : ) )';
    % numbers/fractions of I-cells (punits/nunits)
    imat        = NaN( 4, 2 );
    imat( 1 : 2, 1 ) = [ sum( pidx & iidx ), sum( pidx ) ]';
    imat( 1 : 2, 2 ) = [ sum( ~pidx & iidx ), sum( ~pidx ) ]';
    imat( 3, : ) = imat( 1, : ) ./ imat( 2, : );
    imat( 4, : ) = bino_se_norm( imat( 1, : ), imat( 2, : ) )';
    % number/fractions of null-cells (punits/nunits)
    nmat        = NaN( 4, 2 );
    nmat( 1 : 2, 1 ) = [ sum( pidx & ~eidx & ~iidx ), sum( pidx ) ]';
    nmat( 1 : 2, 2 ) = [ sum( ~pidx & ~eidx & ~iidx ), sum( ~pidx ) ]';
    nmat( 3, : ) = nmat( 1, : ) ./ nmat( 2, : );
    nmat( 4, : ) = bino_se_norm( nmat( 1, : ), nmat( 2, : ) )';

    % number/fractions of connected-cells (punits/nunits)
    cmat        = NaN( 4, 2 );
    cmat( 1 : 2, 1 ) = [ sum( pidx & cidx ), sum( pidx ) ]';
    cmat( 1 : 2, 2 ) = [ sum( ~pidx & cidx ), sum( ~pidx ) ]';
    cmat( 3, : ) = cmat( 1, : ) ./ cmat( 2, : );
    cmat( 4, : ) = bino_se_norm( cmat( 1, : ), cmat( 2, : ) )';
    
    % numbers/fractions of excited units (punits/nunits)
    emat2       = NaN( 4, 2 );
    emat2( 1 : 2, 1 ) = [ sum( pidx & eidx2 ), sum( pidx ) ]';
    emat2( 1 : 2, 2 ) = [ sum( ~pidx & eidx2 ), sum( ~pidx ) ]';
    emat2( 3, : ) = emat2( 1, : ) ./ emat2( 2, : );
    emat2( 4, : ) = bino_se_norm( emat2( 1, : ), emat2( 2, : ) )';
    
    % numbers/fractions of inhibited units (punits/nunits)
    imat2       = NaN( 4, 2 );
    imat2( 1 : 2, 1 ) = [ sum( pidx & iidx2 ), sum( pidx ) ]';
    imat2( 1 : 2, 2 ) = [ sum( ~pidx & iidx2 ), sum( ~pidx ) ]';
    imat2( 3, : ) = imat2( 1, : ) ./ imat2( 2, : );
    imat2( 4, : ) = bino_se_norm( imat2( 1, : ), imat2( 2, : ) )';
    
    % number of connected units (excitatory/inhibitory/excited/inhibited)
    
    % number of post-synaptic peers
    npost           = NaN( 3, 4 );
    npost( :, 1 )   = [ mean( s.nexcPost( eidx & ~pidx ) ) calc_sem( s.nexcPost( eidx & ~pidx ) ) sum( eidx & ~pidx ) ]';
    npost( :, 2 )   = [ mean( s.nexcPost( eidx & pidx ) ) calc_sem( s.nexcPost( eidx & pidx ) ) sum( eidx & pidx ) ]';
    npost( :, 3 )   = [ mean( s.ninhPost( iidx & ~pidx ) ) calc_sem( s.ninhPost( iidx & ~pidx ) ) sum( iidx & ~pidx ) ]';
    npost( :, 4 )   = [ mean( s.ninhPost( iidx & pidx ) ) calc_sem( s.ninhPost( iidx & pidx ) ) sum( iidx & pidx ) ]';
    % E-nunits, E-punits, I-nunits, I-punits
    pval_NPOST      = [ utest( s.nexcPost( eidx & ~pidx ), s.nexcPost( eidx & pidx ) )
        utest( s.ninhPost( iidx & ~pidx ), s.ninhPost( iidx & pidx ) ) ];
    
    % number of pre-synaptic peers
    npre            = NaN( 3, 4 );
    npre( :, 1 )    = [ mean( s.nexcPre( eidx2 & ~pidx ) ) calc_sem( s.nexcPost( eidx2 & ~pidx ) ) sum( eidx2 & ~pidx ) ]';
    npre( :, 2 )    = [ mean( s.nexcPre( eidx2 & pidx ) ) calc_sem( s.nexcPost( eidx2 & pidx ) ) sum( eidx2 & pidx ) ]';
    npre( :, 3 )    = [ mean( s.ninhPre( iidx2 & ~pidx ) ) calc_sem( s.ninhPre( iidx2 & ~pidx ) ) sum( iidx2 & ~pidx ) ]';
    npre( :, 4 )    = [ mean( s.ninhPre( iidx2 & pidx ) ) calc_sem( s.ninhPre( iidx2 & pidx ) ) sum( iidx2 & pidx ) ]';
    % nunits-E, punits-E, nunits-I, punits-I
    H_JB = [ jbtest( s.nexcPre( eidx2 & ~pidx ) ), jbtest( s.nexcPre( eidx2 & pidx ) ) ...
        jbtest( s.ninhPre( iidx2 & ~pidx ) ), jbtest( s.ninhPre( iidx2 & pidx ) ) ];
    if all( H_JB == 0 )
        [ ~, pp1 ]  = ttest2( s.nexcPre( eidx2 & ~pidx ) , s.nexcPre( eidx2 & pidx ) );
        [ ~, pp2 ]  = ttest2( s.ninhPre( iidx2 & ~pidx ), s.ninhPre( iidx2 & pidx ) );
    else
        pp1         = utest( s.nexcPre( eidx2 & ~pidx ), s.nexcPre( eidx2 & pidx ) );
        pp2         = utest( s.ninhPre( iidx2 & ~pidx ), s.ninhPre( iidx2 & pidx ) );
    end
    pval_NPRE   = [ pp1 pp2 ];
            
    % ASG
    asg             = NaN( 3, 4 );
    asg( :, 1 )    = [ nangeomean( s.asg1( eidx & ~pidx ) ) calc_sem( s.asg1( eidx & ~pidx ) ) sum( eidx & ~pidx ) ]';
    asg( :, 2 )    = [ nangeomean( s.asg1( eidx & pidx ) ) calc_sem( s.asg1( eidx & pidx ) ) sum( eidx & pidx ) ]';
    asg( :, 3 )    = [ -nangeomean( -s.asg2( iidx & ~pidx ) ) calc_sem( s.asg2( iidx & ~pidx ) ) sum( iidx & ~pidx ) ]';
    asg( :, 4 )    = [ -nangeomean( -s.asg2( iidx & pidx ) ) calc_sem( s.asg2( iidx & pidx ) ) sum( iidx & pidx ) ]';
    % E-nunits, E-punits, I-nunits, I-punits
    pval_ASG        = [ utest( s.asg1( eidx & ~pidx ), s.asg1( eidx & pidx ) )
        utest( s.asg2( iidx & ~pidx ), s.asg2( iidx & pidx ) ) ];

    fig3( 8 ) = figure;
    subplot( 3, 3, 1 )
    [ bh, eh ] = barwerror( 1 : 4, [ emat( 3, [ 2 1 ] ) imat( 3, [ 2 1 ] ) ], [ emat( 4, [ 2 1 ] ) imat( 4, [ 2 1 ] )] ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'E-nunits', 'E-punits', 'I-nunits', 'I-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Units that are excitatory (E) or inhibitory (I)' )
    
    subplot( 3, 3, 2 )
    [ bh, eh ] = barwerror( 1 : 4, [ emat2( 3, [ 2 1 ] ) imat2( 3, [ 2 1 ] ) ] ...
        , [ emat2( 4, [ 2 1 ] ) imat2( 4, [ 2 1 ] ) ] ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Units that receive excitation or inhibition' )
    
    subplot( 3, 3, 3 )
    [ bh, eh ] = barwerror( 1 : 2, cmat( 3, [ 2 1 ] ), cmat( 4, [ 2 1 ] ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Connected-nunits', 'Connected-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Connected units' )
    
    subplot( 3, 3, 4 )
    [ bh, eh ] = barwerror( 1 : 4, npost( 1, : ), npost( 2, : ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excitatory-nunits', 'Excitatory-punits', 'Inhibitory-nunits', 'Inhibitory-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number of units' )
    title( 'Post-synaptic peers for E and I units' )

    subplot( 3, 3, 5 )
    [ bh, eh ] = barwerror( 1 : 4, npre( 1, : ), npre( 2, : ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number units' )
    title( 'Pre-synaptic peers for E (excitatory) and I (inhibitory) units' )
    
    subplot( 3, 3, 6 )
    [ bh, eh ] = barwerror( 1 : 4, abs( asg( 1, : ) ), asg( 2, : ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for E and I units' )

    subplot( 3, 3, 7 )
    [ bh, eh ] = barwerror( 1 : 2, asg( 1, 1 : 2 ), asg( 2, 1 : 2 ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for E units' )

    subplot( 3, 3, 8 )
    [ bh, eh ] = barwerror( 1 : 2, asg( 1, 3 : 4 ), asg( 2, 3 : 4 ) ...
        , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for I units' )
    
    %-----
    % example of punit and CCH and ACH
    
    if ismac
        filebase = filebaseLookup( 'mDL5', -16 );
    elseif ispc
        filebase = 'G:\mice\mDL5_16\mDL5_16';
    end
    fig3( 1 ) = figure;
    plot_ss( filebase, [ 1 6 ] ); % punit with 5 post-synaptic peers and no pre
    %mono.pairsExc( mono.pairsExc( :, 1 ) == 3, : )
    % punit is on S1 (1.6)
    % show CCH with a post-synaptic INT (S2; 2.11), and a post-synaptic PYR (S3; 3.52)
    % note that [ 3 50 ] is also a punit, with multiple pre-synaptic peers
    fig3( 2 ) = figure;
    plot_ss( filebase, [ 2 11 ] );
    fig3( 3 ) = figure;
    plot_ss( filebase, [ 3 52 ] );
    fig3( 4 ) = figure;
    subplot( 2, 2, 1 )
    plot_s2s( filebase, [ 1 6; 2 11 ], 'plotmode', -4 );
    subplot( 2, 2, 2 )
    plot_s2s( filebase, [ 1 6; 2 12 ], 'plotmode', -4 );
    subplot( 2, 2, 3 )
    plot_s2s( filebase, [ 1 6; 3 52 ], 'plotmode', -4 );
    
    % quantify the strength
	load( [ filebase '.s2s' ], '-mat', 's2s' )
    n12                         = [ 1 6; 2 11 ];
    [ g1, g2, act, sil, s, fig ] = calc_asg( s2s, n12, 'graphics', 1 );
    fig3( 5 ) = fig;
    
    n12                         = [ 1 6; 2 12 ];
    [ g1, g2, act, sil, s, fig ] = calc_asg( s2s, n12, 'graphics', 1 );
    fig3( 6 ) = fig;
    
    n12                         = [ 1 6; 3 52 ];
    [ g1, g2, act, sil, s, fig ] = calc_asg( s2s, n12, 'graphics', 1 );
    fig3( 7 ) = fig;
   
    %-----
    % save the figures
    fig = fig3;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG3_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 4
% firing rate statistics (fano factor)
% ----------------------------------------------------------------------

if fignums( 4 )
    
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
    
    % process the data
    
    
    % plot the figures
    %--------------------------------
    % ACH COM
    fig4( 1 ) = figure; 
    
    ispos                   = sst1.extremum > 0;
    nbins                   = 60;
    byprob                  = 0;
    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ach_com, ispos, nbins, colors_NP, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'Punits'; 
    lh.String{ 2 }          = 'Nunits';

    subplot( 2, 2, 2 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ach_com, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ACH-COM [ms]', byprob );
    lh.String{ 1 }          = 'INT'; 
    lh.String{ 2 }          = 'PYR';
    
    nbins                   = 30;
    byprob                  = 1;
    subplot( 2, 2, 3 )
    sst1i                   = struct_select( sst1, ispos );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1i.ach_com, true( sum( ispos ), 1 ), nbins, colors_NP, 'ACH-COM [ms]', byprob );
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( sst1.ach_com, sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'ACH-COM [ms]', byprob );
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
    title( sprintf( '%s: %0.2g; %s: %0.2g; %s: %0.2g', cnames{ 1 }, pvals( 1 ) ...
        , cnames{ 2 }, pvals( 2 ), cnames{ 3 }, pvals( 3 ) ) )
    
    %--------------------------------
    % FF - baseline rates
    fig4( 3 )               = figure;
    ispos                   = sst1.extremum > 0;
    nbins                   = 60;
    byprob                  = 0;
    subplot( 2, 2, 1 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( log2( sst1.baseline ), ispos, nbins, colors_NP, 'Log Firing rate [spk/s]', byprob );
    lh.String{ 1 }          = 'Punits'; 
    lh.String{ 2 }          = 'Nunits';
    lin2log( 'x', 2, 5 );

    subplot( 2, 2, 2 )
    [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( log2( sst1.baseline ), sst1.shankclu( :, 3 ) == 1, nbins, colors_PI, 'Log Firing rate [spk/s]', byprob );
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
    
    %--------------------------------
    % FF 
    % to do - look at individual images
    fig4( 2 )                   = figure;
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
    
    %-----
    % save the figures
    fig = fig4;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG4_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 5
% Place analyses
% ----------------------------------------------------------------------

if fignums( 5 )

    % load Hadas's data
    L                       = load( [ datadir_hadas fname_hadas '.mat' ], '-mat' );
    
    % partition into two separate structures
    fieldnames              = fields( L );
    nfields                 = length( fieldnames );
    n_all                   = size( getfield( L, 'all_session_tag' ), 1 );
    n_fields                = size( getfield( L, 'session_tag' ), 1 );
    for i                   = 1 : nfields
        
        nrow                = size( getfield( L, fieldnames{ i } ), 1 );
        if nrow == n_all
            eval( sprintf( 's_all.%s = L.%s;', fieldnames{ i }, fieldnames{ i } ) )
        elseif nrow == n_fields
            eval( sprintf( 's_field.%s = L.%s;', fieldnames{ i }, fieldnames{ i } ) )
        end
           
    end
    
    % add the extremum field from sst to the relevant items in s_all
    s_all.all_extremum      = NaN( n_all, 1 );
    u1                      = unique( sst.filebase );
    u2                      = unique( s_all.all_session_tag );
    u                       = unique( [ u1; u2 ] );
    [ ~, f1 ]               = ismember( sst.filebase, u );
    [ ~, f2 ]               = ismember( s_all.all_session_tag, u );
    s1                      = [ f1 sst.shankclu ];
    s2                      = [ f2 s_all.all_shankclu( :, 1 : 2 ) ];
    % if Hadas's data were unique, then we could write:
    %     [ ~, i1, i2 ]           = intersect( s1, s2, 'rows' );
    %     s_all.extremum( i2, : ) = sst.extremum( i1, : );
    % since it is not unique, we use the following:
    us2                     = unique( s2, 'rows' );
    nus2                    = size( us2, 1 );
    for i                   = 1 : nus2
        i1                  = ismember( s1, us2( i, : ), 'rows' );
        i2                  = ismember( s2, us2( i, : ), 'rows' );
        if sum( i1 ) > 0
            s_all.all_extremum( i2, : ) = repmat( sst.extremum( i1, : ), [ sum( i2 ) 1 ] );
        end
    end
    
    % do the same for the s_field:
    s_field.extremum        = NaN( n_fields, 1 );
    u2                      = unique( s_field.session_tag );
    u                       = unique( [ u1; u2 ] );
    [ ~, f2 ]               = ismember( s_field.session_tag, u );
    s2                      = [ f2 s_field.shankclu( :, 1 : 2 ) ];
    us2                     = unique( s2, 'rows' );
    nus2                    = size( us2, 1 );
    for i                   = 1 : nus2
        i1                  = ismember( s1, us2( i, : ), 'rows' );
        i2                  = ismember( s2, us2( i, : ), 'rows' );
        if sum( i1 ) > 0
            s_field.extremum( i2, : ) = repmat( sst.extremum( i1, : ), [ sum( i2 ) 1 ] );
        end
    end

    % combine back for Hadas
    s_combined              = struct_merge( s_all, s_field );
    savename                = [ datadir_hadas fname_hadas '_extremum.mat' ];
    save( savename, 's_combined' )
    
    % process the data
    eidx                    = s_all.all_shankclu( :, 3 ) == 1; % PYR
    iidx                    = s_all.all_shankclu( :, 3 ) == 0; % INT
    pidx                    = s_all.all_extremum > 0; % PUNIT
    eidx( pidx )            = 0; % NUNIT PYR 
    iidx( pidx )            = 0; % NUNIT INT
    colnames                = { 'bits/s', 'bits/spike' };
    row_names               = { 'PYR', 'INT', 'Punits' };
    info_meds               = [ median( s_all.all_info( eidx, : ) )
                                median( s_all.all_info( iidx, : ) )
                                median( s_all.all_info( pidx, : ) ) ];
    info_sem                = [ calc_sem( s_all.all_info( eidx, : ) )
                                calc_sem( s_all.all_info( iidx, : ) )
                                calc_sem( s_all.all_info( pidx, : ) ) ];
                            
    % plot the figures
    fig5( 1 )               = figure;
    for i = 1 : 2
        subplot( 2, 2, i )
        [ bh, eh ] = barwerror( 1 : 3, info_meds( :, i ), info_sem( :, i ) ...
            , colors_NP( 1, : ), 0.8, [ 0 0 0 ], 2 );
        set( gca, 'XTickLabel', row_names )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        ylabel( colnames{ i } )
        title( 'Spatial information' )
    end
    
    % do the same for many other interesting spatial features at the unit
    % and at the field level...
    
    fig5( 2 ) = figure;

    %-----
    % save the figures
    fig = fig5;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG5_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 6
% 
% ----------------------------------------------------------------------

if fignums( 6 )
    
    % process the data
    didx                = ~isnan( sst.depth );
    sst1                = struct_select( sst, didx );
    pidx                = sst1.extremum > 0;
    w                   = sst1.geo_fwhm;
    d                   = sst1.depth * 20; % [um]
    amp                 = abs( sst1.extremum );
    
    % open points:
    % 1. note that mC41 is 15 um spacing
    % 2. also note that in many cases, sst.depth was computed without flipping
    % - should be modified one by one (not necessarily for all)
    
    % determine spatial bins:
    binsize             = 20; 
    minVal              = floor( min( d ) / binsize ) * binsize;
    maxVal              = ceil( max( d ) / binsize ) * binsize;
    rside               = ( 0 + binsize / 2 ) : binsize : ( maxVal + binsize / 2 );
    lside               = ( 0 + binsize / 2 ) : binsize : ( -minVal + binsize / 2 );
    edges               = [ fliplr( -lside ) rside ];
    binC                = ( edges( 1 : end - 1 )  + edges( 2 : end ) )' / 2;
    
    % compute histograms
    hh                  = histc( d( ~pidx ), edges );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h0                  = hh;
    
    hh                  = histc( d( pidx ), edges );
    hh( end - 1 )       = hh( end - 1 ) + hh( end ); 
    hh( end )           = [];
    h1                  = hh;
    
    % plot the figures
    fig6( 1 )           = figure;
    subplot( 2, 2, 1 )
    barh( binC, h0, 1, 'EdgeColor', colors_NP( 1, : ), 'FaceColor', colors_NP( 1, : ) )
    xlabel( 'Nunit count' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 2, 2 )
    barh( binC, h1, 1, 'EdgeColor', colors_NP( 2, : ), 'FaceColor', colors_NP( 2, : ) )
    xlabel( 'Punit count' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 2, 3 )
    frcts               = h1 ./ ( h0 + h1 );
    minCount            = 2;
    frcts( ( h0 + h1 ) <= minCount ) = NaN;
    barh( binC, frcts, 1, 'EdgeColor', colors_NP( 2, : ), 'FaceColor', colors_NP( 2, : ) )
    xlabel( 'Fraction of Punits' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );

    % plot some statstic by depth
    fig6( 2 ) = figure;
    subplot( 2, 2, 1 )
    hold on, 
    ph( 1 ) = plot( w( ~pidx ), d( ~pidx ), '.' ); 
    ph( 2 ) = plot( w( pidx ), d( pidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NP( 1, : ) )    
    set( ph( 2 ), 'color', colors_NP( 2, : ) )    
    xlabel( 'FWHM' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    
    subplot( 2, 2, 2 )
    hold on, 
    ph( 1 ) = plot( amp( ~pidx ), d( ~pidx ), '.' ); 
    ph( 2 ) = plot( amp( pidx ), d( pidx ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_NP( 1, : ) )    
    set( ph( 2 ), 'color', colors_NP( 2, : ) )    
    xlabel( 'Ampplitude [\muV]' )
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
            figname = sprintf( '%s/punits_FIG6_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 7
% 
% ----------------------------------------------------------------------

if fignums( 7 )
    
    % process the data
    
    % plot the figures
    fig7( 1 ) = figure;
    
    fig7( 2 ) = figure;

    %-----
    % save the figures
    fig = fig7;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/punits_FIG7_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

return


%----------------------------------------------------
% internal functions
%----------------------------------------------------


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
