% make_punits_figures           make figures for positive units paper

% 17-may-20 SSo + ES

% revisions
% 31-may-20 (1) completed check_mono derivation of connectivity maps
%           (2) updated figure 3, part 8
% 01-june-20(1) updated datadir and outdir for ispc
% 02-jun-20 (2) added extraction of ASG-e and ASG-i from s2s
%           (3) added plotting of ASG
%           (4) added significance testing for ASG, nPRE, nPOST

% next:
% update figure 3, part 8, with ASG and STC statistics

% figure 1: 
% figure 2: 
% figure 3: connectivity
% figure 4: 
% figure 5: 
% figure 6: 
% figure 7:

function make_punits_figures( fignums, savef, savetype, outdir )

% general constants
Nfigs                   = 10;
colors                  = [ 0 0 0.7; 1 0 0 ];

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
        outdir = '/Users/eranstark/Documents/graphics/punits/data_figures';
        datadir = '/Users/eranstark/Documents/da/punits/';
        % assumes that data are synchronized with odin, if not, write:
        % rsync -avh --progress /Volumes/odin/Shirly/eps/fanofactor/* ~/Documents/da/punits/fanofactor/
        % rsync -avh --progress /Volumes/odin/Shirly/eps/sps/* ~/Documents/da/punits/sps/ 
        % rsync -avh --progress /Volumes/odin/Shirly/eps/sst/* ~/Documents/da/punits/sst/
        % rsync -avh --progress /Volumes/odin/Shirly/eps/s2s/*  ~/Documents/da/punits/s2s/ 
    elseif isunix   % shirly (nanna) Linux
    elseif ispc     % shirly (nanna) Windows
        outdir = 'D:\_Shirly\Lab\eps\data_figures';
        datadir = 'G:\mice\EPS\';
    end
end

% get the database
[ res, sst ]            = shirly_eps_analysis( [], [], [], 1 );
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
% 
% ----------------------------------------------------------------------

if fignums( 2 )
    
    % process the data
    
    % plot the figures
    fig2( 1 ) = figure;
    
    fig2( 2 ) = figure;

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
    H_JB = [ jbtest( s.nexcPost( eidx & ~pidx ) ), jbtest( s.nexcPost( eidx & pidx ) ) ...
    jbtest( s.ninhPost( iidx & ~pidx ) ), jbtest( s.ninhPost( iidx & pidx ) ) ];
    if all( H_JB == 0 )
        [ ~, pp1 ]  = ttest2( s.nexcPost( eidx & ~pidx ) , s.nexcPost( eidx & pidx ) );
        [ ~, pp2 ]  = ttest2( s.ninhPost( iidx & ~pidx ), s.ninhPost( iidx & pidx ) );
    else
        pp1         = utest( s.nexcPost( eidx & ~pidx ), s.nexcPost( eidx & pidx ) );
        pp2         = utest( s.ninhPost( iidx & ~pidx ), s.ninhPost( iidx & pidx ) );
    end
    pval_NPOST   = [ pp1 pp2 ];
    
    
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
    %bh = bar( 1 : 4, [ emat( 3, : ) imat( 3, : ) ], 0.8 );
    %set( bh, 'FaceColor', [ 1 0 1 ], 'EdgeColor', [ 1 0 1 ] )
    [ bh, eh ] = barwerror( 1 : 4, [ emat( 3, [ 2 1 ] ) imat( 3, [ 2 1 ] ) ], [ emat( 4, [ 2 1 ] ) imat( 4, [ 2 1 ] )] ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'E-nunits', 'E-punits', 'I-nunits', 'I-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Units that are excitatory (E) or inhibitory (I), p=' )
    
    subplot( 3, 3, 2 )
    %bh = bar( 1 : 4, [ emat( 3, : ) imat( 3, : ) ], 0.8 );
    %set( bh, 'FaceColor', [ 1 0 1 ], 'EdgeColor', [ 1 0 1 ] )
    [ bh, eh ] = barwerror( 1 : 4, [ emat2( 3, [ 2 1 ] ) imat2( 3, [ 2 1 ] ) ] ...
        , [ emat2( 4, [ 2 1 ] ) imat2( 4, [ 2 1 ] ) ] ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Units that receive excitation or inhibition' )
    
    subplot( 3, 3, 3 )
    %bh = bar( 1 : 4, [ emat( 3, : ) imat( 3, : ) ], 0.8 );
    %set( bh, 'FaceColor', [ 1 0 1 ], 'EdgeColor', [ 1 0 1 ] )
    [ bh, eh ] = barwerror( 1 : 2, cmat( 3, [ 2 1 ] ), cmat( 4, [ 2 1 ] ) ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Connected-nunits', 'Connected-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Fraction of units' )
    title( 'Connected units' )
    
    subplot( 3, 3, 4 )
    [ bh, eh ] = barwerror( 1 : 4, npost( 1, : ), npost( 2, : ) ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excitatory-nunits', 'Excitatory-punits', 'Inhibitory-nunits', 'Inhibitory-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number of units' )
    title( 'Post-synaptic peers for E and I units' )

    subplot( 3, 3, 5 )
    [ bh, eh ] = barwerror( 1 : 4, npre( 1, : ), npre( 2, : ) ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'Number units' )
    title( 'Pre-synaptic peers for E (excitatory) and I (inhibitory) units' )
    
    subplot( 3, 3, 6 )
    [ bh, eh ] = barwerror( 1 : 4, abs( asg( 1, : ) ), asg( 2, : ) ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits', 'Inhibited-nunits', 'Inhibited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for E and I units' )

    subplot( 3, 3, 7 )
    [ bh, eh ] = barwerror( 1 : 2, asg( 1, 1 : 2 ), asg( 2, 1 : 2 ) ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
    colnames = { 'Excited-nunits', 'Excited-punits' };
    set( gca, 'XTickLabel', colnames )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( 'ASG' )
    title( 'Mean ASG for E units' )

    subplot( 3, 3, 8 )
    [ bh, eh ] = barwerror( 1 : 2, asg( 1, 3 : 4 ), asg( 2, 3 : 4 ) ...
        , [ 1 0 1 ], 0.8, [ 0 0 0 ], 2 );
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
% 
% ----------------------------------------------------------------------

if fignums( 4 )
    
    % process the data
    
    % plot the figures
    fig4( 1 ) = figure;
    
    fig4( 2 ) = figure;

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
% 
% ----------------------------------------------------------------------

if fignums( 5 )
    
    % process the data
    
    % plot the figures
    fig5( 1 ) = figure;
    
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
    
    % plot the figures
    fig6( 1 ) = figure;
    
    fig6( 2 ) = figure;

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

% EOF

