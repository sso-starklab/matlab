% make_punits_figures           make figures for positive units paper

% 17-may-20 SSo + ES

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
        datadir = 'D:\_Shirly\Lab\eps';
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
    
    % plot the figures
    
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

