% hfoAnalysisSpikingPlot1       call-back from hfoAnalysisSpiking
%
% call                          fig = hfoAnalysisSpikingPlot1( s, r )
%
% calls                         (circs)       circ_mean, circlin, wheeler
%                               (general)     calc_gain, inrange, make_equal_bins, ParseArgPairs, pol2car, remrnd, replacetok, scaleto, sortcols, uhist
%                               (graph)       addmargins, alines, barwerror, circ, fig_out, hcolorbar, myjet, openbarh, patch_band, scatterwerror, textf
%                               (sets)        parse
%                               (ssp)         firfilt, makegaussfir
%                               (structs)     struct_select
%                               (stats)       bounds, calc_sem, calc_spearman, fisher_exact_test, gtest, moutliers, utest
%
% calls                         (general)     calc_gain, ParseArgPairs, replacetok, scaleto, sortcols
%                               (graph)       alines, imagescbar, myjet, patch_band, textf
%                               (ssp)         firfilt, makegaussfir
%                               (stats)       moutliers
%
% see also                      hfoAnalysisSpiking

% 30-oct-13 ES

% revisions
% 24-mar-21 cleaned up

function fig = hfoAnalysisSpikingPlot1( s, r, varargin )

fig                             = [];

%----------------------------------------------------------------------
% constants
%----------------------------------------------------------------------
sepStyle                        = '--';
blackColor                      = [ 0 0 0 ];
whiteColor                      = [ 1 1 1 ];
colors                          = [ 0 0 0.7; 1 0 0 ];

%----------------------------------------------------------------%
% arguments
%----------------------------------------------------------------%
nargs                           = nargin;
if nargs < 1 || isempty( s )
    return
end
[ toplot ...
    , figname, minEvents, trigMode, exttitle ...
    , smoothphase, pTH ...
    ]                           = ParseArgPairs(...
    { 'toplot' ...
    , 'figname', 'minEvents', 'trigMode', 'exttitle' ...
    , 'smoothphase', 'pTH' ...
    }...
    , { 'gain' ...
    , '', 100, 'induced', '' ...
    , 1, 0.05 ...
    }...
    , varargin{ : } );
if ~isempty( figname )
    [ ~, nm, ext ]              = fileparts( figname );
    tstr                        = [ nm ext ];
else
    tstr                        = '';
end
ufnames                         = unique( s.filebase );
if isempty( tstr ) && length( ufnames ) == 1
    tstr                        = sprintf( '%s', ufnames{ 1 } );
end
tstr                            = [ tstr ' ' trigMode ];
tstr                            = replacetok( tstr, '\_', '_' );
if ~isempty( exttitle )
    if isa( exttitle, 'cell' )
        exttitle                = exttitle{ 1 };
    end
    tstr                        = exttitle;
end

%----------------------------------------------------------------------
% preps
%----------------------------------------------------------------------
if smoothphase
    win                         = makegaussfir( 1, 1 );                     % smooth the histograms
else
    win                         = 1;
end

% remove units with too few events:
kidx                            = s.nevents >= minEvents;
ridx                            = ismember( s.trigfname, unique( s.filename( ~kidx ) ) );
s                               = struct_select( s, kidx );
s                               = struct_select( s, ~ridx );
if exist( 'r', 'var' ) && ~isempty( r )
    r( ridx )                   = [];
end

if nansum( s.phsHists( : ) ) == 0
    return
end

% determine cycle range to be plotted:
gmean                           = nanmean( scale( s.phsHists' ), 2 );
TH                              = 0.05 * max( gmean );
lims                            = s.phsBins( [ find( gmean > TH, 1, 'first' ) find( gmean > TH, 1, 'last' ) ] );
cyclims                         = ceil( max( abs( lims ) ) / ( 2 * pi ) ) * [ -1 1 ] + [ 1 -1 ];

pidx                            = s.shankclu( :, 3 );
nclu                            = length( pidx );

phshists                        = s.phshists';
phsbins                         = s.phsbins( 1, : );
ratesBL                         = s.ratesBL;
phsHists                        = s.phsHists';
phsBins                         = s.phsBins( 1, : );
pval                            = s.pval;
count                           = s.denom;
nbins                           = length( phsbins );

ncycs                           = ( length( s.phsBins ) - 1 ) / nbins;
for i                           = 1 : nclu
    rcycs( i, : )               = nanmean( reshape( count( i, 1 : end -1 ), [ nbins ncycs ] ) );
end
if ~isempty( tstr )
    tstr                        = sprintf( '%s; %d events, %d sec', tstr, length( s.trigs ), round( sum( diff( s.periods, 1,2  ) + 1 )/s.Fs( 1 ) ) );
end
if trigMode( 1 ) == 'i'
    cyclims                     = [ -1 9 ];
elseif isequal( trigMode, 'spontaneous' )
    cyclims                     = [ -4 5 ];
end

% expand phsBinSize for backward compatibility
if size( s.phsBinSize, 1 ) ~= nclu
    s.phsBinSize                = s.phsBinSize * ones( nclu, 1 );
end

%----------------------------------------------------------------------
% plot the phase histograms
fig                             = figure;
for ct                          = 0 : 1
    cidx                        = pidx == ct & s.pval <= pTH;
    if sum( cidx ) == 0
        continue
    end
    
    yticks                      = 1 : sum( cidx );
    
    % collapsed on a single cycle:
    subplot( 3, 2, 3 + 2 * ct ),
    xticks                      = [ phsbins( : ) - 2 * pi; phsbins( : ); phsbins( : ) + 2 * pi ];
    z0                          = phshists( :, cidx );
    z0                          = bsxfun( @rdivide, z0, nansum( rcycs( cidx, : ), 2 )' ); % convert to counts/event
    switch toplot
        case 'gain'
            z0                  = z0 ./ ( ones( nbins, 1 ) * s.phsBinSize( cidx )' );
            z0                  = calc_gain( z0, ones( size( z0, 1 ), 1 ) * ratesBL( cidx )', 3 );
        case 'rate'
            z0                  = z0 / phsBinSize;                          % convert to rates
        case 'occurrence'
            z0                  = z0 * 100;                                 % convert to occurrence (percent)
    end
    
    [ z0, sidx ]                = sortcols( z0 );
    z0                          = repmat( z0, [ 3 1 ] );
    z0                          = firfilt( z0, win );
    [ hh, ah ]                  = imagescbar( xticks, yticks, z0 );
    
    subplot( ah( 1 ) ),
    xlabel( 'Phase [rad]' )
    xlim( [ -1 1 ] * 2 * pi )
    line( xticks, scaleto( cos( xticks ), ylim ), 'color', whiteColor );
    alines( ( -1 : 1 ) * pi, 'x', 'color', blackColor, 'linestyle', sepStyle );
    set( gca, 'xtick', ( -2 : 2 ) * pi, 'xticklabel', remrnd( ( -2 : 2 ) * pi, 0.01 ) );
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    subplot( ah( 2 ) )
    switch toplot
        case { 'gain', 'rate' }
            set( gca, 'xscale', 'log', 'yticklabel', '' )
            axis tight
            alines( [ 1 10 100 ], 'x', 'color', blackColor, 'linestyle', sepStyle );
        case 'occurrence'
            alines( mean( max( z0 ) ), 'x', 'color', blackColor, 'linestyle', sepStyle );
    end
    set( hh( 2 ), 'EdgeColor', colors( ct + 1, : ), 'FaceColor', colors( ct + 1, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    % summarize the histograms:
    subplot( 3, 2, 1 )
    mm                      	= moutliers( z0, @nanmean, 3, [], 2, 2 );
    ss                          = moutliers( z0, @calc_sem, 3, [], 2, 2 );
    xx                          = xticks;
    mm( isnan( mm ) )           = 0;
    ss( isnan( ss ) )           = 0;
    patch_band( xx, mm, ss, colors( ct + 1, : ) );
    if ct == 1
        xlabel( 'Phase [rad]' )
        xlim( [ -1 1 ] * 2 * pi )
        line( xticks, scaleto( cos( xticks ), ylim ), 'color', blackColor );
        alines( ( -1 : 1 ) * pi, 'x', 'color', blackColor, 'linestyle', sepStyle );
        alines( 1, 'y', 'color', blackColor, 'linestyle', sepStyle );
        set( gca, 'xtick', ( -2 : 2 ) * pi, 'xticklabel', remrnd( ( -2 : 2 ) * pi, 0.01 ) );
        set( gca, 'tickdir', 'out', 'box', 'off' )
        ylabel( toplot )
        pos                     = get( gca, 'position' );
        pos0                    = get( ah( 1 ), 'position' );
        set( gca, 'position', [ pos( 1 : 2 ) pos0( 3 ) pos( 4 ) ] );
    end
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    % cycle-resolved plots:
    subplot( 3, 2, 4 + 2 * ct )
    z0                          = phsHists( :, cidx );
    z0                          = z0 ./ count( cidx, : )';
    switch toplot
        case 'gain'
            num                 = ( s.phsHists( cidx, : ) ./ s.denom( cidx, : ) ) ./ ( s.phsBinSize( cidx ) * ones( 1, size( s.denom, 2 ) ) );
            den                 = s.ratesBL( cidx ) * ones( 1, size( s.phsHists, 2 ) );
            z0                  = calc_gain( num, den, 5 )';
        case 'rate'
            z0                  = z0 / phsBinSize;                          % rate [spikes/s]
        case 'occurrence'
            z0                  = z0 * 100;                                 % occurrence (percent of cycles)
    end
    z0                          = z0( :, sidx );                            % keep same order as in the collapsed display
    z0                          = firfilt( z0, win );
    xticks                      = phsBins( : ) / ( 2 * pi );
    xidx                        = inrange( xticks, cyclims );
    [ hh, ah ]                  = imagescbar( xticks( xidx ) , yticks, z0( xidx, : )  );
    subplot( ah( 1 ) )
    xlim( cyclims )
    alines( cyclims( 1 ) : cyclims( 2 ), 'x', 'color', blackColor, 'linestyle', sepStyle );
    xlabel( 'Cycle number' );
    subplot( ah( 2 ) )
    switch toplot
        case { 'gain', 'rate' }
            set( gca, 'xscale', 'log', 'yticklabel', '' )
            axis tight
            alines( [ 1 10 100 ], 'x', 'color', blackColor, 'linestyle', sepStyle );
        case 'occurrence'
            alines( mean( max( z0 ) ), 'x', 'color', blackColor, 'linestyle', sepStyle );
    end
    set( hh( 2 ), 'EdgeColor', colors( ct + 1, : ), 'FaceColor', colors( ct + 1, : ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    % summarize the running plots:
    subplot( 3, 2, 2 )
    mm                          = moutliers( z0, @nanmean, 3, [], 2, 2 );
    ss                          = moutliers( z0, @calc_sem, 3, [], 2, 2 );
    xx                          = phsBins / ( 2 * pi );
    mm( isnan( mm ) )           = 0;
    ss( isnan( ss ) )           = 0;
    patch_band( xx, mm, ss, colors( ct + 1, : ) );
    axis tight
    if ceil( max( ylim ) ) > 0
        ylim( [ 0 ceil( max( ylim ) ) ] )
    end
    if ct == 1
        xlim( cyclims )
        alines( cyclims( 1 ) : cyclims( 2 ), 'x', 'color', blackColor, 'linestyle', sepStyle );
        alines( 1, 'y', 'color', blackColor, 'linestyle', sepStyle );
        xlabel( 'Cycle number' );
        ylabel( toplot )
        pos                     = get( gca, 'position' );
        pos0                    = get( ah( 1 ), 'position' );
        set( gca, 'position', [ pos( 1 : 2 ) pos0( 3 ) pos( 4 ) ] );
    end
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
end
colormap( myjet )
textf( 0.5, 0.975, tstr );

return

% EOF

