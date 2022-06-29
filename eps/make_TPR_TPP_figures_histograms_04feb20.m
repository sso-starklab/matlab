% make_jpn_figures_histograms      internal callback from make_jpn_figures
%
% call                              make_jpn_figures_histograms( param, id, xstr, tstr, param_type)
%
% calls                             circ_hist (circs)
%                                   lin2log (graph)

% 31-dec-19 HS and ES

% revisions
% 01-jan-20 added statistical testing (KW and multiple comparisons)
%           added medians
% 21-jan-20 adapted for JPN project
% 04-Feb-20 readapted for JPN project

function make_TPR_TPP_figures_histograms( param, id, tstr, param_type, varargin)

[ binsize, edges, color_gro, xlimits,xtext,ytext,binc, mass_cent, alpha]        = ParseArgPairs (...
    {'binSize','edges', 'color_gro', 'xlimits','xtext','ytext','binc', 'mass_cent', 'alpha'} ...
    ,{ 0.1, [], [], [], [], [], [], 0, 0.05}, varargin{ : } );

colors                  = [106,27,154;  46,125,50; 173, 20, 87;  0 0.7 0; 0, 0.75, 0.75]/255;

if ~isempty(color_gro)
    colors              = colors (color_gro,:);
end

% preps

switch param_type( 1 : 3 )
    case 'log'
        val         = log2( param );
        if isempty(edges)
            edges       = ( floor( min( val ) ) - binsize/2 ) : binsize : ( ceil( max( val ) ) + binsize/2 );      
        end
        
    case 'lin'
        val         = param;
        if isempty(edges)
            edges       = floor( min( val ) / binsize ) * binsize : binsize : ceil( max( val ) / binsize ) * binsize;
        end
    case 'cir'
        val         = mod( param, 2 * pi );
end


ntypes              = size ( id, 2);
data                = cell (1, ntypes);
n                   = zeros (1, ntypes);
g                   = zeros( size( val ) );

for i = 1 : ntypes
    data{i}         = val ( id(:,i) );
    n(i)            = sum ( id(:,i) );
    g( id(:,i) )    = i;
end


%statistical test
if isequal( param_type( 1 : 3 ), 'cir' ) 
    
    if ntypes == 2
    [ ~, pval ] = wheeler( val( id(:,1) ), val( id(:,2) ), alpha );
    else
        pval = [];
    end
    
elseif ntypes == 2 
    
    pval = ranksum( val( id(:,1) ), val( id(:,2) ) ); % Mann-Whitney U-test
    
else    
    
    [ pval, ~, stats ] = kruskalwallis( val, g, 'off' ); % kruskal wallis
    cc                 = multcompare( stats, 'display', 'off' ); %post hoc multi comparisons
    
end

% computations
hh              = cell (1, ntypes);

if isequal( param_type( 1 : 3 ), 'cir' )
    
    nbins = 10;
    for i = 1 : ntypes
        [ hh{i}, ~, edges ] = circ_hist( data{ i }, nbins );
    end
    if isempty(binc)
        binc                = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
    end
else
    
    if isempty(binc)
        binc                = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
    else
        binSize             = binc(2)-binc(1);
        edges               = binc(1)-binSize/2 : 1 :  binc(end) + binSize/2;
    end
    
    for i = 1 : ntypes
        h               = histc( data{ i }, edges );
        h( end - 1 )    = sum( h( end - 1 : end ) );
        h( end )        = [];
        hh{i}           = h;
    end
    
end
% Each bin includes the left endpoint, but does not include the right endpoint.

% plot
newplot
for i           = 1 : ntypes
    
    if isempty ( hh{i} )
        continue;
    end
    
    subplot( 2, ceil(ntypes/2), i )
    switch param_type( 1 : 3 )
        case 'log'
            lin2log( 'x', 2, 1 );
            med = nanmedian( data{ i } );
            medstr = 2.^med;
            bh          = bar( binc, hh{ i }, 1 );
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
            alines( med, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
            title( sprintf( 'n=%d, med=%0.2g', n(i), medstr) )
            a = 2.^binc;
            if length(binc) < 10
                j = 2;
            elseif length(binc) < 15
                j = 10;
            else
                j = 15;
            end
            set( gca, 'xtick', binc(1:j:end))
            set( gca, 'XTickLabel',round(a(1:j:end),3) )

        case 'cir'
            set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
            set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
            
            
            bh      = polarhistogram(data{ i },edges, 'normalization', 'probability');
            hold on;
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
            med    = wrapTo2Pi ( circ_med( data{ i } ) );
            medstr = med;
            maxcount = max(bh.Values);
            polarplot( [medstr; medstr], [0 0 0; 1 1 1]*maxcount, 'r', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            title( sprintf( 'n=%d, med=%0.3g%c', n(i), medstr/pi,  char(960)) )

        case 'lin'
            bh          = bar( binc, hh{ i }, 1 );
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
            if ~mass_cent
                med = nanmedian( data{ i } );
                alines( med, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
                title( sprintf( 'n=%d, med=%0.2g', n(i), med) )
            else
                A   = data{ i };
                com = calc_com( A, ones(size(A)) );
                com = com(1);
                alines( com, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
                title( sprintf( 'n=%d, center of mass=%0.2g', n(i), com) )
            end
            
    end
    set( gca, 'tickdir', 'out', 'box', 'off' );
    if ~isempty(xlimits)       
        xlim(xlimits);
    end
    if ~isempty(xtext)
        xlabel(xtext);
    end
    if ~isempty(ytext)
        ylabel(ytext);
    end
end



suptitle(tstr)
textf( 0.5, 0.975, sprintf( 'p<=%0.3g', pval ) );

return

% EOF
