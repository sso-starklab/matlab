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

function make_TPR_TPP_histograms2( param, id, tstr, param_type, varargin)

[ binsize, edges, color_gro, xlimits, round_factor, binc, xtext, ytext, ...
    grps, mass_cent, ratio, ylimits, xstr]        = ParseArgPairs (...
    {'binSize','edges', 'color_gro', 'xlimits', 'round_factor', 'binc', 'xtext', 'ytext', ...
    'grps', 'mass_cent', 'ratio', 'ylimits', 'xstr'} ...
    ,{ 0.1, [], [], [], 3, [], [], [], ...
    [], 0, 0, [], []}, varargin{ : } );

colors                  = [106,27,154;  46,125,50; 173, 20, 87;  0 0.7 0; 0, 0.75, 0.75]/255;

if ~isempty(color_gro)
    colors              = colors (color_gro,:);
end
ntypes              = size ( id, 2);

% get minmax vals
min_val             = inf;
max_val             = -inf;
if ~isempty(xlimits)
    min_val             = xlimits(1);
    max_val             = xlimits(2);
else
    for i = 1 : ntypes
        val             = param(id(:,i));
        min_val_t       = min(val);
        max_val_t       = max(val);
        min_val         = min([min_val, min_val_t]);
        max_val         = max([max_val, max_val_t]);
    end
end

% preps
switch param_type( 1 : 3 )
    case 'log'
        val         = log2( param );
        min_val     = floor( log2(min_val) );
        max_val     = ceil( log2(max_val) );
        min_val2     = min_val - binsize/2;
        max_val2     = max_val + binsize/2;
        if isempty(edges)
            edges       = min_val2 : binsize : max_val2;      
        else
            min_val     = edges(1);
            max_val     = edges(end);
%             edges       = log2(edges);
        end
        
    case 'lin'
        val         = param;
        min_val2     = floor( min_val / binsize ) * binsize;
        max_val2     = ceil( max_val / binsize ) * binsize;
        if isempty(edges)
            edges       = min_val2 : binsize : max_val2;
        end
        if ~isempty(binc)
            binsize    = binc(2) -binc(1);
            edges      = (binc(1) - binsize/ 2 ) : binsize : (binc(end) + binsize / 2);
        end
    case 'cir'
        val         = mod( param, 2 * pi );
end


data                = cell (1, ntypes);
n                   = zeros (1, ntypes);
g                   = zeros( size( val ) );

for i = 1 : ntypes
    data{i}         = val ( id(:,i) );
    n(i)            = sum ( id(:,i) );
    g( id(:,i) )    = i;
end

%statistical test
if ntypes == 2 
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
    
else
        
    for i = 1 : ntypes
        h               = histc( data{ i }, edges );
        h( end - 1 )    = sum( h( end - 1 : end ) );
        h( end )        = [];
        if ratio
            tot             = sum(h);
            h               = h/tot;
        end
        hh{i}           = h;
    end
    
end

if isempty(binc)
    binc                = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
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
            a = 2.^binc;
            if length(binc) < 10
                j = 1;
            elseif length(binc) < 15
                j = 10;
            else
                j = 15;
            end
            set( gca, 'xtick', binc(1:j:end))
            set( gca, 'XTickLabel',round(a(1:j:end),round_factor) )
            t = sprintf( sprintf( 'n=%d, med=%0.2g', n(i), medstr) );

        case 'cir'
            set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
            set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )
            
            
            bh      = polarhistogram(data{ i },edges);
            hold on;
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
            med    = wrapTo2Pi ( circ_med( data{ i } ) );
            medstr = med/pi;
            maxcount = max(bh.Values);
            polarplot( [medstr; medstr], [0 0 0; 1 1 1]*maxcount, 'r', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            t = sprintf( sprintf( 'n=%d, med=%0.2g', n(i), medstr) );
            
        case 'lin'
            bh          = bar( binc, hh{ i }, 1 );
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );

            if ~mass_cent
                medstr = nanmedian( data{ i } );
                alines( medstr, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
                t = sprintf( sprintf( 'n=%d, med=%0.2g', n(i), medstr) );
            else
                com = calc_com( data{ i }, ones(size(data{ i })) );
                com = com(1);
                alines( com, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
                t= sprintf( sprintf( 'n=%d, center of mass=%0.2g', n(i), com) );
            end
            
            
    end
    
    if ~isempty(xstr)
                t   = sprintf( '%s, %s', xstr{i}, t);
    end
    title(t);
    set( gca, 'tickdir', 'out', 'box', 'off' );
    xlim([min_val, max_val]);
    
    if ~isempty(ylimits)
        ylim(ylimits);
    end
    
    if ~isempty(xtext)
        xlabel(xtext);
    end
    if ~isempty(ytext)
        ylabel(ytext);
    end
        
end


    
if ~isempty(grps) && ntypes > 2
    rm_idx                    = cc(:,2) == ntypes + 1;
    cc(rm_idx,:)              = [];
    mult_comp_text            = cell(size(cc,1),1);
    s                         = '';
    
    for i = 1 :  size(cc,1) 
        mult_comp_text{i}     = sprintf( '%s vs. %s: %0.2g', grps{ cc( i, 1 ) }, grps{ cc( i, 2 ) }, cc( i, 6 ) );
        s                     = sprintf('%s; %s', s, mult_comp_text{i});
    end
else
    s                         = [];
end

suptitle(tstr)
textf( 0.5, 0.975, sprintf( 'p<=%0.3g%s', pval, s ) );


return

% EOF
