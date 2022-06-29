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

[ binsize, edges, color_gro, xlimits, ylimits, xtext,ytext,binc, mass_cent, ...
    alpha, ratio, rotateXY, nNaNs, raotestF, raytestF, idlabels, Xticks, YR]        = ParseArgPairs (...
    {'binSize','edges', 'color_gro', 'xlimits', 'ylimits', 'xtext','ytext','binc', 'mass_cent', ...
    'alpha', 'ratio', 'rotateXY', 'nNaNs','raotestF','raytestF', 'idlabels','Xticks', 'YR'} ...
    ,{ 0.1, [], [], [], [], [], [], [], 0, ...
    0.05, 0, 0, 1, 0, 1, [], 0, []}, varargin{ : } );

colors          = [0.4940, 0.1840, 0.5560; ...
    0.1840,0.5000,0.2000;...
    0.9290 0.6940 0.1250; ...
    0.3010 0.7450 0.9330; ...
    0.8500 0.3250 0.0980;
    0.369   0.573    0.953; ...
    0.3529    0.8588    0.7020;
    ];

cc = [];

if ~isempty(color_gro)
    colors              = colors (color_gro,:);
end

% preps

switch param_type( 1 : 3 )
    case 'log'
        if ~isempty(xlimits)
            xlimits(1) = log2(xlimits(1));
            xlimits(2) = log2(xlimits(2));
        end
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
    
    if ~isempty(YR)
        cond        = val >= YR(1) & val <= YR(2);
        id(:,i)     = id(:,i) & cond;
    end
    
    data{i}         = val ( id(:,i) );
    n(i)            = sum ( id(:,i) );
    g( id(:,i) )    = i;
    
    if ~nNaNs
        n(i)        = sum ( id(:,i) & ~isnan( val ) );
    end
end


pval = [];

%statistical test
if isequal( param_type( 1 : 3 ), 'cir' )
    
    if ntypes == 2
        [ ~, pval ] = wheeler( val( id(:,1) ), val( id(:,2) ), alpha );
    else
        pval = [];
    end
    
elseif ntypes == 2
    
    pval = ranksum( val( id(:,1) ), val( id(:,2) ) ); % Mann-Whitney U-test
    
elseif ntypes > 2
    
    val                = val(g ~= 0);
    g                  = g(g ~= 0);
    [ pval, ~, stats ] = kruskalwallis( val, g, 'off' ); % kruskal wallis
    cc                 = multcompare( stats, 'display', 'off' ); %post hoc multi comparisons
%     
end

% computations
hh              = cell (1, ntypes);

if isequal( param_type( 1 : 3 ), 'cir' )
    
    nbins = 20;
    for i = 1 : ntypes
        if param_type=='circ'
            [ hh{i}, ~, ~ ] = circ_hist( data{ i }, nbins );
        else
        [ hh{i}, ~, edges ] = circ_hist( data{ i }, nbins );
        end
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
        N               = histcounts (data{i}, edges);
        if ratio
            N           = N / ( sum(N) );
        end
        hh{i}           = N;
    end
    
end
% Each bin includes the left endpoint, but does not include the right endpoint.
if ntypes > 1
    nrows               = 2;
else
    nrows               = 1;
end


% plot
for i           = 1 : ntypes
    
    if isempty ( hh{i} )
        continue;
    end
    
    subplot( nrows, ceil(ntypes/2), i )
    switch param_type( 1 : 3 )
        case 'log'
            
            lin2log( 'x', 2, 1 );
            med = nanmedian( data{ i } );
            medstr = 2.^med;
            bh          = bar( binc, hh{ i }, 1 );
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
            alines( med, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
            if ~isempty(idlabels)
                title( sprintf( '%s: n=%d, med=%0.3g', idlabels{i}, n(i), medstr) )
            else
                title( sprintf( 'n=%d, med=%0.3g', n(i), medstr) )
            end
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
            med    = wrapTo2Pi ( circ_med ( data{ i } ) );
            [phi,R]      = circ_mean (data{ i } ) ;
            phi = wrapTo2Pi(phi);
            maxcount = max(bh.Values);
            polarplot( [phi; phi], [0 R*maxcount], 'k', 'lineWidth',2);
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            tit = sprintf( 'n=%d, med=%0.3g%c, R=%.2f', n(i), phi/pi,  char(960), R);
            if raotestF
                [~, rao_p] = rao_test( data{ i });
                tit    = sprintf('%s, RAOp=%.2f', tit,rao_p);
            end
            if raytestF
                ray_p = ray_test( data{ i });
                tit    = sprintf('%s, RAYp=%.2f', tit,ray_p);
            end
            title( tit )
            
        case 'lin'
            
            if ~rotateXY
                bh          = bar( binc, hh{ i }, 1 );
                medaxis     = 'x';
            else
                bh          = barh( binc, hh{ i }, 1 );
                medaxis     = 'y';
            end
            
            
            set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
            if ~mass_cent
                med = nanmedian( data{ i } );
                range(1) = prctile(data{i},25) ;
                range(2) = prctile(data{i},75) ;
                SEM      = calc_sem(data{i});
                alines( med, medaxis, 'color', [ 0 0 0 ], 'linestyle', '--' );
                if ~isempty(idlabels)
                    title( sprintf( '%s: n=%d, med=%0.2g[%.2g-%.2g]', idlabels{i}, n(i), med, range(1), range(2)) )
                else
                    title( sprintf( 'n=%d, med=%0.2g[%.2g-%.2g]', n(i), med, range(1), range(2)) )
                    
                end
            else
                A   = data{ i };
                [com, SD] = calc_com( A, ones(size(A)) );
                SEM       = SD / sqrt(n(i));
                com = com(1);
                alines( com, medaxis, 'color', [ 0 0 0 ], 'linestyle', '--' );
                if ~isempty(idlabels)
                    title( sprintf( '%s: n=%d, COM=%0.2f, SEM=%0.2f', idlabels{i}, n(i), com,SEM) )
                else
                    title( sprintf( 'n=%d, COM=%0.2f, SEM=%0.2f', n(i), com, SEM) )
                end
            end
            
    end
    set( gca, 'tickdir', 'out', 'box', 'off' );
    if ~isempty(xlimits)
        xlim(xlimits);
    end
    if~isempty(ylimits)
        ylim(ylimits);
    end
    if ~isempty(xtext)
        xlabel(xtext);
    end
    if ~isempty(ytext)
        ylabel(ytext);
    end
    
    
    if   ~isequal(param_type(1:3),'cir')
        axis square
    end
    
    if Xticks
        xticks(round(binc,2))
    end
    
end



suptitle(tstr)

if ~isempty(pval)
    textf( 0.5, 0.975, sprintf( 'p<=%0.3g', pval ) );
end

if ~isempty(cc)
    
    ncomp       = size(cc,1);
    
    fprintf('\n');
    
    for i = 1 : ncomp
        
        if ~isempty(idlabels)
            fprintf('%s vs. %s: %.3f \n', idlabels{cc(i,1)}, idlabels{cc(i,2)}, cc(i,end));
        else
            fprintf('%d vs. %d: %.3f \n', cc(i,1), cc(i,2), cc(i,end));
            
        end
        
    end
    fprintf('\n');
    
    
end

return

% EOF
