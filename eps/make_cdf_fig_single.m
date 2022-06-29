% prm - values
% cond - matrix [n m] n length of prm
%                     m number of group
% tstr - title
% xtext - x-axes legend
% ytext - y-axes legend
% ticks - custum c ticks
% xstr - group names in legend
% com - calculation of center od mass instead of median

function make_cdf_fig_single(prm, cond, varargin)



[ tstr, xtext, ytext, ticks, xstr, xlimits, ...
    logscale, logscale2, com, color_gro, ...
    zeroNaNs]        = ParseArgPairs (...
    {'tstr','xtext', 'ytext', 'ticks', 'xstr', 'xlimits', ...
    'logscale','logscale2', 'com','color_gro', ...
    'zeroNaNs'} ...
    ,{ [],[], [], [], [], [], ...
    1, 0, 0, [], ...
    0 }, varargin{ : } );

if logscale2
    logscale = 0;
end
    
colors          = [0.4940, 0.1840, 0.5560; ...
    0.1840,0.5000,0.2000;...
    0.9290 0.6940 0.1250; ...
    0.3010 0.7450 0.9330; ...
    0.8500 0.3250 0.0980; 
    0.369   0.573    0.953; ...
    0.3529    0.8588    0.7020;
    1 0 0 ;
    0 0 1];

if ~isempty(color_gro)
    colors              = colors (color_gro,:);
end


if isempty(cond)
    cond    = true(size(prm));
end

ntypes = size(cond,2);
if isempty(xstr)
    xstr = cell(ntypes,1);
end
if zeroNaNs
    prm(isnan(prm)) = 0;
end

L               = cell(ntypes,1);
g                   = zeros( size( prm ) );

for i = 1 :ntypes
    
    data        = prm( cond(:,i) );
    n           = sum( cond(:,i) & ~isnan(prm) );
    p25         = prctile(data,25);
    p75         = prctile(data,75);
    g( cond(:,i) )    = i;

    if ~com
        med         = nanmedian(data);
    else
        med     = calc_com( data(:), ones(size(data)) );
        sem      = calc_sem(data(:));
    end
    
    if logscale
        Ldata   = log10(data);
        Lmed    = log10(med);
    elseif logscale2
        Ldata   = log2(data);
        Lmed    = log2(med); 
    else
        Ldata   = data;
        Lmed    = med;
    end
   
        
    [cdf,x]     = ecdf( Ldata );
    p(i) = stairs(x,cdf,'LineWidth',2,'Color',colors(i,:) );
    alines(0.5,'y','LineStyle','--','Color','k');
    alines( Lmed,'x','LineStyle','--','Color',colors(i,:));
    
    hold on;
    
    if ~com
        L{i} = sprintf('%s (m=%.3g[%.3g-%.3g];n=%d)', xstr{i}, med, p25, p75, n );
    else
        L{i} = sprintf('%s (m=%.2g+-%.2g; n=%d)', xstr{i}, med, sem,  n );
    end
    
end

legend(p,L);

if isempty(ticks)
    minX    = nanmin(prm);
    maxX    = nanmax(prm);
    ticks   = linspace (minX,maxX,5);
end

if logscale
    xticks(log10(ticks));
    xticklabels(ticks);
    xlim([log10(ticks(1)), log10(ticks(end))]);
elseif logscale2
        xticks(log2(ticks));
    xticklabels(ticks);
    xlim([log2(ticks(1)), log2(ticks(end))]);
end



set( gca, 'tickdir', 'out', 'box', 'off' );
axis square;


if ~isempty(xtext)
    xlabel(xtext);
end

if ~isempty(ytext)
    ylabel(ytext);
end

if ~isempty(xlimits)
    if logscale
        xlim(log10(xlimits));
    elseif logscale2
        xlim(log2(xlimits));
    else
        xlim(xlimits);
    end
end

%statistical test
if ntypes == 1
    title(tstr);
    return;
elseif ntypes == 2 
    [pval,~,stats] = ranksum( prm( cond(:,1) ), prm( cond(:,2) ) ); % Mann-Whitney U-test
    U = stats.ranksum;
else   
    prm(g==0)          = [];
    g(g==0)            = [];
    [ pval, ~, stats ] = kruskalwallis( prm, g, 'off' ); % kruskal wallis
    cc                 = multcompare( stats, 'display', 'off' ); %post hoc multi comparisons
end

title(sprintf('%s; p=%.2g',tstr, pval));

fprintf('\n \n');
if ~isempty(xstr) && ntypes > 2
    mult_comp_text            = cell(size(cc,1),1);
    s                         = '';
    for i = 1 :  size(cc,1) 
        fprintf( '%s vs. %s: %0.2g \n', xstr{ cc( i, 1 ) }, xstr{ cc( i, 2 ) }, cc( i, 6 ) );
    end
end

return;