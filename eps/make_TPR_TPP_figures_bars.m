% make_hifi_figures_bars        internal callback from make_hifi_figures
%
% call                          make_jpn_figures_bars( param, id, stimod, xstr, ALPHA )
%
%                               log_prm     logical parameter
%                               id          group logical vectors (a column per each group)
%                               xstr        name of each group
%                               tstr        name string
%                               alpha   
% 
% calls                         binomial_inlines, inline_stats, gtest (stats)
%                               make_cl_on_bar (graph)

% 31-dec-19 HS and ES

% revisions
% 01-jan-20 added statistical testing (pairwise g-test; comparison to chance level)
%           added line for chance level
% 04-Feb-20 adapted for JPN project
% 30-Dec-20 changed chi square test to g test

function make_TPR_TPP_figures_bars( log_prm, id, xstr, tstr, alpha, varargin )


[ color_gro, ytext, ylimits]        = ParseArgPairs (...
    {'color_gro', 'ytext', 'ylimits'} ...
    ,{[], [], [] }, varargin{ : } );

colors                  = [106,27,154;  46,125,50; 173, 20, 87;  0 0.7 0; 0, 0.75, 0.75]/255;

if ~isempty(color_gro)
    colors              = colors (color_gro,:);
end

% inlines
[ ~, ~, bino_se_norm ]  = binomial_inlines;
[ binp, ~, ~, chi2xy ]  = inline_stats;

% preps
[rows1, ~]                   = size(log_prm);
[rows2, n_grp]           = size(id);
counts                  = NaN( n_grp, 2 );

for i = 1 : n_grp
    
    if ( (rows1 < 2) && (rows2 < 2) )
            counts(i,:)     = [log_prm(i), id(i)];
    elseif ( (rows1 > 2) && (rows2 < 2) )
            counts(i,:)     = [sum( log_prm(:,i) ), id(i)];
    else
        counts( i, : )      = [ sum( log_prm & id( :,i ) ), sum( id(:,i) ) ];
    end
 
end


% computations
f                       = counts( :, 1 ) ./ counts( :, 2 );
se                      = bino_se_norm( counts( :, 1 ), counts( :, 2 ) );
pvals                   = binp( counts( :, 1 ), counts( :, 2 ), alpha );

% % pair-wise comparisons
C                       = combnk(1 : n_grp,2);
combs                   =  size(C,1);
cc                      = NaN (combs, 3);
cc (:,1:2)              = C;

for i                   = 1 : combs
    j                   = cc( i, 1 );
    k                   = cc( i, 2 );
   % cc( i, 3 )          = chi2xy( counts( j, 1 ), counts( j, 2 ), counts( k, 1 ), counts( k, 2 ) );
    cc( i, 3 )          = gtest( [counts( j, 1 ), counts( j, 2 ); counts( k, 1 ), counts( k, 2 ) ]);

end

% plot
hold on
for i                   = 1 : n_grp
    bh                  = bar( i, counts( i, 1 ) ./ counts( i, 2 ) );
    set( bh, 'EdgeColor', colors( i, : ), 'FaceColor', colors( i, : ) );
    str                 = sprintf( '%d/%d', counts( i, 1 ), counts( i, 2 ) );
    text( i, 1.05, str );
    make_cl_on_bar( i, f( i ) + se( i ), f( i ) - se( i ), [ 0 0 0 ], 0.25, 2 );
    xstrp{ i }          = sprintf( '%s (%0.2g)', xstr{ i }, pvals( i ) );
end
set( gca, 'tickdir', 'out', 'box', 'off' );
title( tstr )
set( gca, 'xtick', 1 : n_grp )
set( gca, 'XTickLabel', xstrp )
xlim( [ 0.5, n_grp+0.5 ] )
ylim( [ 0, 1.1 ] )
alines( alpha, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );

if ~isempty(ylimits)
    ylim(ylimits);
end

if ~isempty(ytext)
    ylabel(ytext);
else
    ylabel('Fraction');
end

for i = 1: combs
    str                 = sprintf( '%s vs. %s: %0.2g', xstr{ cc( i, 1 ) }, xstr{ cc( i, 2 ) }, cc( i, 3 ) );
    text( i, 0.8, str );    
end


return

% EOF