% make_punits_figures_one_hist
%
% see shirly_eps_analysis for preliminary version

% 07-jun-20 ES & SSo

function [ myus, sds, pval, bins_x, hx0, hx1, lh ] = make_punits_figures_one_hist( x, tidx, nbins, colors, str, byprob )

% acquire data
x0                          = x( tidx == 0 );
x1                          = x( tidx == 1 );

% compute mean, SDs, and pval
myus                        = [ nanmean( x0 ) nanmean( x1 ) ];
sds                         = [ nanstd( x0 ) nanstd( x1 ) ];
pval                        = utest( x0, x1);

% generate histogram
bords_x                     = minmax( [ x0; x1 ] );
edges_x                     = linspace( bords_x( 1 ), bords_x( 2 ), nbins + 1 );
hx0                         = histc( x0, edges_x );
hx1                         = histc( x1, edges_x );
hx0( end - 1 )              = hx0( end - 1 ) + hx0( end );
hx0( end )                  = [];
hx1( end - 1 )              = hx1( end - 1 ) + hx1( end );
hx1( end )                  = [];
bins_x                      = ( edges_x( 1 : end - 1 ) + edges_x( 2 : end ) ) / 2;

if byprob
    hx0                     = hx0 / sum( hx0 );
    hx1                     = hx1 / sum( hx1 );
end
       
% plot histogram
newplot
bhx0                        = stairs( bins_x, hx0, 'color', colors( 1, : ), 'linewidth', 1 );
hold on
bhx1                        = stairs( bins_x, hx1, 'color', colors( 2, : ), 'linewidth', 1 );
lh                          = legend( [ bhx0 bhx1 ], { 'Negative', 'Positive' } );
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( '%0.2g/%0.2g; p=%0.2g', myus( 1 ), myus( 2 ), pval ) );
xlabel( str )

return

% EOF
