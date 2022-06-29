% make_punits_figures_add_iso_ellipses

% see shirly_eps_analysis for preliminary version

% 25-jun-20 ES & SSo

function h = make_punits_figures_add_iso_ellipses( param1, param2, labels, colors )

% constants
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
