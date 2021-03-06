% nCx examples:
filename = 'mS234_03';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [2 42]); % biphasic
filename = 'mP23_06';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [4 11]); % punit
% CA1 examples:
filename = 'mDS2_13';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [1 16]); % biphasic
filename = 'mC41_53';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [3 18]); % punit

load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic_lean_28apr22.mat', '-mat')
is17 = abs(sst.max (17,:))> abs(sst.max(16,:));
is17 = is17';
    for w = 1 : length (sst.filebase)              % if the peak sample is 17 and not 16
        if is17 (w)
            if sst.max(17,w) > 0                   % Punit
                maxw = max(sst.max(17,w));
                sst.extremum (w) = maxw;
            else                                   % Nunit
                minw = min(sst.max(17,w));
                sst.extremum (w) = minw;
            end
        end
    end 
isbip           = ~isnan(sst.bpi) & ~isinf( sst.bpi );
ispos           = sst.extremum > 0 & ~isbip;

colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Biphasic
colors_IP_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR

sumNCX = sst.region==1;
sumCA1 = sst.region==3;
dens15 = sst.probeDens == 0;
dens20 = sst.probeDens == 1|sst.probeDens == 2;
span = sst.geo_fwhm.*(15*dens15+20*dens20);

condNCX = [~ispos & ~isbip & sumNCX, ispos & sumNCX, isbip & sumNCX];
sumNCX  = [sum(~ispos & ~isbip & sumNCX), sum(ispos & sumNCX), sum(isbip & sumNCX)];

condCA1 = [~ispos & ~isbip & sumCA1, ispos & sumCA1, isbip & sumCA1];

   slice_names                 = {'Nunits', 'Punits', 'Bipolar'};
for spi                     = 1 : 2
    subplot( 1, 2, spi )
    switch spi
        case 1
            % pie chart for the fraction of the punits out of the total units
            sums            = [ sum( ~ispos & ~isbip ) sum( ispos ) sum( isbip )];
        case 2
            % pie chart for the fraction of the punits out of the total units
            sums            = [ sum( sst.nspks( ~ispos & ~isbip ) ) sum( sst.nspks( ispos ) ) sum( sst.nspks( isbip ) )];
    end
    ph                      = pie( sums, slice_names );
    rez{ 2, 2 }.slice_names( spi, : )   = slice_names;
    rez{ 2, 2 }.sums( spi, : )          = sums;
    for i                   = 1 : 3
        set( ph( 2 * i - 1 ), 'FaceColor', colors_NPB( i, : ), 'EdgeColor', colors_NPB( i, : ) )
        set( ph( 2 * i ), 'String', sprintf( '%s (%d)', slice_names{ i }, sums( i ) ) )
    end
    title( sprintf( '%0.2g%%', round( sums( 2 ) / sum( sums ) * 1000 ) / 10 ) )
end


xstr = {'Nunits', 'Punits', 'Biphasic'};
% prm - values
% cond - matrix [n m] n length of prm
%                     m number of group
% tstr - title
% xtext - x-axes legend
% ytext - y-axes legend
% ticks - custum c ticks
% xstr - group names in legend
% com - calculation of center od mass instead of median

% nCx
figure,
subplot(2,2,1)
make_cdf_fig_single( 1 ./ sst.fmax * 1000 , condNCX, 'tstr', 'nCx Width [ms]', 'xstr', xstr, 'logscale', 0, 'xlimits', [0 1.5])
subplot(2,2,2)
make_cdf_fig_single( span, condNCX, 'tstr', 'nCx FWHM [um]', 'xstr', xstr, 'logscale', 0, 'xlimits', [0 110])
subplot(2,2,3)
make_cdf_fig_single( sst.frate, condNCX, 'tstr', 'nCx firing rate [spk/s]', 'xstr', xstr, 'logscale', 0)
set(gca,'XScale','log')
set(gca,'XTick',[0.01 0.1 1 10])
subplot(2,2,4)
make_cdf_fig_single( sst.ach_com, condNCX, 'tstr', 'nCx ACH-COM [ms]', 'xstr', xstr, 'logscale', 0, 'xlimits', [10 40])
% CA1
figure,
subplot(2,2,1)
make_cdf_fig_single( 1 ./ sst.fmax * 1000 , condCA1, 'tstr', 'CA1 Width [ms]', 'xstr', xstr, 'logscale', 0, 'xlimits', [0 1.5])
subplot(2,2,2)
make_cdf_fig_single( span, condCA1, 'tstr', 'CA1 FWHM [um]', 'xstr', xstr, 'logscale', 0, 'xlimits', [0 110])
subplot(2,2,3)
make_cdf_fig_single( sst.frate, condCA1, 'tstr', 'CA1 firing rate [spk/s]', 'xstr', xstr, 'logscale', 0)
set(gca,'XScale','log')
set(gca,'XTick',[0.01 0.1 1 10])
 subplot(2,2,4)
make_cdf_fig_single( sst.ach_com, condCA1, 'tstr', 'CA1 ACH-COM [ms]', 'xstr', xstr, 'logscale', 0, 'xlimits', [10 40])


% to include multimodals in the pies: run index lines in 'multipolar.mat', then this script:
colors_NMPB                 = [ 128 130 132; 145 47 114; 0 163 126; 42 101 176] / 255;
slice_names                 = {'Nunits', 'Multi modal', 'Punits', 'Biphasic'};
for rgi                     = 1 : 2
    subplot( 1, 2, rgi )
    switch rgi
        case 1
            % pie chart for the fractions out of the total units in nCX
            sums            = [ sum( NCX_Nun ) sum( NCX_bidx| NCX_pidx| NCX_bpidx) sum( ispos&sst.region==1 ) sum( isbip&sst.region==1 )];
        case 2
            % pie chart for the fractions out of the total units in CA1
            sums            = [ sum( CA1_Nun ) sum( CA1_bidx| CA1_pidx| CA1_bpidx ) sum( ispos&sst.region==3 ) sum( isbip&sst.region==3 )];
    end
    ph                      = pie( sums, slice_names );
    rez{ 2, 2 }.slice_names( rgi, : )   = slice_names;
    rez{ 2, 2 }.sums( rgi, : )          = sums;
    for i                   = 1 : 4
        set( ph( 2 * i - 1 ), 'FaceColor', colors_NMPB( i, : ), 'EdgeColor', colors_NMPB( i, : ) )
        set( ph( 2 * i ), 'String', sprintf( '%s (%d)', slice_names{ i }, sums( i ) ) )
    end
end
xstr = {'Nunits', 'Multi modal', 'Punits', 'Biphasic'};