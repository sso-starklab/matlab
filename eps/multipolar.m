% 26-apr-22 script for figure 1 
% need to add binomial inlines for bar graphs

% nCx examples:
filename = 'mDS2_07';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [2 6]); % multi
figure, plot_ss (filebase, [2 18]); % single

% CA1 examples:
filename = 'mC41_41';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [3 19]); % multi
figure, plot_ss (filebase, [3 14]); % single


% stats:
load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic_TH_23jun22.mat', '-mat')

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

isvbip0 = ~isnan(sst.vB) & ~isinf(sst.vB);
isvbip =logical( sum(isvbip0'));
bidx = find (~isbip & isvbip' & ~ispos);
pidx = find (~isbip & sst.upol ==0 & ~ispos);
sum(ismember(bidx,pidx));
bidx = ~isbip & isvbip' & ~ispos;                  % indeces of Nunits that also have biphasic spikes but not positive spikes
pidx = ~isbip & sst.upol ==0 & ~ispos & ~isvbip';  % indeces of Nunits that also have positive spikes but not biphasic spikes
bpidx = ~isbip & ~ispos & isvbip' & sst.upol ==0;    % indeces of Nunits that also have positive and biphasic spikes
sum1 = sum(bidx+pidx);
sum2 = length(find(pidx&bidx));
sumtot = sum1-sum2;                 % number of multipolar
sumNunits = sum(~isbip& ~ispos);    % number of Nunits in total
sumSingMod = sumNunits-sumtot;      % number of Single Modal Nunits ("pure Nunits")
sumDualMod = sum(bidx&~pidx&~bpidx)+sum(~bidx&pidx&~bpidx); % number of dual Modal Nunits (Nunits that also have biphasic/positive spikes)
sumTripMod = sum(bpidx);            % number of triple Modal Nunits (Nunits that also have biphasic and positive spikes)

% fig2
sumPyr = ~isbip& ~ispos &sst.pyr;
sumInt = ~isbip& ~ispos &~sst.pyr;

pyr_Nun = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sst.pyr;
int_Nun = ~isbip& ~ispos&~bidx&~pidx&~bpidx&~sst.pyr;
pyr_bidx = bidx&~pidx&~bpidx&sst.pyr;
int_bidx = bidx&~pidx&~bpidx&~sst.pyr;
pyr_pidx = pidx&~bidx&~bpidx&sst.pyr;
int_pidx = pidx&~bidx&~bpidx&~sst.pyr;
pyr_bpidx = bpidx&sst.pyr;
int_bpidx = bpidx&~sst.pyr;

% fig3
dens15 = sst.probeDens == 0;
dens20 = sst.probeDens == 1|sst.probeDens == 2;
totdens15 = sum(dens15&~isbip& ~ispos);
totdens20 = sum(dens20&~isbip& ~ispos);
dens15_Nun = ~isbip& ~ispos&~bidx&~pidx&~bpidx&dens15;
dens20_Nun = ~isbip& ~ispos&~bidx&~pidx&~bpidx&dens20;
dens15_bidx = bidx&~pidx&~bpidx&dens15;
dens20_bidx = bidx&~pidx&~bpidx&dens20;
dens15_pidx = pidx&~bidx&~bpidx&dens15;
dens20_pidx = pidx&~bidx&~bpidx&dens20;
dens15_bpidx = bpidx&dens15;
dens20_bpidx = bpidx&dens20;

% fig4
sumCA1 = ~isbip& ~ispos &sst.region==3;
sumNCX = ~isbip& ~ispos &sst.region==1;

CA1_Nun = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sst.region==3;
NCX_Nun = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sst.region==1;
CA1_bidx = bidx&~pidx&~bpidx&sst.region==3;
NCX_bidx = bidx&~pidx&~bpidx&sst.region==1;
CA1_pidx = pidx&~bidx&~bpidx&sst.region==3;
NCX_pidx = pidx&~bidx&~bpidx&sst.region==1;
CA1_bpidx = bpidx&sst.region==3;
NCX_bpidx = bpidx&sst.region==1;
sumSingMod2 = sum(CA1_Nun)+sum(NCX_Nun);
sumDualMod2 = sum(CA1_bidx)+sum(NCX_bidx)+sum(CA1_pidx)+sum(NCX_pidx);
sumTripMod2 = sum(CA1_bpidx)+sum(NCX_bpidx);

% fig5
[sites, ~]= cellfun(@size,sst.mean);
sum1_8 = ~isbip& ~ispos &sites<9;
sum9_11 = ~isbip& ~ispos &sites>8&sites<12;
sum12_16 = ~isbip& ~ispos &sites>11&sites<17;
sum17_32 = ~isbip& ~ispos &sites>16&sites<33;

Nun1_8 = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sum1_8;
Nun9_11 = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sum9_11;
Nun12_16 = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sum12_16;
Nun17_32 = ~isbip& ~ispos&~bidx&~pidx&~bpidx&sum17_32;

bidx1_8 = bidx&~pidx&~bpidx&sum1_8;
bidx9_11 = bidx&~pidx&~bpidx&sum9_11;
bidx12_16 = bidx&~pidx&~bpidx&sum12_16;
bidx17_32 = bidx&~pidx&~bpidx&sum17_32;
pidx1_8 = pidx&~bidx&~bpidx&sum1_8;
pidx9_11 = pidx&~bidx&~bpidx&sum9_11;
pidx12_16 = pidx&~bidx&~bpidx&sum12_16;
pidx17_32 = pidx&~bidx&~bpidx&sum17_32;
bpidx1_8 = bpidx&sum1_8;
bpidx9_11 = bpidx&sum9_11;
bpidx12_16 = bpidx&sum12_16;
bpidx17_32 = bpidx&sum17_32;

% fig6
% find the span of each unit by going over sst.va then do abs(va) and then
% as when abs(va)>th. we choose TH = 40uV. the span is going to be the idx
% of the last site minus the idx of the first site that passed the TH.
Va = abs(sst.vA);
TH = 40;
idxVA = Va>TH;
[~,B]=max(idxVA,[],2);
idxfind = [];
for i = 1:size(idxVA,1)
    try
    idxfind(i) = find(idxVA(i,:),1,'last');
    catch
    idxfind(i) = NaN;
    end
end
idxfind=idxfind';
span_sites = idxfind-B;
span_sites(span_sites==0)=1;


% now span with FWHM
USF = 32;
fwhm = [];
for i = 1:length(sst.geo_fwhm)
   sw              = sum(  sst.mean{ i }.^2, 2 );                                                % spatial signature of each waveform
    % geometrical dispersion by FWHM
    zu              = fft_upsample( sw( : ), USF, 1 );
    hm              = ( max( zu ) - min( zu ) ) / 2 + min( zu );
    idxFW = zu >= hm;
    idxlast = find(idxFW,1,'last');
    idxfirst = find(idxFW,1,'first');
    span_sitesFW = idxlast-idxfirst+1;
    fwhm (i)           = span_sitesFW / USF;                                        % [sites]    
end

%span = sst.geo_fwhm.*(15*dens15+20*dens20);
% span = span_sites.*(15*dens15+20*dens20);
span = fwhm'.*(15*dens15+20*dens20);
spanuse = span (~isbip&~ispos);
span1 = ~bidx&~pidx&~bpidx&~ispos&~isbip;
span2b = bidx&~pidx&~bpidx;
span2p = pidx&~bidx&~bpidx;
span2 = span2p|span2b;
span3 = bpidx;
spangrp = span1*1 + span2*2 + span3*3;

% for graphics purposes
bidx0 = find(bpidx);

figure
for i = 1 : length( bidx0 )

    filename = sst.filebase{ bidx0( i ) };
    filebase = filebase_lookup(filename);
    clf
    plot_ss( filebase, sst.shankclu( bidx0( i ),:));
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx0 ), sst.filebase{ bidx0( i ) }, sst.shankclu( bidx0( i ) , 1 ) ...
        , sst.shankclu( bidx0( i ) , 2 ), sst.bpi( bidx0( i )  ) );
    title( replacetok( str, '\_', '_' ) )
    pause, 
end

%%
fig1 = figure;
% barwerror([1 2 3],[sumSingMod./sumNunits...
%     sumDualMod./sumNunits  sumTripMod./sumNunits], [], colors_NPB(1,:));
barwerror([1 2],[sumSingMod./sumNunits...
    (sumDualMod+sumTripMod)./sumNunits ], [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Fraction of multimodal units');
ylabel('Fraction');
% set(gca, 'XTickLabel',{'Single Modal' 'Dual Modal' 'Triple Modal'})
set(gca, 'XTickLabel',{'Single Modal' 'Multi Modal'})
ylim ([0 1]);

fig11 = figure;
barwerror([1 2],[sumSingMod (sumDualMod+sumTripMod) ], [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('NUmber of multimodal units');
ylabel('Count');
set(gca, 'XTickLabel',{'Single Modal' 'Multi Modal'})
text(0.8, 8000 ,sprintf('n=%d',sumSingMod));
text(1.9, 8000 ,sprintf('n=%d',(sumDualMod+sumTripMod)));

fig12 = figure;
subplot (1,2,1)
barwerror([1 2],[sum(sst.region==1&(pyr_Nun|int_Nun)) sum(sst.region==1&(bidx|pidx|bpidx)) ], [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Number of multimodal units in nCx');
ylabel('Count');
set(gca, 'XTickLabel',{'Single Modal' 'Multi Modal'})
text(0.8, 2400 ,sprintf('n=%d',sum(sst.region==1&(pyr_Nun|int_Nun))));
text(1.9, 2400 ,sprintf('n=%d',sum(sst.region==1&(bidx|pidx|bpidx))));
subplot (1,2,2)
barwerror([1 2],[sum(sst.region==3&(pyr_Nun|int_Nun)) sum(sst.region==3&(bidx|pidx|bpidx)) ], [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Number of multimodal units in CA1');
ylabel('Count');
set(gca, 'XTickLabel',{'Single Modal' 'Multi Modal'})
text(0.8, 4800 ,sprintf('n=%d',sum(sst.region==3&(pyr_Nun|int_Nun))));
text(1.9, 4800 ,sprintf('n=%d',sum(sst.region==3&(bidx|pidx|bpidx))));

fig2 = figure;
% barwerror([1 2 3 4 5 6],[sum(pyr_Nun)./sumSingMod sum(int_Nun)./sumSingMod ...
%     (sum(pyr_bidx)+sum(pyr_pidx))./sumDualMod (sum(int_bidx)+sum(int_pidx))./sumDualMod ...
%     sum(pyr_bpidx)./sumTripMod sum(int_bpidx)./sumTripMod], [], colors_NPB(1,:));
% set( gca, 'tickdir', 'out', 'box', 'off' );
% title('Fraction of unit types');
% ylabel('Fraction');
% set(gca, 'XTickLabel',{'SingMod pyr' 'SingMod int' 'DualMod pyr' 'DualMod int' 'TripMod pyr' 'TripMod int'})
% ylim ([0 1]);
% alines( sum(sumPyr)./sumNunits, 'y', 'color',colors_IP_light(2,:) , 'linestyle', '--' );
% alines( sum(sumInt)./sumNunits, 'y', 'color',colors_IP_light(1,:) , 'linestyle', '--' );

barwerror([1 2 3 4],[sum(pyr_Nun)./sumSingMod sum(int_Nun)./sumSingMod ...
    (sum(pyr_bidx)+sum(pyr_pidx)+sum(pyr_bpidx))./(sumDualMod+sumTripMod) (sum(int_bidx)+sum(int_pidx)+sum(int_bpidx))./(sumDualMod+sumTripMod)] ...
    , [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Fraction of unit types');
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod pyr' 'SingMod int' 'MultiMod pyr' 'MultiMod int'})
ylim ([0 1]);
alines( sum(sumPyr)./sumNunits, 'y', 'color',colors_IP_light(2,:) , 'linestyle', '--' );
alines( sum(sumInt)./sumNunits, 'y', 'color',colors_IP_light(1,:) , 'linestyle', '--' );

fig21 = figure;
subplot (1,2,1)
barwerror([1 2 3 4],[sum(pyr_Nun&sst.region==1)./sum(NCX_Nun) (sum(pyr_bidx&NCX_bidx)+sum(pyr_pidx&NCX_pidx)+sum(pyr_bpidx&NCX_bpidx))./(sum(NCX_bidx|NCX_pidx|NCX_bpidx))...
    sum(int_Nun&sst.region==1)./sum(NCX_Nun) (sum(int_bidx&NCX_bidx)+sum(int_pidx&NCX_pidx)+sum(int_bpidx&NCX_bpidx))./(sum(NCX_bidx|NCX_pidx|NCX_bpidx))] ...
    , [sem_PYRNcxMM sem_PYRNcxMM sem_INTNcx sem_INTNcxMM], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Fraction of unit types in nCx');
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod pyr' 'SingMod int' 'MultiMod pyr' 'MultiMod int'})
ylim ([0 1]);
alines( sum(sumPyr&NCX_Nun)./sum(NCX_Nun), 'y', 'color',colors_IP_light(2,:) , 'linestyle', '--' );
alines( sum(sumInt&NCX_Nun)./sum(NCX_Nun), 'y', 'color',colors_IP_light(1,:) , 'linestyle', '--' );
subplot (1,2,2)
barwerror([1 2 3 4],[sum(pyr_Nun&sst.region==3)./sum(CA1_Nun) (sum(pyr_bidx&CA1_bidx)+sum(pyr_pidx&CA1_pidx)+sum(pyr_bpidx&CA1_bpidx))./(sum(CA1_bidx|CA1_pidx|CA1_bpidx))...
    sum(int_Nun&sst.region==3)./sum(CA1_Nun) (sum(int_bidx&CA1_bidx)+sum(int_pidx&CA1_pidx)+sum(int_bpidx&CA1_bpidx))./(sum(CA1_bidx|CA1_pidx|CA1_bpidx))] ...
    , [sem_PYRCA1 sem_PYRCA1MM sem_INTCA1 sem_INTCA1MM], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Fraction of unit types in CA1');
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod pyr' 'SingMod int' 'MultiMod pyr' 'MultiMod int'})
ylim ([0 1]);
alines( sum(sumPyr&CA1_Nun)./sum(CA1_Nun), 'y', 'color',colors_IP_light(2,:) , 'linestyle', '--' );
alines( sum(sumInt&CA1_Nun)./sum(CA1_Nun), 'y', 'color',colors_IP_light(1,:) , 'linestyle', '--' );



fig3 = figure;
barwerror([1 2 3 4 5 6],[sum(dens15_Nun)./sumSingMod sum(dens20_Nun)./sumSingMod ...
    (sum(dens15_bidx)+sum(dens15_pidx))./sumDualMod (sum(dens20_bidx)+sum(dens20_pidx))./sumDualMod ...
    sum(dens15_bpidx)./sumTripMod sum(dens20_bpidx)./sumTripMod], [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Modal type as a function of probe density');
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod 15um' 'SingMod 20um' 'DualMod 15um' 'DualMod 20um' 'TripMod 15um' 'TripMod 20um'})
ylim ([0 1]);
alines( totdens15./sumNunits, 'y', 'color',colors_IP_light(2,:) , 'linestyle', '--' );
alines( totdens20./sumNunits, 'y', 'color',colors_IP_light(1,:) , 'linestyle', '--' );

fig4 = figure;
barwerror([1 2 3 4 5 6],[sum(CA1_Nun)./sumSingMod2 sum(NCX_Nun)./sumSingMod2 ...
    (sum(CA1_bidx)+sum(CA1_pidx))./sumDualMod2 (sum(NCX_bidx)+sum(NCX_pidx))./sumDualMod2 ...
    sum(CA1_bpidx)./sumTripMod2 sum(NCX_bpidx)./sumTripMod2], [], colors_NPB(1,:));
set( gca, 'tickdir', 'out', 'box', 'off' );
title('Modal type as a region');
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod CA1' 'SingMod NCX' 'DualMod CA1' 'DualMod NCX' 'TripMod CA1' 'TripMod NCX'})
ylim ([0 1]);
alines( sum(sumCA1)./(sum(sumCA1)+sum(sumNCX)), 'y', 'color', colors_IP_light(2,:) , 'linestyle', '--' );
alines( sum(sumNCX)./(sum(sumCA1)+sum(sumNCX)), 'y', 'color', colors_IP_light(1,:) , 'linestyle', '--' );

fig5 = figure;
title('Modal type as a function of number of sites');
subplot(2,2,1)
barwerror([1 2 3],[sum(Nun1_8)./sum(sum1_8) (sum(bidx1_8)+sum(pidx1_8))./sum(sum1_8) ...
    sum(bpidx1_8)./sum(sum1_8)], [], colors_NPB(1,:));
alines( sumSingMod./sumNunits, 'y', 'color', colors_IP_light(1,:) , 'linestyle', '--' );
alines(sumDualMod./sumNunits, 'y', 'color', colors_NPB(2,:) , 'linestyle', '--' );
alines( sumTripMod./sumNunits, 'y', 'color', colors_NPB(3,:) , 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' );
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod' 'DualMod' 'TripMod'})
title('Fraction of modal type for 1-8 sites');
ylim ([0 1]);
subplot(2,2,2)
barwerror([1 2 3],[sum(Nun9_11)./sum(sum9_11) (sum(bidx9_11)+sum(pidx9_11))./sum(sum9_11) ...
    sum(bpidx9_11)./sum(sum9_11)], [], colors_NPB(1,:));
alines( sumSingMod./sumNunits, 'y', 'color', colors_IP_light(1,:) , 'linestyle', '--' );
alines(sumDualMod./sumNunits, 'y', 'color', colors_NPB(2,:) , 'linestyle', '--' );
alines( sumTripMod./sumNunits, 'y', 'color', colors_NPB(3,:) , 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' );
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod' 'DualMod' 'TripMod'})
    title('Fraction of modal type for 9-11 sites');
ylim ([0 1]);
subplot(2,2,3)
barwerror([1 2 3],[sum(Nun12_16)./sum(sum12_16) (sum(bidx12_16)+sum(pidx12_16))./sum(sum12_16) ...
    sum(bpidx12_16)./sum(sum12_16)], [], colors_NPB(1,:));
alines( sumSingMod./sumNunits, 'y', 'color', colors_IP_light(1,:) , 'linestyle', '--' );
alines(sumDualMod./sumNunits, 'y', 'color', colors_NPB(2,:) , 'linestyle', '--' );
alines( sumTripMod./sumNunits, 'y', 'color', colors_NPB(3,:) , 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' );
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod' 'DualMod' 'TripMod'})
title('Fraction of modal type for 12-16 sites');  
ylim ([0 1]);
subplot(2,2,4)
barwerror([1 2 3],[sum(Nun17_32)./sum(sum17_32) (sum(bidx17_32)+sum(pidx17_32))./sum(sum17_32) ...
    sum(bpidx17_32)./sum(sum17_32)], [], colors_NPB(1,:));
alines( sumSingMod./sumNunits, 'y', 'color', colors_IP_light(1,:) , 'linestyle', '--' );
alines(sumDualMod./sumNunits, 'y', 'color', colors_NPB(2,:) , 'linestyle', '--' );
alines( sumTripMod./sumNunits, 'y', 'color', colors_NPB(3,:) , 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' );
ylabel('Fraction');
set(gca, 'XTickLabel',{'SingMod' 'DualMod' 'TripMod'})
title('Fraction of modal type for 17-32 sites');   
ylim ([0 1]);

fig6 = figure;
% boxplot (span(~isbip&~ispos), spangrp(~isbip&~ispos))
%          tstr                    = 'Unit span [um]';
%          grps                    = {'SingMod' 'DualMod' 'TripMod'};
%          id                      = [ span1,  span2, span3];
%          prm                     = span;
%          color_gro               = [1 2 4 ];
%          ytext                   = '';
%          xtext                   = '';
% make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
%     'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)

         tstr                    = 'Unit span [um]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1,  span2|span3];
         prm                     = span;
         color_gro               = [ 2 4 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
xlim ([0 200])

fig7 = figure;
subplot(2,2,1)
         tstr                    = 'Pyr span in NCx[um]';
         grps                    = {'SingMod' 'DualMod' 'TripMod'};
         id                      = [ span1&sst.pyr&NCX_Nun,  (span2&sst.pyr&NCX_bidx)|(span2&sst.pyr&NCX_pidx),span3&sst.pyr&NCX_bpidx];
         prm                     = span;
         color_gro               = [1 2 4 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,2)
         tstr                    = 'Pyr span in CA1[um]';
         grps                    = {'SingMod' 'DualMod' 'TripMod'};
         id                      = [ span1&sst.pyr&CA1_Nun,  (span2&sst.pyr&CA1_bidx)|(span2&sst.pyr&CA1_pidx), span3&sst.pyr&CA1_bpidx];
         prm                     = span;
         color_gro               = [1 2 4 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,3)
         tstr                    = 'Int span in NCx[um]';
         grps                    = {'SingMod' 'DualMod' 'TripMod'};
         id                      = [ span1&~sst.pyr&NCX_Nun,  (span2&~sst.pyr&NCX_bidx)|(span2&~sst.pyr&NCX_pidx), span3&~sst.pyr&NCX_bpidx];
         prm                     = span;
         color_gro               = [1 2 4 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,4)
         tstr                    = 'Int span in NCx[um]';
         grps                    = {'SingMod' 'DualMod' 'TripMod'};
         id                      = [ span1&~sst.pyr&CA1_Nun,  (span2&~sst.pyr&CA1_bidx)|(span2&~sst.pyr&CA1_pidx), span3&~sst.pyr&CA1_bpidx];
         prm                     = span;
         color_gro               = [1 2 4 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)

fig71 = figure;
subplot(2,2,1)
         tstr                    = 'PYR span in NCx[um]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&sst.pyr&NCX_Nun,  (span2&sst.pyr&NCX_bidx)|(span2&sst.pyr&NCX_pidx)|(span3&sst.pyr&NCX_bpidx)];
         prm                     = span;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
xlim ([0 150])
subplot(2,2,2)
         tstr                    = 'PYR span in CA1[um]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&sst.pyr&CA1_Nun,  (span2&sst.pyr&CA1_bidx)|(span2&sst.pyr&CA1_pidx)|(span3&sst.pyr&CA1_bpidx)];
         prm                     = span;
         color_gro               = [4 1 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
xlim ([0 150])
subplot(2,2,3)
         tstr                    = 'INT span in NCx[um]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&~sst.pyr&NCX_Nun,  (span2&~sst.pyr&NCX_bidx)|(span2&~sst.pyr&NCX_pidx)|(span3&~sst.pyr&NCX_bpidx)];
         prm                     = span;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
xlim ([0 150])
subplot(2,2,4)
         tstr                    = 'INT span in CA1[um]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&~sst.pyr&CA1_Nun,  (span2&~sst.pyr&CA1_bidx)|(span2&~sst.pyr&CA1_pidx)|(span3&~sst.pyr&CA1_bpidx)];
         prm                     = span;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
xlim ([0 150])

fig9 = figure;
subplot(2,2,1)
         tstr                    = 'PYR firing rate in NCx [spk/s]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&sst.pyr&NCX_Nun,  (span2&sst.pyr&NCX_bidx)|(span2&sst.pyr&NCX_pidx)|(span3&sst.pyr&NCX_bpidx)];
         prm                     = sst.frate;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,2)
         tstr                    = 'PYR firing rate in CA1 [spk/s]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&sst.pyr&CA1_Nun,  (span2&sst.pyr&CA1_bidx)|(span2&sst.pyr&CA1_pidx)|(span3&sst.pyr&CA1_bpidx)];
         prm                     = sst.frate;
         color_gro               = [4 1 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,3)
         tstr                    = 'INT firing rate in NCx [spk/s]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&~sst.pyr&NCX_Nun,  (span2&~sst.pyr&NCX_bidx)|(span2&~sst.pyr&NCX_pidx)|(span3&~sst.pyr&NCX_bpidx)];
         prm                     = sst.frate;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,4)
         tstr                    = 'INT firing rate in CA1 [spk/s]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&~sst.pyr&CA1_Nun,  (span2&~sst.pyr&CA1_bidx)|(span2&~sst.pyr&CA1_pidx)|(span3&~sst.pyr&CA1_bpidx)];
         prm                     = sst.frate;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)

fig8 = figure;
subplot(2,2,1)
         tstr                    = 'PYR ACH-COM in NCx [ms]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&sst.pyr&NCX_Nun,  (span2&sst.pyr&NCX_bidx)|(span2&sst.pyr&NCX_pidx)|(span3&sst.pyr&NCX_bpidx)];
         prm                     = sst.ach_com;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,2)
         tstr                    = 'PYR ACH-COM in CA1 [ms]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&sst.pyr&CA1_Nun,  (span2&sst.pyr&CA1_bidx)|(span2&sst.pyr&CA1_pidx)|(span3&sst.pyr&CA1_bpidx)];
         prm                     = sst.ach_com;
         color_gro               = [4 1 ];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,3)
         tstr                    = 'INT ACH-COM in NCx [ms]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&~sst.pyr&NCX_Nun,  (span2&~sst.pyr&NCX_bidx)|(span2&~sst.pyr&NCX_pidx)|(span3&~sst.pyr&NCX_bpidx)];
         prm                     = sst.ach_com;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
subplot(2,2,4)
         tstr                    = 'INT ACH-COM in CA1 [ms]';
         grps                    = {'SingMod' 'MultiMod'};
         id                      = [ span1&~sst.pyr&CA1_Nun,  (span2&~sst.pyr&CA1_bidx)|(span2&~sst.pyr&CA1_pidx)|(span3&~sst.pyr&CA1_bpidx)];
         prm                     = sst.ach_com;
         color_gro               = [4 1];
         ytext                   = '';
         xtext                   = '';
make_cdf_fig_single(prm, id, 'xstr',grps, 'tstr', tstr, ...
    'xtext',xtext, 'color_gro', color_gro, 'ytext', ytext,'logscale',0,'logscale2',0)
