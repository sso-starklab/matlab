% re plot the data for Cortical Sensing

pulke.rule1 = [0.756757	0.910448	0.828571];%succ rate
pulke.testnum1 = [74	134	70];% testing trials
pulke.sem1 = [0.05121963117	0.03831829117	0.04435756439]; %SEM
for i= 1:length(pulke.rule1)
    pv(i) = myBinomTest(pulke.rule1(i)*pulke.testnum1(i), pulke.testnum1(i), 0.5,'one');
end
pulke.pval1 = pv;
nn1=i;

pulke.rule2 = [0.4642857143	0.3866666667	0.5208333333	0.6011904762	0.585	0.6319444444	0.5897435897	0.6805555556	0.6607142857	0.7272727273	0.6666666667	0.7058823529];
pulke.testnum2 = [224	300	96	168	200	144	156	144	168	132	144	204];
pulke.sem2 = [0.04502990482	0.02935608075	0.04747815303	0.0403975427	0.04636576904	0.034701873	0.03461936666	0.03815138068	0.04555678078	0.03689011716	0.03255614334	0.02412135777];
for i= 1:length(pulke.rule2)
    pv2(i) = myBinomTest(pulke.rule2(i)*pulke.testnum2(i), pulke.testnum2(i), 0.5,'one');
end
pulke.pval2 = pv2;
nn2=i;

fig1 = figure;
errorbar(1:nn1,pulke.rule1,pulke.sem1)
hold on
errorbar(nn1+1:nn1+nn2,pulke.rule2,pulke.sem2)

pv_tot = [pv pv2 ];
for i = 1:length(pv_tot)
if pv_tot(i) <=0.01
plot(i,1,'*k')
end
end
ylim ([0 1])
alines([nn1+0.5,nn1+nn2+0.5],'x','LineStyle','--')
alines([0.5],'y','LineStyle','--')
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel('Session number')
ylabel('Success rate')
title('mA271')

figname = '/amir1/cortical_sensing/figs/Pulke_visual';
save_print(figname)


%pulke cortical
pulke.rule6 = [0.5625	0.5333];
pulke.testnum6 = [80	120];
pulke.sem6 = [0.07544510777	0.03956837836];
for i= 1:length(pulke.rule6)
    pv6(i) = myBinomTest(pulke.rule6(i)*pulke.testnum6(i), pulke.testnum6(i), 0.5,'one');
end
pulke.pval6 = pv;
nn6=i;

pulke.rule7 = [0.5157	0.5428	0.5667	0.6714	0.6714];
pulke.testnum7 = [96	140	150	70	70];
pulke.sem7 = [0.03239417719	0.04285714286	0.04542567626	0.06801360408	0.05654448613];
for i= 1:length(pulke.rule7)
    pv7(i) = myBinomTest(pulke.rule7(i)*pulke.testnum7(i), pulke.testnum7(i), 0.5,'one');
end
pulke.pval7 = pv7;
nn7=i;


fig2 = figure;
errorbar(1:nn6,pulke.rule6,pulke.sem6)
hold on
errorbar(nn6+1:nn6+nn7,pulke.rule7,pulke.sem7)
pv_tot = [pv6 pv7];
for i = 1:length(pv_tot)
if pv_tot(i) <=0.01
plot(i,1,'*k')
end
end
alines([2.5],'x','LineStyle','--')
alines([0.5],'y','LineStyle','--')
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel('Session number')
ylabel('Success rate')
title('mC41 - visual rules')
xlim([0 8])

axis square

figname = '/amir1/cortical_sensing/figs/Pulke_cortical';
save_print(figname)

%%
%now Rachel
 %rachel visual:
rachel.rule1 = [0.7	0.6833333333	0.56	0.56	0.5666666667	0.575	0.7428571429	0.7333333333];
rachel.testnum1 = [40	60	96	50	120	40	70	90];
rachel.sem1 = [0	0.08333333333	0.03651483717	0.05099019514	0.04322831096	0.1030776406	0.0685118789	0.02886751346];
for i= 1:length(rachel.rule1)
    pv(i) = myBinomTest(rachel.rule1(i)*rachel.testnum1(i), rachel.testnum1(i), 0.5,'one');
end
rachel.pval1 = pv;
nn1=i;

rachel.rule2 = [0.6384615385	0.67	0.725	0.7142857143];
rachel.testnum2 = [130	100	80	70];
rachel.sem2 = [0.03825800849	0.04955356249	0.03133915853	0.02608202655];
for i= 1:length(rachel.rule2)
    pv2(i) = myBinomTest(rachel.rule2(i)*rachel.testnum2(i), rachel.testnum2(i), 0.5,'one');
end
rachel.pval2 = pv2;
nn2=i;

rachel.rule3 = [0.775	0.7384615385	0.7857142857];
rachel.testnum3 = [80	130	70];
rachel.sem3 = [0.025	0.03496969666	0.07046975518];
for i= 1:length(rachel.rule3)
    pv3(i) = myBinomTest(rachel.rule3(i)*rachel.testnum3(i), rachel.testnum3(i), 0.5,'one');
end
rachel.pval3 = pv3;
nn3=i;

fig1 = figure;
errorbar(1:nn1,rachel.rule1,rachel.sem1)
hold on
errorbar(nn1+1:nn1+nn2,rachel.rule2,rachel.sem2)
errorbar(nn1+nn2+1:nn1+nn2+nn3,rachel.rule3,rachel.sem3)
pv_tot = [pv pv2 pv3];
for i = 1:length(pv_tot)
if pv_tot(i) <=0.01
plot(i,1,'*k')
end
end
alines([8.5 12.5],'x','LineStyle','--')
alines([0.5],'y','LineStyle','--')
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel('Session number')
ylabel('Success rate')
title('mA154 - visual rules')

figname = '/amir1/cortical_sensing/figs/Rachel_visual';
save_print(figname)

