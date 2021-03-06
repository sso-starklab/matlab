idxNCX = find(units_used.region==1);
idxpostNCX = ismember(units_used.numPost,idxNCX);
idxNCX2 = find(pyr_int.units_used.region==1);
idxpostNCX2 = ismember(pyr_int.units_used.numPost,idxNCX2);
idxCA1 = find(units_used.region==3);
idxpostCA1 = ismember(units_used.numPost,idxCA1);
idxCA12 = find(pyr_int.units_used.region==3);
idxpostCA12 = ismember(pyr_int.units_used.numPost,idxCA12);

figure;
 subplot (1,2,1)
    [x,y,medPI] = CDF_for_GWN(pyr_int.units_used.lagpairs(idxpostNCX2),0);
    ph = plot(x,y);
    hold on
    [xx,yy,medPB]= CDF_for_GWN(units_used.lagpairs(idxpostNCX),0);
     ph = plot(xx,yy);
    legend(sprintf('PYR-INT pairs: n=%d,med = %0.2g',sum(idxpostNCX2),medPI),sprintf('PYR-BIP pairs: n=%d,med = %0.2g',sum(idxpostNCX),medPB));
    set( gca, 'tickdir', 'out', 'box', 'off' );
    alines(1.6,'x','lineStyle','--', 'color','k')
    alines(0.5,'y','lineStyle','--', 'color','k')
    xlabel('lag between pairs [ms]')
    axis square
    title ('nCX')
    xlim([0 5])
 subplot (1,2,2)
    [x,y,medPI] = CDF_for_GWN(pyr_int.units_used.lagpairs(idxpostCA12),0);
    ph = plot(x,y);
    hold on
    [xx,yy,medPB]= CDF_for_GWN(units_used.lagpairs(idxpostCA1),0);
    ph = plot(xx,yy);
    legend(sprintf('PYR-INT pairs: n=%d,med = %0.2g',sum(idxpostCA12),medPI),sprintf('PYR-BIP pairs: n=%d,med = %0.2g',sum(idxpostCA1),medPB));
    set( gca, 'tickdir', 'out', 'box', 'off' );
    alines(1.6,'x','lineStyle','--', 'color','k')
    alines(0.5,'y','lineStyle','--', 'color','k')
    xlabel('lag between pairs [ms]')
     axis square
    title ('CA1')     
    xlim([0 5])



num_pairs= length(s.pval);
for i = 1:num_pairs
    idx1 = contains(sst.filebase,s.filebase{i});
    idx2 = ismember(sst.shankclu,s.shankclu1(i,:),'rows');
    idx3 = ismember(sst.shankclu,s.shankclu2(i,:),'rows');
    numpre = find(idx1&idx2);
    numpost = find(idx1&idx3);
    s.map(i,1:2) = [numpre numpost];
end    
    
    
% shankclu1
idx1_s_bip = isbip(s.map(:,1));
idx1_s_pos = ispos(s.map(:,1));
idx1_s_reg = sst.region(s.map(:,1));
idx1_s_pyr = sst.uType(s.map(:,1));
COM1 = sst.geo_com(s.map(:,1));

% shankclu2
idx2_s_bip = isbip(s.map(:,2));
idx2_s_pos = ispos(s.map(:,2));
idx2_s_reg = sst.region(s.map(:,2));
idx2_s_pyr = sst.uType(s.map(:,2));
COM2 = sst.geo_com(s.map(:,2));


% pairs pyr-bip
pairs_PB = idx1_s_pyr ==1 & idx2_s_bip &s.what==1;
pairs_PB_ncx = pairs_PB & idx1_s_reg==1;
pairs_PB_ca1 = pairs_PB & idx1_s_reg==3;

% pairs int-bip
pairs_IB = idx1_s_pyr ==0 & idx2_s_bip &s.what==1;
pairs_IB_ncx = pairs_IB & idx1_s_reg==1;
pairs_IB_ca1 = pairs_IB & idx1_s_reg==3;

% pairs pyr-int
pairs_PI = idx1_s_pyr ==1 & idx2_s_pyr==0 &s.what==1;
pairs_PI_ncx = pairs_PI & idx1_s_reg==1;
pairs_PI_ca1 = pairs_PI & idx1_s_reg==3;

% pairs pyr-pyr
pairs_PP = idx1_s_pyr ==1 & idx2_s_pyr ==1 &s.what==1;
pairs_PP_ncx = pairs_PP & idx1_s_reg==1;
pairs_PP_ca1 = pairs_PP & idx1_s_reg==3;

% pairs pun-nun
pairs_PuN = idx1_s_pos ==1 & (idx1_s_pyr |~idx1_s_pyr)&s.what==1;
pairs_PuN_ncx = pairs_PuN & idx1_s_reg==1;
pairs_PuN_ca1 = pairs_PuN & idx1_s_reg==3;

% pairs pun-bip
pairs_PuB = idx1_s_pos ==1 & idx2_s_bip &s.what==1;
pairs_PuB_ncx = pairs_PuB & idx1_s_reg==1;
pairs_PuB_ca1 = pairs_PuB & idx1_s_reg==3;

% for inter unit distance we need to check two cases:
% the distance is  sqrt(shankclu(1,1)^2 + COM^2)

pairs = pairs_PB;

% for pairs_PB
CoM_pre = COM1(pairs);
CoM_post = COM2(pairs);

shank1 = s.shankclu1(pairs,1);
shank2 = s.shankclu2(pairs,1);

intershank = abs(shank2-shank1)*200; %in um
innershank = abs(CoM_post-CoM_pre)*20; % in um

dist_tot = sqrt(intershank.^2 + innershank.^2);

%lags
lags = s.lag04(pairs)*1000;
figure; scatter(dist_tot,lags)
lsline
ylabel('Lag [ms]')
xlabel('Distance [um]')
[ cc, pval ] = calc_spearman( dist_tot,lags,1000)    



spk.filebase = 'null';
l=length(s.filebase);
%calc_asg for all pairs
for i = 21223 :  l 
        afilename = s.filebase{i};
        filebase = filebase_lookup(afilename);
        if ~contains(spk.filebase,filebase) || ~contains(s.filebase{i},s.filebase{i-1})
            spk                             = load_spikes( filebase );
            par                             = LoadXml( filebase );
            SpikesFs                        = par.SampleRate;
            % remove spikes during stimuli
            vals                            = LoadStims( filebase );
            if ~isempty(vals)
                uvals                           = uniteranges( vals( :, 1 : 2 ) );
                ridx                            = inranges( spk.res, uvals );
                spk.clu( ridx )                 = [];
                spk.res( ridx )                 = [];
            end
        end
        shankclu1                       = s.shankclu1(i,:);
        shankclu2                       = s.shankclu2(i,:);
        % keep only the spike times of the two relevant units
        clunum1                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
        clunum2                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
        st1                             = spk.res( spk.clu == clunum1 );
        st2                             = spk.res( spk.clu == clunum2 );
        % call 
        [ g1, g2, act, sil, si]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 0, 'BinSizeMS', 0.4, 'halfWidthMS', 30);% w/ deconvolution and a smaller bin size
        s.stg(i) = g1;
        if act
            s.what(i) = 1;
        elseif sil
           s.what(i) = -1;
        end
        % extract the STC as a convolution kernel
        try
        maxlag                      = ( size( si.cch, 1 ) - 1 ) / 2;
        tidx                     = si.g1base( : ) - maxlag;
        kval                     = si.gcch( si.g1base( : ) );
        kval(isnan(kval))=0;
        z                         = zeros( tidx( 1 ) - 1, 1 );
        STC                       = [ z; kval ];                              % [spks/s]
        kP                          = STC * si.dt;                               % [spks/bin]
        ktP                         = ( ( 1 : tidx( end ) )' - 1 ) * si.dt;      % [ms]
        nk                       = length( kP );
        z0                          = zeros( nk - 1, 1 );
        k                       =  [ z0; kP ];                              % zeropad the kernel to be symmetric (but causal)
        t0                       = -( nk - 1 : -1  : 1 )' * si.dt;
        kt                       = [ t0; ktP ];
        % keep
        STCs                  = k;
        STC_time               = kt;
        catch
        STCs                   = 0;
        STC_time               = 0;   
        end
        s.STCs{i} = STCs;
        s.STC_time{i} = STC_time;
        fprintf('pair#%d\n',i)
end

allSTC = s.STCs(pairs_PB);
alltimes = s.STC_time(pairs_PB);
l = length(allSTC);
I=[];
II=[];
sizeSTC=[];
idx1 =[];
idx2=[];
%plot STC
for i = 1:l
     if i==21222
        continue
    end
    sizeSTC(i) = length( allSTC{i});
end
[maxSTC II]= max(sizeSTC);
Bins = alltimes{II};

for i = 1:l
    if i==21222
        continue
    end
    ST = allSTC{i};
    [M I(i)] = max(ST);
    times = alltimes{i};
    maxtime(i) = times(I(i));
    if length(ST) <maxSTC
        addbins = 0.5*(maxSTC-length(ST));
        binstoadd = zeros(addbins,1);
        ST = [binstoadd;ST;binstoadd];
    end
    STCsfull(:,i) = ST;
end

figure;
[idx1 idx2] = sort(I);
imagesc( Bins*1000,1 : l , scale(STCsfull(:,idx2))')
        colormap( myjet )
                axis xy
                xlim([-0.8 5.6])
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel('Time[ms]');
ylabel('Synaptic transmission');
title ('STCs');
axis square

% find lag
lag = maxtime;
figure; histogram(lag,'BinWidth', 0.0002)    
asg = sum(STCsfull);
title ('lag');

s.lag04 = lag;

figure;
[x,y,medPI] = CDF_for_GWN(s.lag04(pairs_PI)*1000,0);
ph = plot(x,y);
hold on
[xx,yy,medPB]= CDF_for_GWN(s.lag04(pairs_PB)*1000,0);
 ph = plot(xx,yy);
legend(sprintf('PYR-INT pairs: n=%d,med = %0.2g',sum(pairs_PI),medPI),sprintf('PYR-BIP pairs: n=%d,med = %0.2g',sum(pairs_PB),medPB));
set( gca, 'tickdir', 'out', 'box', 'off' );
alines(1.6,'x','lineStyle','--', 'color','k')
alines(0.5,'y','lineStyle','--', 'color','k')
xlabel('lag between pairs [ms]')
axis square
xlim([0 5])

figure;
 subplot (1,2,1)
    [x,y,medPI] = CDF_for_GWN(s.lag04(pairs_PI_ncx)*1000,0);
    ph = plot(x,y);
    hold on
    [xx,yy,medPB]= CDF_for_GWN(s.lag04(pairs_PB_ncx)*1000,0);
     ph = plot(xx,yy);
    legend(sprintf('PYR-INT pairs: n=%d,med = %0.2g',sum(pairs_PI_ncx),medPI),sprintf('PYR-BIP pairs: n=%d,med = %0.2g',sum(pairs_PB_ncx),medPB));
    set( gca, 'tickdir', 'out', 'box', 'off' );
    alines(1.6,'x','lineStyle','--', 'color','k')
    alines(0.5,'y','lineStyle','--', 'color','k')
    xlabel('lag between pairs [ms]')
    axis square
    title ('nCX')
    xlim([0 5])
 subplot (1,2,2)
    [x,y,medPI] = CDF_for_GWN(s.lag04(pairs_PI_ca1)*1000,0);
    ph = plot(x,y);
    hold on
    [xx,yy,medPB]= CDF_for_GWN(s.lag04(pairs_PB_ca1)*1000,0);
    ph = plot(xx,yy);
    legend(sprintf('PYR-INT pairs: n=%d,med = %0.2g',sum(pairs_PI_ca1),medPI),sprintf('PYR-BIP pairs: n=%d,med = %0.2g',sum(pairs_PB_ca1),medPB));
    set( gca, 'tickdir', 'out', 'box', 'off' );
    alines(1.6,'x','lineStyle','--', 'color','k')
    alines(0.5,'y','lineStyle','--', 'color','k')
    xlabel('lag between pairs [ms]')
     axis square
    title ('CA1')     
    xlim([0 5])

% F_resampling_test -     
[h,p]             = F_resampling_test(s.lag04(pairs_PI_ncx)*1000 , s.lag04(pairs_PB_ncx)*1000)
[h,p]             = F_resampling_test(s.lag04(pairs_PI_ca1)*1000 , s.lag04(pairs_PB_ca1)*1000)

x = s.lag04(pairs_PI)*1000;
y = s.lag04(pairs_PB)*1000;

figure, boxplot( [ x; y ], [ ones( size( x ) ); 2 * ones( size( y ) ) ], 'notch', 'on' ), ylim( [ -0.5 6 ] )
 

[ ~, pval_ftest  ]      = F_resampling_test( x, y, [], 0, 1000 );           % equal variances
[ ~, pval_kstest ]      = kstest2( x, y );                                  % equal distributions
[ median( x ), bounds( x, 0.5 )' ] % median and IQR
[ median( y ), bounds( y, 0.5 )' ]

% for pairs_PB
CoM_pre = COM1;
CoM_post = COM2;

shank1 = s.shankclu1(:,1);
shank2 = s.shankclu2(:,1);

s = addprobedensity (s);

intershank(s.probeDens==0) = abs(shank2(s.probeDens==0)-shank1(s.probeDens==0))*200; %in um
innershank(s.probeDens==0) = abs(CoM_post(s.probeDens==0)-CoM_pre(s.probeDens==0))*15; % in um
intershank(s.probeDens==1) = abs(shank2(s.probeDens==1)-shank1(s.probeDens==1))*200; %in um
innershank(s.probeDens==1) = abs(CoM_post(s.probeDens==1)-CoM_pre(s.probeDens==1))*20; % in um
intershank(s.probeDens==2) = abs(shank2(s.probeDens==2)-shank1(s.probeDens==2))*200; %in um
innershank(s.probeDens==2) = abs(CoM_post(s.probeDens==2)-CoM_pre(s.probeDens==2))*20; % in um
intershank(s.probeDens==3) = abs(shank2(s.probeDens==3)-shank1(s.probeDens==3))*250; %in um
innershank(s.probeDens==3) = abs(CoM_post(s.probeDens==3)-CoM_pre(s.probeDens==3))*20; % in um
idx_nbs1 = s.probeDens==3 & shank1==1 & shank2==4;
idx_nbs2 = s.probeDens==3 & shank2==1 & shank1==4;
idx_nbs3 = s.probeDens==3 & shank1==2 & shank2==3;
idx_nbs4 = s.probeDens==3 & shank2==2 & shank1==3;
intershank(idx_nbs1|idx_nbs2|idx_nbs3|idx_nbs4) = abs(shank2(idx_nbs1|idx_nbs2|idx_nbs3|idx_nbs4)-shank1(idx_nbs1|idx_nbs2|idx_nbs3|idx_nbs4))*30; %in um

dist_tot = sqrt(intershank.^2 + innershank.^2);
s.dx = dist_tot;

% inhibitory int that are excitatory for biphasic
filebases = unique(s.filebase);
for file = 1:length(filebases)
    flag = ismember(s.filebase, filebases{file});
    sf = struct_select( s, flag );
    presyn = sf.shankclu1;
    [ii,jj,kk]=unique(presyn,'rows','stable');
    f = histc (kk, 1:numel(jj));
    s.freq{file,1} = [f ii ];
end

for file = 1:length(filebases)
    if any(s.freq{1,file}(:,1)>1)
        