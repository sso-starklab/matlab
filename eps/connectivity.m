function [act,g1,S] EPS_connectivity 
load ('/probox1/mice/EPS/struct_punits_biphasic_lean_28apr22.mat', '-mat')
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
%%%

isvbip0 = ~isnan(sst.vB) & ~isinf(sst.vB);
isvbip =logical( sum(isvbip0'));
bidx = ~isbip & isvbip' & ~ispos;                  % indeces of Nunits that also have biphasic spikes but not positive spikes
pidx = ~isbip & sst.upol ==0 & ~ispos & ~isvbip';  % indeces of Nunits that also have positive spikes but not biphasic spikes
bpidx = ~isbip & ~ispos & isvbip' & sst.upol ==0;    % indeces of Nunits that also have positive and biphasic spikes
pnidx = ~isbip & ispos & sst.upol ==0;              % indeces of Punits that also have negative and biphasic spikes
pyridx = ~isbip & ~ispos & sst.pyr ==1              ;% Nunits that are PYR
intidx = ~isbip & ~ispos & sst.pyr ==0              ;% Nunits that are INT




% pndix example: mA234_08 2_7, mC41_40 5_6
% we will start with Punits that also have negative and biphasic spikes
sst_con1 = struct_select(sst,pnidx);
%go over all units in sst_con
numunits = length(sst_con1.shankclu);
abs_min = NaN(numunits,1);
abs_min_sample = abs_min;
sample_max = abs_min;
for nn = 1:numunits
    try
    w = sst_con1.mean{nn,1};
    % positive peak sample:
    sample_max(nn) = find(sst_con1.max(:,nn) == sst_con1.extremum(nn));
    % positive peak ch:
    w2 = w(:,sample_max(nn));
    [~,ch_max] = max(abs(w2));
    
    % upsample to fix for bias
    w3                           = fft_upsample( w', 4, 1 );
    % compute
    [ nn2, mm ]                    = size( w3 );                                % number of samples for baseline estimation
    bvals                       = mean( w3( 1 : floor( nn2 / 3 ), : ), 1 );    % baseline
    % remove bias
    mean_amp = bvals';
    wt                              = (w-mean_amp)';
    %find extremum per ch
    numch = size(wt,2);
    nP = NaN(numch,1);
    MnP = nP;
    for i =1:numch
        [ eidx, evals, etype ]       	= local_max( wt(:,i), 'ext' );
        [MnP(i),I] = min(evals);
        nP(i) = eidx(I);
    end
    %remove ch that are punits of bpi
    vB1 = sst_con1.vB(nn,1:numch);
    nanidx = isnan(vB1) | isinf(vB1);
    ridx =  ~nanidx;
    ridx(ch_max) = 1;
    abs_min(nn) = min(MnP(~ridx));
    idxtemp = MnP==abs_min(nn);
    abs_min_sample(nn) =  nP(idxtemp);
    catch
    end
end

% now that we found the punits samples and the negetive spikes min and
% sample we can use threshold to remove the negetive spikes that are
% smaller than 40uW
TH = -40;
kthidx = abs_min<=TH;
delay = abs_min_sample(kthidx) - sample_max(kthidx);
figure;
[x,y,med] = CDF_for_GWN(delay,0);
x=x*0.05; %in ms
ph = plot(x,y);
set( gca, 'tickdir', 'out', 'box', 'off' );
alines(med*0.05,'x')
alines(0.5,'y')
xlabel('Delay between pos and neg peaks [ms]')
ylabel('CDF')
legend(sprintf('n=%d,med = %0.2g',sum(kthidx),med*0.05));


% we will move on to Nunits that also have positive spikes
sst_con2 = struct_select(sst,pidx);
%go over all units in sst_con
numunits = length(sst_con2.shankclu);
abs_min2 = NaN(numunits,1);
abs_min2_sample = abs_min2;
sample_max2 = abs_min2;
for nn = 1:numunits
    try
    w = sst_con2.mean{nn,1};
    % negative peak sample:
    sample_max2(nn) = find(sst_con2.max(:,nn) == sst_con2.extremum(nn));
    % negative peak ch:
    w2 = w(:,sample_max2(nn));
    [~,ch_max] = max(abs(w2));
    
    % upsample to fix for bias
    w3                           = fft_upsample( w', 4, 1 );
    % compute
    [ nn2, mm ]                    = size( w3 );                                % number of samples for baseline estimation
    bvals                       = mean( w3( 1 : floor( nn2 / 3 ), : ), 1 );    % baseline
    % remove bias
    mean_amp = bvals';
    wt                              = (w-mean_amp)';
    %find extremum per ch
    numch = size(wt,2);
    nP = NaN(numch,1);
    MnP = nP;
    for i =1:numch
        [ eidx, evals, etype ]       	= local_max( wt(:,i), 'ext' );
        [MnP(i),I] = max(evals);
        nP(i) = eidx(I);
        
    end
    %remove ch that are punits of bpi
    vB1 = sst_con2.vB(nn,1:numch);
    nanidx = isnan(vB1) | isinf(vB1);
    ridx =  ~nanidx;
    ridx(ch_max) = 1;
    abs_min2(nn) = max(MnP(~ridx));
    idxtemp = MnP==abs_min2(nn);
    abs_min2_sample(nn) =  nP(idxtemp);
    catch
    end
end

% now that we found the Nunits samples and the positive spikes max and
% sample, we can use threshold to remove the positive spikes that are
% smaller than 40uV
TH = 40;
kthidx2 = abs_min2>=TH;
delay2 = -(abs_min2_sample(kthidx2) - sample_max2(kthidx2));
figure;
[x,y,med] = CDF_for_GWN(delay2,0);
x=x*0.05; %in ms
ph = plot(x,y);
set( gca, 'tickdir', 'out', 'box', 'off' );
alines(med*0.05,'x')
alines(0.5,'y')
xlabel('Delay between pos and neg peaks [ms]')
ylabel('CDF')
legend(sprintf('n=%d,med = %0.2g',sum(kthidx2),med*0.05));

totdelay = [delay ;delay2];
figure;
[x,y,med] = CDF_for_GWN(totdelay,0);
x=x*0.05; %in ms
ph = plot(x,y);
set( gca, 'tickdir', 'out', 'box', 'off' );
alines(med*0.05,'x')
alines(0.5,'y')
xlabel('Delay between pos and neg peaks [ms]')
ylabel('CDF')
pv = signrank([abs_min2_sample(kthidx2);sample_max(kthidx)],[sample_max2(kthidx2);abs_min_sample(kthidx)]);
legend(sprintf('n=%d,med = %0.2g',length(totdelay),med*0.05));

IQR = bounds(totdelay,0.5)*0.05;


numUnits = length(isbip);
% kidx = NaN(length (nidx),1);
% lag =  NaN(length (nidx),1);
% wid =  NaN(length (nidx),1);
% biphasic that are  (lag<2) with a soma and wid<1
l=0;
kidx =[];
lag = [];
wid = [];
postidx=[];
preidx=[];
for n = 1 : numUnits % go over all units and treat them as though they are presynaptic
    if ~ispos(n) && ~isbip(n)
        sst_b = NaN(length (sst.lagf {n}) ,1);
        for b = 1:length (sst.lagf {n}) %for each unit (n) look at all the postsynaptic units (b)
            if sst.lagf {n,1}(b)<2 && sst.lagf {n,1}(b)>-0.5  && sst.widf {n,1}(b)<=1 %&& any(sst.pcchU{n,1}(51:53,b) < 0.001)
               afilename = sst.filebase{n};
               idxname = find(ismember(sst.filebase,afilename));
               sst_b(b) = idxname(b); %this is the idx of n inside the whole sst.
               if isbip(sst_b(b))
                  l=l+1;
                  postidx(l) = sst_b(b);
                  preidx(l) = n;
                  lag(l) = sst.lagf {n,1}(b);
                  wid(l) = sst.widf {n,1}(b);
    %               isbipn2(b) = isbipn;
               end
            end
        end 
    end
end


% figure
for i = 1 :  l 
        afilename = sst.filebase{preidx(i)};
        filebase = filebase_lookup(afilename);
            shankclu1                       = sst.shankclu(preidx(i),:);
            shankclu2                       = sst.shankclu(postidx(i),:);
            if i==1 || ~contains(sst.filebase{preidx(i)},sst.filebase{preidx(i-1)})
                spk                             = load_spikes( filebase );
                par                             = LoadXml( filebase );
                SpikesFs                        = par.SampleRate;
                % remove spikes during stimuli
                vals                            = LoadStims( filebase );
                uvals                           = uniteranges( vals( :, 1 : 2 ) );
                ridx                            = inranges( spk.res, uvals );
                spk.clu( ridx )                 = [];
                spk.res( ridx )                 = [];
            end
            % keep only the spike times of the two relevant units
            clunum1                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
            clunum2                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
            st1                             = spk.res( spk.clu == clunum1 );
            st2                             = spk.res( spk.clu == clunum2 );
            % call 
            [ g1(i), g2, act(i), sil, s ]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 1, 'BinSizeMS', 0.05 );% w/ deconvolution and a smaller bin size
            S{i} = s;
            %             suptitle(sprintf('%s, clu1 = %d.%d, clu2 = %d.%d',afilename, shankclu1(1),shankclu1(2),shankclu2(1),shankclu2(2)))
%             pause, 
            fprintf('unit#%d\n',i)
end

%%%%%
% punits that are (lag<0) with a soma and wid<1
n2idx = 1:length(ispos);
k2idx = NaN(length (nidx),1);

for p = 1 : length (n2idx)
    for n = 1:length (sst.lagf {n2idx(p)})
        if sst.lagf {n2idx(p),1}(n)<0 
           a2filename = sst.filebase{n2idx(p)};
           idx2name = find(ismember(sst.filebase,a2filename));
           sst_n = idx2name(n);
           isposn = ispos(sst_n);
           if isposn
              k2idx (p) = sst_n;
              lag(p) = sst.lagf {n2idx(p),1}(n);
              wid(p) = sst.widf {n2idx(p),1}(n);
           end
        end
    end      
end


% find the nunit to bipha with lag less than 1

nunpos_idx = lag<= 1 & wid <= 1;
nun_idx2 = find(nunpos_idx);
nunpos_idx2 = nunpos_idx & ~ispos;
pos_idx = k2idx(nunpos_idx2);
nun_idx2 = find(nunpos_idx2);

shank_post  = sst.shankclu(pos_idx,:);
shank_pre = sst.shankclu(nun_idx2,:);