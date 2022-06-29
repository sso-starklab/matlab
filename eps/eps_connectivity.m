
%multi modal analysis is only for multimodal nunits with pspikes and
%multimodal punits with nspikes


function [units_used,lag,widfwhm,asg,delay,delay2,totdelay,pv,IQR] = eps_connectivity (Nunit_type,unit_type, varargin)
pv = [];
IQR = [];
delay = [];
delay2 = [];
totdelay = [];

nargs                       = nargin;
if nargs < 2 || isempty( Nunit_type ) || isempty(unit_type)
    Nunit_type = 'PYR';
    unit_type = 'BIP';
end

[ graphics, Multimodal...
    ,substructFlag, inhibition]      = ParseArgPairs( ...
    { 'graphics', 'Multimodal'...
    ,substructFlag, inhibition } ...
    , { 0, 0 }...
    , varargin{ : } );



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

sst.uType = zeros(length(sst.pyr),1);
sst.uType(sst.pyr)  = 1;
sst.uType(~sst.pyr)  = 0;
sst.uType(isbip) = 3;
sst.uType(ispos) = 2;

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
allidx = ~isbip & ~ispos;                           %all Nunits
Nunit_types = {'All','PYR','INT'};
unit_types = {'INT','BIP','POS'};

idxNunitType = find(contains(Nunit_types,Nunit_type));
idxunitType = find(contains(unit_types,unit_type));
switch  idxunitType
    case 1
       idxpost =  intidx;
    case 2
        idxpost = isbip;
     case 3
        idxpost = ispos;
end

switch  idxNunitType
    case 1 % all nunits
       idxpre =  2;
    case 2
        idxpre = 1;
     case 3
        idxpre = 0;
end


if Multimodal
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
    if graphics
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
    end

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
    totdelay = [delay ;delay2];
    pv = signrank([abs_min2_sample(kthidx2);sample_max(kthidx)],[sample_max2(kthidx2);abs_min_sample(kthidx)]);
    IQR = bounds(totdelay,0.5)*0.05;

    if graphics
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


        figure;
        [x,y,med] = CDF_for_GWN(totdelay,0);
        x=x*0.05; %in ms
        ph = plot(x,y);
        set( gca, 'tickdir', 'out', 'box', 'off' );
        alines(med*0.05,'x')
        alines(0.5,'y')
        xlabel('Delay between pos and neg peaks [ms]')
        ylabel('CDF')
        legend(sprintf('n=%d,med = %0.2g',length(totdelay),med*0.05));
    end
else

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


    datadir         = '/probox1/mice/EPS/';
    % check mono for each session
    units = struct_select(sst,idxpost);

    usess = unique(units.filebase);
    for i = 1:length (usess)
            asess       = usess{ i };
            uidx        = ismember( sst.filebase, asess );
            ssti            = struct_select( sst, uidx );      
            L1          = load( [ datadir 's2s/' ssti.filebase{ 1 } '.s2s' ], '-mat' );
            L0          = load( [ datadir 'sst/' ssti.filebase{ 1 } '.sst' ], '-mat' );
            gidx        = ismember( L0.sst.shankclu, ssti.shankclu, 'rows' );
            if idxpre == 1
                idxused      =  L0.sst.pyr==1;
            elseif idxpre == 0
                idxused      =  L0.sst.pyr==0; 
            else
                idxused      =  L0.sst.pyr==1 | L0.sst.pyr==0; 
            end
            L1.s2s.shankclu(:,3) = idxused;
            SS(i).mono        = check_mono( L1.s2s ,gidx);
            %frate of PV
            SS(i).frate = L0.sst.frate(gidx);
    end

    % find presynpatic Nunit for each unit
    for i = 1:length(units.shankclu)
       for j = 1:length(SS)
           [a b c] = fileparts(SS(j).mono.filebase);
           if isequal(b, units.filebase{i})
              idx1 =  SS(j).mono.shankclu(:,1) == units.shankclu(i,1);
              idx2 = SS(j).mono.shankclu(:,2) == units.shankclu(i,2);
              idx = idx1&idx2;

              unitNum = find(idx>0);

              PYRpeersidx = SS(j).mono.pairsExc(:,2) == unitNum;
              PYRnum = SS(j).mono.pairsExc(PYRpeersidx,1);
              prePYR = SS(j).mono.shankclu(PYRnum,1:2);
              units.prePyr(i).shankclu = prePYR;
           end
       end
    end

    for post = 1:length(units.shankclu)  
        if ~isempty(units.prePyr(post).shankclu)
            asess       = units.filebase{ post };
            uidx        = ismember( sst.filebase, asess );
            sst_temp = struct_select(sst,uidx);
            ashankclu = units.prePyr(post).shankclu(:,1:2);
            l=0;
            for i = 1: size(ashankclu,1)
               idxtemp = find(ismember(sst_temp.shankclu,ashankclu(i,:),'rows'));
               if idxpre == 1
                    idxused      =  sst_temp.uType(idxtemp)==1;
                elseif idxpre == 0
                    idxused      =   sst_temp.uType(idxtemp)==0;
                else
                    idxused      =  sst_temp.uType(idxtemp)==1 | sst_temp.uType(idxtemp)==0; 
               end
%                if ~idxused
%                 units.prePyr(post).shankclu(i,:) = [];
%                end
               if idxused
                   l=l+1;
                units.prePyr(post).Realshankclu(l,:) = units.prePyr(post).shankclu(i,:);
               end               
            end
        end
    end
    % for n = 1 : numUnits % go over all units and treat them as though they are presynaptic
    % %     if ~ispos(n) && ~isbip(n)
    %      if idxpre(n)
    %         sst_b = NaN(length (sst.lagf {n}) ,1);
    %         for b = 1:length (sst.lagf {n}) %for each unit (n) look at all the postsynaptic units (b)
    %             if sst.lagf {n,1}(b)<2 && sst.lagf {n,1}(b)>-0.5  && sst.widf {n,1}(b)<=1 %&& any(sst.pcchU{n,1}(51:53,b) < 0.001)
    %                afilename = sst.filebase{n};
    %                idxname = find(ismember(sst.filebase,afilename));
    %                sst_b(b) = idxname(b); %this is the idx of n inside the whole sst.
    % %                if isbip(sst_b(b))
    %                 if idxpost(sst_b(b))
    %                   l=l+1;
    %                   postidx(l) = sst_b(b);
    %                   preidx(l) = n;
    %                   lag(l) = sst.lagf {n,1}(b);
    %                   wid(l) = sst.widf {n,1}(b);
    %     %               isbipn2(b) = isbipn;
    %                end
    %             end
    %         end 
    %     end
    % end
    
    % the following loop creates a [m n] cell array S, where m is the size
    % of post-synaptic units, and n is the size of each post-synaptic
    % pre-synaptic units
    % dilute the units structure to keep only units with realPrePYR
    for i = 1:length(units.prePyr)
        idxdilute(i) = isempty(units.prePyr(i).Realshankclu) ;
    end
    idxdilute = logical(idxdilute)';
    keepidx = ~idxdilute;
    units_used = struct_select(units,keepidx);
        
    l = length(units_used.prePyr);
    spk.filebase = 'null';
    % figure
    for i = 1 :  l 
        afilename = units_used.filebase{i};
        filebase = filebase_lookup(afilename);
        ashankclu = units_used.prePyr(i).Realshankclu(:,1:2);
        if ~contains(spk.filebase,filebase) || ~contains(units_used.filebase{i},units_used.filebase{i-1})
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
        for j = 1: size(ashankclu,1)
                shankclu2                       = units_used.shankclu(i,:);
                shankclu1                       = units_used.prePyr(i).Realshankclu(j,1:2);
                % keep only the spike times of the two relevant units
                clunum1                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
                clunum2                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
                st1                             = spk.res( spk.clu == clunum1 );
                st2                             = spk.res( spk.clu == clunum2 );
                % call 
                [ g1, g2, act, sil, s]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 0, 'BinSizeMS', 0.1 );% w/ deconvolution and a smaller bin size
                units_used.prePyr(i).g1(j) = g1;
                units_used.prePyr(i).act(j) = act;
                units_used.prePyr(i).s{j} = s;
                % extract the STC as a convolution kernel
                try
                maxlag                      = ( size( s.cch, 1 ) - 1 ) / 2;
                tidx                     = s.g1base( : ) - maxlag;
                kval                     = s.gcch( s.g1base( : ) );
                kval(isnan(kval))=0;
                z                         = zeros( tidx( 1 ) - 1, 1 );
                STC                       = [ z; kval ];                              % [spks/s]
                kP                          = STC * s.dt;                               % [spks/bin]
                ktP                         = ( ( 1 : tidx( end ) )' - 1 ) * s.dt;      % [ms]
                nk                       = length( kP );
                z0                          = zeros( nk - 1, 1 );
                k                       =  [ z0; kP ];                              % zeropad the kernel to be symmetric (but causal)
                t0                       = -( nk - 1 : -1  : 1 )' * s.dt;
                kt                       = [ t0; ktP ];
                % keep
                STCs                  = k;
                STC_time               = kt;
                catch
                STCs                   = 0;
                STC_time               = 0;   
                end
                units_used.prePyr(i).STCs{j} = STCs;
                units_used.prePyr(i).STC_time{j} = STC_time;
                fprintf('unit#%d\n',i)
           end
    end
    
    m=0;
    for i = 1 :  l 
        ashankclu = units_used.prePyr(i).Realshankclu(:,1:2);
        for j = 1:size(ashankclu,1)
            m=m+1;
            allSTC{m} = units_used.prePyr(i).STCs{j};
           alltimes{m} =  units_used.prePyr(i).STC_time{j};
        end
    end
    
    %plot STC
    for i = 1:m
        sizeSTC(i) = length( allSTC{i});
    end
    [maxSTC II]= max(sizeSTC);
    Bins = alltimes{II};

    for i = 1:m
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
    imagesc( Bins*1000,1 : m , scale(STCsfull(:,idx2))')
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
    figure; histogram(lag,'BinWidth', 0.0001)    
    asg = sum(STCsfull);
    title ('lag');

    % find asg
    figure;
    [~,edges] = histcounts(log10(asg));
    histogram(asg,10.^edges)
    set(gca, 'xscale','log')
	title ('ASG');

    % find width
    Bins_ms = Bins*1000; 
    for i = 1:size(STCsfull,2)
        aa = find(STCsfull(:,i)>0,1,'first'); 
        if ~isempty(aa)
            first_up(i) = aa;
        else
            first_up(i) = 1;
        end  
        [M1(i),I1(i)] = max( STCsfull(first_up(i):end,i)) ;
        [M2(i),I2(i)] = min( STCsfull(first_up(i):end,i)) ;
    end
    hm              = (M1- M2) / 2 + M2;
    fwhm            = sum( STCsfull >= hm ) ;  
    bin_size = 0.1; %ms
    widfwhm = fwhm*bin_size; %ms


end
return
%%%%%
%EoF

%find all INT that are excite BPI:
[units_used_INT_BIP,lag_INT_BIP,widfwhm_INT_BIP,asg_INT_BIP,~,~,~,~,~] = eps_connectivity ('INT','BIP');

%now use this INT to find if there are inhibitory to other units.



%  m=0;
% l = length(units_used_INT_BIP.prePyr);
% 
%     for i = 1 :  l 
%         ashankclu = units_used_INT_BIP.prePyr(i).Realshankclu(:,1:2);
%         for j = 1:size(ashankclu,1)
%             m=m+1;
%             allSTC{m} = units_used_INT_BIP.prePyr(i).STCs{j};
%             alltimes{m} =  units_used_INT_BIP.prePyr(i).STC_time{j};
%             shankpre(m,:) = units_used_INT_BIP.shankclu(i,:);
%             shankpost(m,:) = units_used_INT_BIP.prePyr(i).shankclu(j,:);
%             filenames{m} = units_used_INT_BIP.filebase{i};
%             pair_num(m) = m;
%            
%         end
%     end
% 
% %find max ASG
% [~, I_ASG] = max(asg_INT_BIP);
% shankclu1 = shankpre(I_ASG,:);
% shankclu2 = shankpost(I_ASG,:);
% filename = filenames{I_ASG} ;
% filebase = filebase_lookup(filename,1);
%  spk                             = load_spikes( filebase );
% par                             = LoadXml( filebase );
%  SpikesFs                        = par.SampleRate;
% % remove spikes during stimuli
% vals                            = LoadStims( filebase );
% if ~isempty(vals)
%     uvals                           = uniteranges( vals( :, 1 : 2 ) );
%     ridx                            = inranges( spk.res, uvals );
%     spk.clu( ridx )                 = [];
%     spk.res( ridx )                 = [];
% end
%  clunum1                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
% clunum2                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
% st1                             = spk.res( spk.clu == clunum1 );
% st2                             = spk.res( spk.clu == clunum2 );
% % call 
% [ g1, g2, act, sil, s]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 1, 'BinSizeMS', 0.1 );% w/ deconvolution and a smaller bin size
% 
% 
% 



