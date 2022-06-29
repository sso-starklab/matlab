
load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic.mat', '-mat')

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

pidx = find( ispos );
kidx = NaN(length (pidx),1);
% Punits that are "synchronized" (lag<1) with a soma and wid<1
for p = 1 : length (pidx)
    for n = 1:length (sst.lag {pidx(p)})
        if sst.lag {pidx(p),1}(n)<1
            kidx (p) = n;
        end
    end      
end

bidx = find( isbip );
kidx = NaN(length (bidx),1);
% biphasic that are  (lag<2) with a soma and wid<1
for b = 1 : length (bidx)
    for n = 1:length (sst.lag {bidx(b)})
        if any(sst.lag {bidx(b),1}(n)<1)
            kidx (b) = n;
        end
    end      
end

nidx = 1:length(isbip);
kidx = NaN(length (nidx),1);
lagRuff =  NaN(length (nidx),1);
widRuff=  NaN(length (nidx),1);
% biphasic that are  (lag<2) with a soma and wid<1
for b = 1 : length (nidx)
    for n = 1:length (sst.lagf {nidx(b)})
        if sst.lagf {nidx(b),1}(n)<2
           afilename = sst.filebase{nidx(b)};
           idxname = find(ismember(sst.filebase,afilename));
           sst_n = idxname(n);
           isbipn = isbip(sst_n);
           if isbip(sst_n)
              kidx (b) = sst_n;
              lagRuff(b) = sst.lagf {nidx(b),1}(n);
              widRuff(b) = sst.widf {nidx(b),1}(n);
%               isbipn2(b) = isbipn;
           end
        end
    end      
end

% find the nunit to bipha with lag less than 1
pupikidx = lagRuff== 1 & widRuff <= 1;
pupik = kidx(pupikidx);
pupik2 = find(pupikidx);

shank_post  = sst.shankclu(pupik,:);
shank_pre = sst.shankclu(pupik2,:);
isbip(pupik);
%exmpale
i= 4;
afilename = sst.filebase{pupik2(i)};
filebase = filebase_lookup(afilename);
figure; plot_ss(filebase,shank_pre(i,:))
figure; plot_ss(filebase,shank_post(i,:))

figure
count = 0;
for i = 1 : length( shank_post )
        afilename = sst.filebase{pupik2(i)};
        filebase = filebase_lookup(afilename);
        if any(sst.pcchU{pupik2(i),1}(51:52,sst.shankcluOrig(pupik(i),2))<0.05)
%             figure; plot_ss(filebase,shank_pre(i,:))
%             figure; plot_ss(filebase,shank_post(i,:))
%             str = sprintf( '%d/%d, %s, %d.%d', i,length( shank_post ), sst.filebase{ pupik2(i)}, sst.shankclu( pupik2(i) , 1 ), sst.shankclu( pupik2(i) , 2 )) ;
%             title( replacetok( str, '\_', '_' ) )
%             pause, 
            count = count+1;
        end
end


ridx = zeros (length(bidx),1);
for b = 1 : length (bidx)
    if max(sst.max(:,bidx(b)))< 150
        ridx (b) = 1;
    end
end

bidx (find(ridx)) = [];
bidx1 = bidx (sst.bpi(bidx)<-1 &sst.bpi(bidx)>-2 );

bidx1 = find(isbip & eidx2); % excited bipolars
bidx = isbip & eidx2;
bidx1 = find(s.asg1(bidx)>0.5);

bidx1 = find(sst.upol ==0);

isvbip0 = ~isnan(sst.vB) & ~isinf(sst.vB);
isvbip =logical( sum(isvbip0'));
bidx1 = find (~isbip & isvbip' & ~ispos);
bidx2 = find (~isbip & sst.upol ==0 & ~ispos);
sum(ismember(bidx1,bidx2))
bidx1 = ~isbip & isvbip' & ~ispos;  % indeces of Nunits that also have biphasic spikes
bidx2 = ~isbip & sst.upol ==0 & ~ispos; % indeces of Nunits that also have positive spikes
sum1 = sum(bidx1+bidx2);
sum2 = length(find(bidx2&bidx1));
sumtot = sum1-sum2;

bidx1 = find (~isbip & ~isvbip' & ~ispos & sst.upol ==0);

figure
for i = 1 : length( bidx1 )

    filename = sst.filebase{ bidx1( i ) };
    filebase = filebase_lookup(filename);
    clf
    plot_ss( filebase, sst.shankclu( bidx1( i ),:));
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx1 ), sst.filebase{ bidx1( i ) }, sst.shankclu( bidx1( i ) , 1 ) ...
        , sst.shankclu( bidx1( i ) , 2 ), sst.bpi( bidx1( i )  ) );
    title( replacetok( str, '\_', '_' ) )
    pause, 
end

figure
for i = 1 : length( bidx1 )

    filename = sst_p.filebase{ bidx1( i ) };
    filebase = filebase_lookup(filename);
    clf
    plot_ss( filebase, sst_p.shankclu( bidx1( i ),:));
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx1 ), sst_p.filebase{ bidx1( i ) }, sst_p.shankclu( bidx1( i ) , 1 ) ...
        , sst_p.shankclu( bidx1( i ) , 2 ), sst_p.bpi( bidx1( i )  ) );
    title( replacetok( str, '\_', '_' ) )
    pause, 
end

    
for i = 22 : length( bidx1 )
    filename = sst.filebase{ bidx1( i ) };
    filebase = filebase_lookup(filename);
    figure
    plot_ss( filebase, sst.shankclu( bidx1( i ),:));
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx1 ), sst.filebase{ bidx1( i ) }, sst.shankclu( bidx1( i ) , 1 ) ...
        , sst.shankclu( bidx1( i ) , 2 ), sst.bpi( bidx1( i )  ) );
    figure
    plot_ss( filebase, sst_abvTH.shankclu( bidx1( i ),:));
    str = sprintf( '%d/%d, %s, %d.%d, %0.3g', i, length( bidx1 ), sst_abvTH.filebase{ bidx1( i ) }, sst_abvTH.shankclu( bidx1( i ) , 1 ) ...
        , sst_abvTH.shankclu( bidx1( i ) , 2 ), sst_abvTH.bpi( bidx1( i )  ) );
    title( replacetok( str, '\_', '_' ) )
    pause, 
end