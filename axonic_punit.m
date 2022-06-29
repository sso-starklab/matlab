% the goal is to find a unequivocally axonic punit: if a punit is ~synchronized (synchrony with a small delay) with a
% nunit, and both have similar (e.g. excitatory) connections with a third unit, and the
% connectivity between the nunits can be predicted from the connectivity between {nunit1
% and punit}, and {punit and nunit2}. (convolution of the STCs or multiplication of the ASGs)
function statsP = axonic_punit (idx, sst)
for i = 1 : length (idx)
    
    try
    filename = sst.filebase{idx(i)};
    
    statsP.filename{i}= filename;
    
    filebase = filebase_lookup(filename); 
    
    shankclu = sst.shankclu (idx(i),:);
    statsP.Pshankclu(i,:)= shankclu;
    
    mono = check_mono(filebase);
    catch
        continue
    end
    statsP.monoshankclu{i} = mono.shankclu;
    statsP.monoexc{i} = mono.pairsExc;
    statsP.monosync{i} = mono.pairsSync;
    idx0 = mono.shankclu==shankclu; 
    idx1 = idx0(:,1)&idx0(:,2);
    pnum = find(idx1);
    pexc = mono.pairsExc (:,1)== pnum;
    % now find the units that were excited by the punit
    excited = mono.pairsExc(pexc,2);
    % now find other units that excite the excited units
    aa = size(mono.pairsExc,1);
    bb = size(mono.pairsSync,1);
    cc = max(aa,bb);
    
    spp= NaN(cc,cc);
    spp2= NaN(cc,cc);
    spp3= NaN(cc,cc);

    for j = 1 : length (excited)
        aexcited = excited (j);
        aidx = mono.pairsExc( :,2)== aexcited;
        exciting = mono.pairsExc( aidx,1);
        sp = [];
        for l=1:length(exciting)
         syncidx = mono.pairsSync(:,1) == exciting(l);
         syncnum = mono.pairsSync(syncidx,2);
         for k =1 : length(syncnum)
             if ismember(syncnum(k),pnum) %if this happens then the exciting unit is in sync with the punit!!!!
               sp(l) =  exciting(l);
               continue
             else
                  sp(l) = NaN;
             end
         end
        end
        sp2 = [];
        for l=1:length(exciting)
         syncidx = mono.pairsSync(:,2) == exciting(l);
         syncnum = mono.pairsSync(syncidx,1);
         for k =1 : length(syncnum)
             if ismember(syncnum(k),pnum) %if this happens then the exciting unit is in sync with the punit!!!!
               sp2(l) =  exciting(l);
               continue
             else
                  sp2(l) = NaN;
             end
         end
        end
               
         sp3 = [];
        for l=1:length(exciting)
         excidx = mono.pairsExc(:,1) == exciting(l);
         excnum = mono.pairsExc(excidx,2);
         for k =1 : length(excnum)
             if ismember(excnum(k),pnum) %if this happens then the exciting unit is exciting the punit!!!!
               sp3(l) =  exciting(l);
               continue
             else
                  sp3(l) = NaN;
             end
         end
        end
        spp(j,1:length(sp))= sp;
        spp2(j,1:length(sp2))= sp2;
        spp3(j,1:length(sp3))= sp3;

    end
    nanidx = ~isnan(spp);
    numsync = sum(sum(nanidx));
    statsP.spp {i}= spp;
    statsP.numsync(i) = numsync;
    
    nanidx2 = ~isnan(spp2);
    numsync2 = sum(sum(nanidx2));
    statsP.spp2 {i}= spp2;
    statsP.numsync2(i) = numsync2;
    
    nanidx3 = ~isnan(spp3);
    numexc = sum(sum(nanidx3));
    statsP.spp3 {i}= spp3;
    statsP.numexc(i) = numexc;
end
return