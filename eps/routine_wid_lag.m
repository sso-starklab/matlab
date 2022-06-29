% add s2s.


load ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic.mat', '-mat')
% for the full sst:

nunits1 = length (sst.shankclu);
ilevel = 'B';
wid_f     = cell (nunits1,1);
lag_f     = cell (nunits1,1);
pvalsUpper = cell (nunits1,1);
pvalsLower = cell (nunits1,1);

% ii: mainUnitD
% ii1: mainUnitS
% jj: secUnitD
% jj1: secUnitS


for ii = 1 : nunits1
    filename = sst.filebase{ii};
    filebase = filebase_lookup(filename,1);
    if ii==1 || ~isequal(filename, sst.filebase{ii-1})
        fname               = sprintf( '%s.s2s', filebase  );
        L1                  = load( fname, '-mat' );
        gidx                = check_cluster_quality( filebase, ilevel );
        tmp                 = L1.s2s.t;
        tmp2                = L1.s2s.gcch;
        tmp22               = L1.s2s.pvalsUpper;
        tmp23               = L1.s2s.pvalsLower;
        tmp3                = tmp2(:, gidx, gidx);
        tmp32               = tmp22(:, gidx, gidx);
        tmp33               = tmp23(:, gidx, gidx);
        L1.s2s              = struct_select( L1.s2s, gidx );
        L1.s2s.t            = tmp;
        L1.s2s.gcch         = tmp3;
        L1.s2s.pvalsUpper   = tmp32;
        L1.s2s.pvalsLower   = tmp33;
        ii1 = 0;
    end
    jj = 1;
    ii1 = ii1+1;
    nunitsession = size(L1.s2s.shankclu,1);
    for jj1 = 1 : nunitsession
        [ wid{ii,1}(jj), lag{ii,1}(jj) ] = calc_cch_width( L1.s2s.gcch(:,ii1,jj1), L1.s2s.t);
        [ wid_f{ii,1}(jj), lag_f{ii,1}(jj) ] = calc_cch_width( L1.s2s.gcch(:,ii1,jj1), L1.s2s.t, 1);
        pvalsUpper{ii,1}(:,jj1) = tmp32 (:,ii1,jj1);
        pvalsLower{ii,1}(:,jj1) = tmp33 (:,ii1,jj1);
        jj = jj+1;
    end
end

sst.wid = wid;
sst.lag = lag;
sst.widf = wid_f;
sst.lagf = lag_f;
sst.pcchU = pvalsUpper;
sst.pcchL = pvalsLower;


for i0 = 1 : nunits1
    filename = sst.filebase{i0};
    filebase = filebase_lookup(filename,1);
    if i0==1 || isequal(filename, sst.filebase{i0-1})
        sst.shankcluOrig (i0,2) = j0;
        j0 = j0+1;
    else
        j0=1;
        sst.shankcluOrig (i0,2) = j0;
        j0 = j0+1;
    end
end

save ('/media/shirly/C22865A128659567/mice/EPS/struct_punits_biphasic', 'sst')

%----------------------%

pvalsUpper = cell (nunits1,1);
pvalsLower = cell (nunits1,1);
for iii = 1 : nunits1
    filename = sst.filebase{iii};
    filebase = filebase_lookup(filename,1);
    if iii==1 || ~isequal(filename, sst.filebase{iii-1})
        fname               = sprintf( '%s.s2s', filebase  );
        L1                  = load( fname, '-mat' );
        pvalsU               = L1.s2s.pvalsUpper;
        pvalsL               = L1.s2s.pvalsLower;
    end
    
    pvalsUpper{iii,1}(jj) = pvalsU(:,sst.shankcluOrig(iii),:);
end


L2 = load('/media/shirly/C22865A128659567/mice/mA234/dat/mA234_01/mA234_01.s2s', '-mat')
