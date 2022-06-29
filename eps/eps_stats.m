make_punits_figures( 2);
dbstop at 731
names = {'mV99', 'mF84', 'mC41', 'mK01', 'mF79','mF93', 'mP23','mP101',...
    'mF105','mO251', 'mDL5','mF108', 'mB142','mA234', 'mS234','mDS1', 'mDS2','mP20'};

for i = 1 : length (names)
    idx1            = contains (sst.filebase, names{i});
    sst_temp        = struct_select(sst,idx1);
    num_units       = sum(idx1);
    num_punits      = sum(ispos(idx1));
    num_bipolar     = sum(isbip(idx1));
    num_nunits      = num_units-num_punits-num_bipolar;
    num_pyr         = sum(sst_temp.pyr&~isbip(idx1)&~ispos(idx1));
    num_int         = sum(~sst_temp.pyr&~isbip(idx1)&~ispos(idx1));
    num_sess        = length(unique(sst_temp.filebase));
    
    mice{i,1}.name        = names {i};
    mice{i,1}.num_sess    = num_sess;
    mice{i,1}.num_nunits  = num_nunits;
    mice{i,1}.num_pyr     = num_pyr;
    mice{i,1}.num_int     = num_int;
    mice{i,1}.num_punits  = num_punits;
    mice{i,1}.num_bipolar = num_bipolar;

end