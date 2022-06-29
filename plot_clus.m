function plot_clus (idx, sst)
for i = 1 : length (idx)
    filename = sst.filebase{idx(i)};
    filebase = filebase_lookup(filename);  
    shankclu = sst.shankclu (idx(i),:);
    try
    figure(i), plot_ss (filebase, shankclu(:,1:2));
    catch
        fprintf ('filebase %s shankclu %d %d not in probox\n', filename, shankclu);
    end
end
return