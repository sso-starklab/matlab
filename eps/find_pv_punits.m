function s_light = find_pv_punits (sst1, pidx)
j=1;
for i = 1: length (sst1.filebase)
    if isequal (sst1.opsinType{i}, 'camkii::chr2') && sst1.act (i)==1 && pidx(i)==1
        disp (sst1.filebase(i))
        disp (sst1.shankclu(i,:))
        s_light.filebase{j,:} = sst1.filebase{i};
        s_light.shankclu(j,:) = sst1.shankclu(i,:);
        try
        filebase = filebase_lookup(s_light.filebase{j});
        figure, plot_ss(filebase, s_light.shankclu(j,1:2));
        catch
        end
        j=j+1;
    end
end
return

% EOF
