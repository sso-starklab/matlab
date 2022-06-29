function filebase = filebase_lookup(filename)
aa = find(filename(:) == '_');
mouse_name = filename(1:aa-1);
if isunix
    if size(mouse_name,2) < 5 
        filebase = sprintf('/media/shirly/C22865A128659567/mice/%s/dat/%s/%s',mouse_name,filename,filename);
        return
    else
        filebase = sprintf('/media/shirly/C22865A128659567/mice/%s/dat/%s/%s',mouse_name,filename,filename); 
        return
    end
elseif ispc
    if size(mouse_name,2) < 5 
        filebase = sprintf('G:\\mice\\%s\\dat\\%s\\%s',mouse_name,filename,filename);
        return
    else
        filebase = sprintf('G:\\mice\\%s\\dat\\%s\\%s',mouse_name,filename,filename); 
        return
    end
end