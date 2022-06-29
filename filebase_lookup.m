function filebase = filebase_lookup(filename, database)

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( filename )
    return
end
if nargs < 2 || isempty( database )
    database = 1;
end

aa = find(filename(:) == '_');
mouse_name = filename(1:aa-1);

if database == 1
    if isunix
        if size(mouse_name,2) < 5 
            filebase = sprintf('/probox1/mice/%s/dat/%s/%s',mouse_name,filename,filename);
            return
        else
            filebase = sprintf('/probox1/mice/%s/dat/%s/%s',mouse_name,filename,filename); 
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
end

if database == 2
    if isunix
        if size(mouse_name,2) < 5 
            filebase = sprintf('/probox2/mice/%s/dat/%s/%s',mouse_name,filename,filename);
            return
        else
            filebase = sprintf('/probox2/mice/%s/dat/%s/%s',mouse_name,filename,filename); 
            return
        end
    elseif ispc
        if size(mouse_name,2) < 5 
            filebase = sprintf('F:\\mice\\%s\\dat\\%s\\%s',mouse_name,filename,filename);
            return
        else
            filebase = sprintf('F:\\mice\\%s\\dat\\%s\\%s',mouse_name,filename,filename); 
            return
        end
    end
end