function Nlist_all    = sortSpikeAutomatic5(filebase, varargin)

resFilebase             = fileparts(filebase);
grouping                = false;
[resFilebase, grouping] = ParseArgPairs( {'resFilebase','grouping'}, {resFilebase,grouping}, varargin{ : } );

par            = LoadXml( filebase );
Nshanks        = size(par.SpkGrps,2);

fileStr         = split(filebase, '/');
sessionName     = fileStr{length(fileStr)};
newFile         = extractBefore(filebase, sessionName);
path1           = [newFile ,sessionName];
path_npy        = [path1, '/npy_files/'];

cd(path1);
if exist(path_npy, 'dir')
    rmdir npy_files s;
end

mkdir npy_files;


fname           = [path_npy,'/info'];
fid             = fopen( fname, 'w' );
[ ~ ]           = fprintf( fid, '%s\n', sessionName);
[ ~ ]           = fprintf( fid, '%d\n', Nshanks);
[ ~ ]           = fclose(fid);

 for i =1:Nshanks
     nchannels       = length(par.SpkGrps(1,i).Channels);
     shankNum        = num2str(i);
     prepareData(filebase,path1,shankNum,nchannels);
 end
 

path2 = '/media/shirly/Data/_Shirly/Lab/matlab/for_AUSS/sort/';        % filebase of the python files
cd(path2);
%! source /home/tali/anaconda3/etc/profile.d/conda.sh && conda activate auss && python sort_shank.py
system(sprintf('source /home/shirly/anaconda3/etc/profile.d/conda.sh && conda activate auss && python sort_shank.py %s',path_npy))

for i= 1:Nshanks
    shank_num = num2str(i);
    clu       = readNPY([path_npy,sessionName,'.clu.',shank_num, '_2.npy']);
    if ~grouping
        % organizing clu according to sorting (without grouping)
        Nlist     = cluOrgenize_py1(clu,filebase,resFilebase,sessionName,shank_num);
    else
        % writing clu after grouping
        Nlist    = cluOrgenize_py2(clu,filebase,shank_num);
    end
    Nlist_all{i} = Nlist;

end


end   