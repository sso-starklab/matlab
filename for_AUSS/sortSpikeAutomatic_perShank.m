% shankNum - int, *not* str

function Nlist    = sortSpikeAutomatic_perShank(filebase, shankNum, varargin)

resFilebase             = fileparts(filebase);
grouping                = false;
tsc                = false;
[resFilebase, grouping, tsc] = ParseArgPairs( {'resFilebase','grouping','tsc'}, {resFilebase,grouping,tsc}, varargin{ : } );

    
par            = LoadXml( filebase );
%Nshanks        = size(par.SpkGrps,2);

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
[ ~ ]           = fprintf( fid, '%d\n', shankNum);
[ ~ ]           = fclose(fid);

nchannels       = length(par.SpkGrps(1,shankNum).Channels);
shankNumStr     = num2str(shankNum);
prepareData2(filebase,path1,shankNumStr,nchannels,'tsc',tsc);

 
path2 = '/media/shirly/Data/_Shirly/Lab/matlab/for_AUSS/sort/';        % filebase of the python files
cd(path2);
! source /home/shirly/anaconda3/etc/profile.d/conda.sh && conda activate auss && python sort_one_shank.py



clu       = readNPY([path_npy,sessionName,'.clu.',shankNumStr, '_2.npy']);
if ~grouping
    % organizing clu according to sorting (without grouping)
    Nlist     = cluOrgenize_py1(clu,filebase,resFilebase,sessionName,shankNumStr);
else
    % writing clu after grouping
    Nlist    = cluOrgenize_py2(clu,filebase,shankNumStr);
end


end 