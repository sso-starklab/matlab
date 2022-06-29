% if tsc == false (defaulte) - make sure there are no trash units in clu

function  []  = prepareData2(filebase,path,ShankNum,nchannels,varargin)

tsc                = false;
[tsc] = ParseArgPairs( {'tsc'}, {tsc}, varargin{ : } );

%currentPath               = pwd;
%p                         = mfilename('fullpath');
%path                      = fileparts(p);
path                      = [path,'/npy_files'];
cd(path);
fileStr       = split(filebase, '/');
sessionName   = fileStr{length(fileStr)};


if tsc
    cleanClu    = tsc3(filebase,ShankNum,nchannels);
else
    cleanClu    = load([filebase,'.clu.',ShankNum]);
    cleanClu(1) = [];
end

LoadFn      = [filebase,'.res.',ShankNum];
res         = load(LoadFn);

writeNPY(cleanClu,[sessionName,'.clu.',ShankNum,'.npy'])
writeNPY(res,[sessionName,'.res.',ShankNum,'.npy'])

[mspk,sspk,nspk_vec] = mWaveForm(filebase, nchannels, cleanClu, ShankNum);
id                   = (cleanClu==0 | cleanClu==1);
Uvec                 = unique(cleanClu);
Uid                  = (Uvec==0 | Uvec==1);
cleanClu(id)         = [];
res(id)              = [];
mspk(:,:,Uid)        = [];
sspk(:,:,Uid)        = [];
nspk_vec(Uid')       = [];


if nchannels < 8
    [mspk,sspk]      = addChanels(mspk, sspk, nchannels);
end

cc                   = CCG(res,cleanClu,20,20,20000,unique(cleanClu),'count');

writeNPY(mspk,[sessionName,'.mspk.',ShankNum,'.npy'])
writeNPY(sspk,[sessionName,'.sspk.',ShankNum,'.npy'])
writeNPY(nspk_vec,[sessionName,'.nspk_vec.',ShankNum,'.npy'])
writeNPY(cc,[sessionName,'.cc.',ShankNum,'.npy'])

end

