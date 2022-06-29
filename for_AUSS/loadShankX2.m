
% This function get path for one session and the shank bunber e.g. 
% path   = 'odin\Recordings\mF108\dat\mF108_01\'
% shank  = '1'
% The function returns the matrix X which is a matrix of featurs for eac cluster in the clu A files
% features are:
% mean WF using 8
% std of WF using 8 channels
% auto CCH
% the function also returns the Y label for each clusters
% 1 is D and below (zero or one clusters)
% 2 is C groups
% 3 is B and above

% rev changed the filebase input to 'odin\Recordings\mF108\dat\mF108_01\mF108_01\' 13_12_20 TR
% New Inputs:
% filebase        full path and session name e.g '/odin/Recordings/mF105/dat/mF105_10/mF105_10';
% shank           shanknumber as a string e.g '4'
% channel_num     number of recordings saits in the shank e.g 11



function  X    = loadShankX2(filebase,shank,channel_num)


[mainfb,endfb]  = fileparts(filebase);
LoadFn          = [mainfb,'/auto/',endfb,'.clu.',shank];
clu            = load(LoadFn);

%clu                 = load([filebase,'.clu.',shank]);
clu(1)              = [];
res                 = load([filebase,'.res.',shank]);

vecOfCulsters       = unique(clu);
xtag                = [];
[mean_spk,std_spk]  = mWaveForm(filebase, channel_num, clu, shank);

if channel_num < 8
    [mean_spk,std_spk]      = addChanels(mean_spk, std_spk, channel_num);
end


for j =1:numel(unique(clu))
    
    group           = vecOfCulsters(j);
    if group~=0 &&  group~=1
    idxcluA         = clu == group;
   
    % computing mean WF and std 
    [mean_spk2,ind] = trimSpk_8ch(mean_spk(:,:,j));
    std_spk2        = std_spk(ind,:,j);
    
   % computing Auto CCH
    cch1            = CCG(res(idxcluA),clu(idxcluA),20,30,20000,group,'count');
    cch             = cch1./max(cch1);
    xtag(j,:)       = [reshape(mean_spk2',1,[]),reshape(std_spk2',1,[]),cch',sum(cch1)];
    end
 
end
% idx0              = find(sum(xtag,2)~=0);
X                   = xtag;

return
