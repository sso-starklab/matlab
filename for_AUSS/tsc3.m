% tsc2 is the second version of the tsc function

% The tsc2 function is made to help the user with the frist few steps in
% the manual curation spike sorting.
%
% tsc2: takes as an input filebase, shank number, and numbers of channels in the specific shank and returns
% a re-write a clu file which units are odered form "good" to "bad".
% 
% to reorder the units in the clu file the function use an adaboost traind classifaier
% the inputs to the clasiffaier are:
% mean WF
% std of meanWF
% normelized auto CCG 

% Inputs:
% filebase        full path and session name e.g '/odin/Recordings/mF105/dat/mF105_10/mF105_10';
% ShankNum        shanknumber as a string e.g '4'
% channel_num     number of recordings saits in the shank e.g 11

%LS 12 March 20

function cleanClu =  tsc3(filebase,ShankNum,channel_num)
tic
[mainfb,endfb]  = fileparts(filebase);
LoadFn          = [mainfb,'/auto/',endfb,'.clu.',ShankNum];
cluS            = load(LoadFn);

% set empty featurs variable 
%nclu         = cluS(1);
cluS(1)      = [];
nclu         = length(unique(cluS));
temp         = [];
ch_feat      = [];
temp_feat    = [];
u            = unique(cluS);
cleanClu     = [];
% get features for the classifaier
X            = loadShankX2(filebase,ShankNum,channel_num);    

load XGtsc;

[tag,pred]  = predict(XGtsc2,X);
% idx         = pred(:,2)<=0.5;
% Tlist       = [u(idx),pred(idx,1)];
[i,m]        = max(pred,[],2);
pred2        = pred;

Npred        = zeros(size(pred,1),1);
idx          = m==3;

Npred(idx)   = pred2(idx,3)+((pred2(idx,3)-pred2(idx,1)))+10;
idx          = m==2;
idx2         = idx & pred2(:,3)>=pred2(:,1);
Npred(idx2)  = pred2(idx2,2)+5 +(pred2(idx2,3)-pred2(idx2,1)) ;

idx3         = idx & pred2(:,3)<pred2(:,1);
Npred(idx3)  = pred2(idx3,2)+2 +(pred2(idx3,1)-pred2(idx3,3));

idx          = m==1;
Npred(idx)   =1-( pred2(idx,1)+(pred2(idx,1)-pred2(idx,3)));
pred         = Npred;
% orgenize new clu file
[sortPred,idxP]     = sort(pred,'descend');
newClu              = zeros(length(cluS),1);
cleanClu            = zeros(length(cluS),1);
for i =1:nclu
   if u(idxP(i,1))>1 
   ind          = cluS == u(idxP(i,1));
   newClu(ind)  = i+1;
       if pred(idxP(i,1))>2
        cleanClu(ind)  = i+1;
       elseif pred(idxP(i,1))<=2
        cleanClu(ind)  = 0;
       end
   end
end

% write a new clu file
clufname        = [filebase,'.clu.',ShankNum];
fid             = fopen( clufname, 'w' );
[ ~ ]           = fprintf( fid, '%d\n', nclu );
[ ~ ]           = fprintf( fid, '%d\n', newClu );
[ ~ ]           = fclose(fid);

fprintf(['finshed procecing ' num2str(toc) ' sec' '\n' 'cutoff' num2str(sum(pred>=5)+2) '\n'])



 end








