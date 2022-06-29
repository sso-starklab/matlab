%organizing clu according to sorting (without merging)
%filebase1 - data folder
%filebase2 - writing folder


function Nlist   = cluOrgenize_py1(clu_sorted,filebase1,filebase2,sessionName,ShankNum)
Nlist      =[];
%  Load clu and create new clu
clu       = load([filebase1,'.clu.',ShankNum]);
clu(1)    = [];
newClu    = zeros(length(clu),1);
vecUnits  = unique(clu_sorted);

if length(vecUnits)>1
% Orgenize clu that units that need be mereged will have folowing numbers
for i =1:length(vecUnits)
    U   = vecUnits(i);
    if U > 1
        idx       = clu_sorted == U;
        Umerge    = unique(clu(idx));
        C         = [];
        for j  = 1 :length(Umerge)
            idx2          = clu == Umerge(j);
            Uname         = max(newClu) + 1;
            newClu(idx2)  = Uname;
            C             = [C,Uname];
        end
        Nlist{i}          = C + 1;
    end
end

newClu   = newClu+1;


% Orgenize the trash units
idxZ     = clu_sorted <= 1;
Zunits   = unique(clu(idxZ));

for i =1:length(Zunits)
     U           = Zunits(i);
     idx         = clu == U;    
     newClu(idx) = max(clu) + U + 100;
end

nclu      = length(unique(newClu));
clufname  = [filebase2,'/',sessionName,'.clu.',ShankNum];
fid2      = fopen( clufname, 'w' );
[ ~ ]     = fprintf( fid2, '%d\n', nclu);
[ ~ ]     = fprintf( fid2, '%d\n', newClu);
[ ~ ]     = fclose(fid2);
end
end
