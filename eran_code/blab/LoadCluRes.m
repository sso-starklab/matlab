%
% A simple matlab function to load from the .clu and .res files the vectors G and T
% from many electrodes 
% Usage: [T,G,Map,Par]=LoadCluRes(FileBase,ElGpsToLoad, ClusToLoad)
% Enter the FileBase without the extension .clu or .res
% ElGpsToLoad - list of electrode groups to load (load them all by default)
% ClusToLoad allows to specifiy the vector of El numbers in ElGpsToLoad and
% matching length vector of Clu numbers (original indexing) such that only
% those are loaded
% alternatively ClusToLoad can specify ElClu pairs , and then ElGspToLoad
% can be just unique list of electrodes to load where these ElClu pairs are
% T is in samples, G goes from 1 to total number of clusters (excludes
% noise and artifacts )
% Map is a matrix displaying the correspondance between new cluster numbers (first column) and inital
% shank # (second column) and cluster # (third column)
% Pascale production, Anton just helped and made few additions :))
%
% 16-jan-12 (1) additional argument clustr enables loading from non *.clu.* files
%           (2) additional argument minCluToKeep enables selecging a subset
%           (3) FileBase overloading (par structure directly)

function [T,G,Map,Par]=LoadCluRes( FileBase, varargin )

if isa( FileBase, 'char' ) && exist( fileparts( FileBase ), 'dir' )
    Par = LoadPar([FileBase '.xml']);
elseif isa( FileBase, 'struct' )
    Par = FileBase;
    FileBase = Par.FileName;
else
    error( sprintf( '%s: input type mismatch!', upper( mfilename ) ) )
end

[ElGpsToLoad, ClusToLoad, clustr, minCluToKeep ] = DefaultArgs(varargin,{[1:Par.nElecGps],[],'clu',2 });

G=[];
T=[];
Map=[];
maxG=0;
clustr = sprintf( '.%s.', clustr );
clustr0 = '.clu.';

% Loop over x=ElGpsToLoad from LoadPar

for x=ElGpsToLoad(:)'
    % for each ElGp, load clu and res
    clufile = [FileBase clustr num2str(x)];
    if ~FileExists( clufile )
        clufile = [FileBase clustr0 num2str(x)];
        if ~FileExists( clufile )
            continue
        end
    end
    g = LoadClu( clufile );
    t = LoadRes( [FileBase '.res.' num2str(x)] );
    
    % Removes clusters artifact and noise clusters (0 & 1)
    indx = (g >= minCluToKeep );
    if sum(indx)==0
        continue;
    end;
    g=g(indx);
    t=t(indx);

    % creates vector of initial g and renames cluster # since 1 to n
    [ugini,b,g]=unique(g); % ugini: vector of initial unique values of g
    g=maxG+g;
    ug=unique(g);

    % concatenates all the g and t
    G=[G;g];
    maxG=max(G);
    T=[T;t];

    % Creates a "map" matrix
    shk=zeros(length(ug),1)+x; % electrode #
    map=[ug,shk,ugini];
    Map=[Map;map];
end
%sort the spikes not to have surprizes later -A
[T si] = sort(T);
G=G(si);

%now more fancy : if one specifies [El Clu] pairs to load
if ~isempty(ClusToLoad)
    if length(ClusToLoad)~=length(ElGpsToLoad) && size(ClusToLoad,2)~=2 
        warning('length(ClusToLoad)~=length(ElGpsToLoad)');
        return;
    end
   if size(ClusToLoad,2)~=2 
       myi = find(ismember(Map(:,2:3),[ElGpsToLoad(:), ClusToLoad(:)],'rows'));
   else
       myi = find(ismember(Map(:,2:3),ClusToLoad,'rows'));
   end
   gi = ismember(G, myi);    
   G = G(gi);
   T = T(gi);
   [uclu dummy G] = unique(G);
   Map = [[1:length(myi)]' Map(myi,2:3)];
   
end
    
    