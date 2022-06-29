% writing clu after sorting
% filebase1 -  writing folder

function Nlist = cluOrgenize_py2(clu_sorted,filebase1,shank_num)

nclu      = length(unique(clu_sorted));
clufname  = [filebase1,'.clu.',shank_num];
fid2      = fopen( clufname, 'w' );
[ ~ ]     = fprintf( fid2, '%d\n', nclu);
[ ~ ]     = fprintf( fid2, '%d\n', clu_sorted);
[ ~ ]     = fclose(fid2);
    
Nlist = [];
end