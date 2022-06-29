% [ Res ] = LoadRes(FileName)
%
% A simple matlab function to load a .res file

% 16-jan-12 ES

function Res = LoadRes( FileName )
 
fid = fopen( FileName, 'r' );
 
if fid == -1
    error( [ 'Could not open file ' FileName ] );
end
 
Res = fscanf( fid, '%d' );
fclose( fid );
 
return
