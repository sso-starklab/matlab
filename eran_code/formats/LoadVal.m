% [ Val ] = LoadVal( FileName )
%
% A simple matlab function to load a .val file (5 column ASCII file)
%
% Note: the format of a *.val.* file is
% [ trigger_onset trigger_offset onset_slope offset_slope mean_value ]
% onset/offset are in samples (of the original file, e.g. 20000)
% slopes and mean are in Volts (or could be A2DU)

% 27-jan-12 ES

function Val = LoadVal( FileName )

fid = fopen( FileName, 'r' );
 
if fid == -1
    error( [ 'Could not open file ' FileName ] );
end
 
Val = fscanf( fid, '%f', [ 5 inf ] )';
fclose( fid );

return

%EOF

filebase = '/Volumes/Data/phaser3/mouse365/25nov11/dat/es25nov11.013/es25nov11.013';
val = LoadVal( valfname );

