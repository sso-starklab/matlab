% LoadInd               load data from *.ind.# file
%
% call                  [ ind, nChannels ] = LoadInd( filename )
%
% gets                  filename            full file name with path
%
% returns               ind                 the contents: [ nClu x nChannels ]
%                       nChannels           number of channels per clu
%
% calls                 nothing
% 
% conventions
%   this file starts with a scalar, indicating the number of channels per unit
%   for each unit, there is a distinct line with the channel numbers for
% which the waveform (in the *spk* file) and the features (in the *fet* file)
% correspond to. 
%   if there is a noise cluster, it must also correspond to one of the rows
% of the ind file.

% 04-sep-19 ES

function [ ind, nChannels ] = LoadInd( filename )

fp              = fopen( filename, 'r');
if fp == -1
    error( ['Could not open file ' filename] );
end

nChannels       = fscanf( fp, '%d', 1 );
ind             = fscanf( fp, '%f', [ nChannels, inf ] )';
fclose(fp);

return

% EOF
