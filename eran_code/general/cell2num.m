% CELL2NUM      assign an integer to a string in a cell array.
%
% nums = cell2num( cell, key )
%
% see also      num2cellK

% 28-may-13 ES

function nums = cell2num( cell, key )

sx = size( cell );
cell = cell(:);
if nargin < 2 || isempty( key )
    key = unique( cell );
end
key = key(:);
Nc = length(cell);
Nk = length(key);
nums = zeros(Nc,1);

for i = 1:Nc
    for j = 1:Nk
        if strcmp(cell{i},key{j})
            nums(i) = j; 
            break
        end
    end
end

nums = reshape( nums, sx );

return

% EOF
