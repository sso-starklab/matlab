% mstrfind      multi strfind
% 
% [ tf idx ] = mstrfind( str, strs )
%
% str       a single string
% strs      cell array of strings
% 
% tf        0 if none are found, 1 if any
% idx       indices of found strings. indices in str are not returned

% 11-jul-13 ES

function [ tf, idx ] = mstrfind( str, strs )

tf          = false( 1 );
idx         = [];
if ~isa( str, 'char' )
    return
end
if isa( strs, 'char' )
    strs    = { strs };
end
for i       = 1 : length( strs )
    if ~isempty( strfind( str, strs{ i } ) )
        tf  = 1;
        idx = [ idx; i ];
    end
end

return

% EOF