% num2strs          convert numeric array to cell array of strings
%
% call              c = num2strs( mat, format )
%
% gets              mat         the numeric array (any number of elements/dimensions)
%                   format      argument to num2str
%
% see also          num2str

% 24-jun-20 ES & SSo

function c = num2strs( mat, format )

nargs                   = nargin;
if nargs < 1 || isempty( mat )
    c                   = [];
    return
end
if nargs < 2 || isempty( format )
    format              = [];
end

sz                      = size( mat );
c                       = cell( sz );
n                       = numel( mat );
for i                   = 1 : n
    x                   = mat( i );
    if isempty( format )
        y               = num2str( x );
    else
        y               = num2str( x, format );
    end
    c{ i }              = y;
end

return

% EOF