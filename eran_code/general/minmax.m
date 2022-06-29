% minmax        entire array

% 18-oct-12 ES

function y = minmax( x )

if nargin > 0
    y = [ min( x( : ) ) max( x( : ) ) ];
else
    y = [];
end

return