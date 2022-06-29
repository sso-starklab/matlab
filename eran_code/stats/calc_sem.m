% CALC_SEM          matrix sem by dimension dim.
%
% call              Y = CALC_SEM( X, DIM )
% 
% gets              X       data
%                   DIM     dimension {1} - sem of columns

% 08-mar-03 ES

function y = calc_sem( x, dim )

if all(all(isnan(x))),
    if nargin < 2 || isempty( dim ), dim = 1; end
    y = nanmean(x,dim);
    return 
end 
if ~exist( 'dim', 'var' ) || isempty( dim )
    if any( size( x ) == 1 )
        x = x(:);
    end
    dim = 1;
end

y = nanstd( x, [], dim ) ./ sqrt( sum( ~isnan( x ), dim ) );

return