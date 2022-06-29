% moutliers     mean/median after removing outliers
%
% out = moutliers( v, fun, nstd, bymad, dim, args )
%
% to compute e.g. SEM without outliers (100 SD from the mean), type:
% >> moutliers( v, @calc_sem, 100 )
%
% notes:
% 1. outliers is over all dimesions
% 2. fun may be by dimension. it is called as follows:
%       fun( v, args{ : } )
% thus for instance if fun is @std and dim=2, args has to be 
%       [], 2
%
% example:
% xx is a 3D array. to remove outliers (10 SD) from the 3rd dimesion and
% average, call with:
%
%   xxm = moutliers( xx, @nanmean, 10, 0, 3, 3 );
%   
% see also outliers

% 28-may-13 ES

function out = moutliers( v, fun, nstd, bymad, dim, varargin )

nargs = nargin;
if nargs < 2 || isempty( fun ), fun = @mean; end
if nargs < 3 || isempty( nstd ), nstd = 10; end
if nargs < 4 || isempty( bymad )
    if isequal( fun, @median ) || isequal( fun, 'median' ) ...
            || isequal( fun, @nanmedian ) || isequal( fun, 'nanmedian' )
        bymad = 1;
    else
        bymad = 0;
    end
end
if nargs < 5 || isempty( dim ), dim = 1; end

vhat = outliers( v, nstd, bymad, 1, dim );
out = feval( fun, vhat, varargin{ : } );

return

% EOF
