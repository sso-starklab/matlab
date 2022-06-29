% outliers          a given number of SDs
%
% idx = outliers( x, nstd, bymad, rflag, dim )
%
% x         data; columns (rows,sheets..) of x are treated separately
% nstd      {3} number of deviations to treat as an outlier
% bymad     {0} by default, use mean/SD; if 1, use median/MAD
% rflag     {0} by default, return the indices of the outliers
%               if 1 and x is a vector, dilute the outliers
%               if 1 and x is a high-dim array, replace outliers with NaN
% dim       {1} dimension to work on

% see also      moutliers

% 18-nov-12 ES

% revisions
% 20-jun-13 dim added

function idx = outliers( x, nstd, bymad, rflag, dim )

sx = size( x );
rshp = 0;
if sum( sx( 2 : end ) == 1 )
    rshp = 1;
    x = x( : );
end
nargs = nargin;
if nargs < 2 || isempty( nstd )
    nstd = 3;
end
if nargs < 3 || isempty( bymad )
    bymad = 0;
end
if nargs < 4 || isempty( rflag )
    rflag = 0;
end
if nargs < 5 || isempty( dim )
    dim = 1;
end
if rflag == 1 && size( x, 2 ) > 1
    rflag = NaN;
end

if bymad
    mx = nanmedian( x, dim );
    sd = mad( x, 1, dim );
else
    mx = nanmean( x, dim );
    sd = nanstd( x, [], dim );
end
% b = [ 1 1 ]' * mx + nstd * [ -1 1 ]' * sd;
% idx = bsxfun( @lt, x, b( 1, : ) ) | bsxfun( @gt, x, b( 2, : ) );
sxhat = sx;
sxhat( : ) = 1; 
sxhat( dim ) = sx( dim );
mx1 = repmat( mx, sxhat );
w1 = nstd * repmat( sd, sxhat );
limLo = mx1 - w1;
limHi = mx1 + w1;
idx = x > limHi | x < limLo;

if rshp
    idx = reshape( idx, sx );
end
if rflag == 1
    x( idx ) = [];
    idx = x;
elseif isnan( rflag )
    x( idx ) = NaN;
    idx = x;
end

return


% EOF

x = randn( 1000, 3 );
idx = outliers( x, 3 );