% sumrows               sums matrix rows by subsets
%
% call                  smat = sumrows( mat, rows );
%
% gets                  mat         matrix
%                       rows        source row indices, assumed to be sorted (verified internally)
%
% returns               smat        summed matrix, of dimensions
%                                       max( rows ) x size( mat, 2 )
%                                   each row corresponds to
%                                       smat( j, : ) = sum( mat( rows == j, : ), 1 )
%
% example:
% >> sumrows( [ 1 2 3; 4 5 6; 7 8 9; 10 11 12 ], [ 1 1 3 3 ] )
% ans =
%      5     7     9
%      0     0     0
%     17    19    21
%
% calls                 colvec
%
% see also              bincols

% NOTE: realized that MATLAB has a routine accumarray.m built in, check
% whether that can somehow substitute sumrows.m...

% 04-mar-13 ES
 
% revisions
% 04-nov-20 cleaned up

function smat = sumrows( mat, rows )

if nargin < 2 || isempty( mat ) || isempty( rows )
    error( '%s: missing arguments', upper( mfilename ) )
end

mat                             = colvec( mat );
[ m, n ]                        = size( mat );
rows                            = rows( : );
if ~issorted( rows )
    rows                        = sort( rows );
end
nrows                           = length( rows );
if m > nrows
    mat( nrows + 1, : )         = [];
end

borders                         = diff( [ rows; rows( end ) + 1 ] ) > 0;
cmat                            = cumsum( mat );
ridx                            = unique( rows );
smat                            = zeros( max( rows ), n );
smat( ridx, : )                 = diff( [ zeros( 1, n ); cmat( borders, : ) ] );

return

% EOF
