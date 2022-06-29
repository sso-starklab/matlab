% F_resampling_test      an a-parametric test for the null hypothesis (H0) that A and B have equal variances
%
% CALL                   [h,p]             = F_resampling_test( A, B, alpha, tail, nreps )
%
% GETS                   A                 value A for all datapoints
%                        B                 value B for all datapoints
% OPTIONAL           
%                        nreps             {500} number of permutations
%                        tail              {0} 0 - H1: variances are not equal
%                                              1 - H1: variance of A is larger than variance of B
%                                             -1 - H1: variance of A is smaller than variance of B
% CALLS                  mixmat
%
% RETURNS            
%                        h                  0 - do not reject H0; 1 - reject H0 
%                        p                  p value 
%
% see also               ftest

% written by             HS + ES  17-04-20

% revisions
% 21-apr-20 (1) full support of RAM usage optimization
%           (2) modified test statistic for consistency with ftest

function [h,p]              = F_resampling_test( A, B, alpha, tail, nreps )

MAXSIZE                     = 1e7;

%--------------------------------------------------------------------%
% check inputs
%--------------------------------------------------------------------%


% argument handling
nargs                       = nargin;
if nargs < 2 || isempty( A ) || isempty( B )
   error( 'missing arguments' )
end
if nargs < 3 || isempty( alpha )
    alpha               = 0.05; 
end
if alpha <= 0 || alpha >= 1 
    error( 'alpha must be within the range 0-1' )
end
if nargs < 4 || isempty( tail )
    tail                    = 0;
end
if ~ismember( tail, [ -1 0 1 ] )
    error( 'tail must be within the set {-1,0,1}' )
end
if nargs < 5 || isempty( nreps )
    nreps                   = 500;
end


%--------------------------------------------------------------------%
% Generate permuted matrix
%--------------------------------------------------------------------%

% get sizes
A                           = A( ~isnan( A ), : );
B                           = B( ~isnan( B ), : );
sizeA                       = size(A,1);
sizeB                       = size(B,1);

% concatenate data
MAT                         = [A; B];

% permute the indexing vector
mat_size                    = ( nreps + 1 ) * ( sizeA + sizeB );

if mat_size < MAXSIZE
    m                       = ( 1 : ( sizeA + sizeB ) )';
    M                       = mixmat( repmat( m, [ 1 nreps ] ), 1, 1 );
    M                       = [ m M ];
else
    m                       = [ false( sizeA, 1 ); true( sizeB, 1 ) ];
end


%--------------------------------------------------------------------%
% Calculate statistic per each permutation
%--------------------------------------------------------------------%

if mat_size < MAXSIZE
    idxA                    = M( 1 : sizeA, : );
    idxB                    = M( sizeA + ( 1 : sizeB ), : );
    matA                    = MAT( idxA );
    matB                    = MAT( idxB );
    varA                    = var( matA );
    varB                    = var( matB );
    s                       = log( varB ./ varA );
else
    s                       = zeros(nreps+1,1);
    for i = 1 : nreps + 1
        if i == 1
            idx             = m;
        else
            idx             = mixmat(m,1,1);
        end
        A_temp              = MAT(~idx);
        B_temp              = MAT(idx);
        varA                = var( A_temp );
        varB                = var( B_temp );
        s( i )              = log( varB ./ varA );
    end
end

%--------------------------------------------------------------------%
% Calculate p value
%--------------------------------------------------------------------%
s0                          = s( 1 );
sP                          = s( 2 : end );

p1                          = ( sum( s0 >= sP ) + 1 ) / ( nreps + 1 );
p2                          = ( sum( s0 <= sP ) + 1 ) / ( nreps + 1 );

switch tail
    case 0
        p                   = 2 * min( p1, p2 );
    case 1
        p                   = p1;
    case -1
        p                   = p2;
end

% binary test result
h                       = p <= alpha;

return

% EOF
