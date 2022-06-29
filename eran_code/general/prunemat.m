% prunemat          prune a matrix of pairwise indices
%
% pmat = prunemat( mat )
%
% mat: matrix with each column indicating a set. each row includes a matched
% pair of elements plus zeros for the other sets. the rows of the pruned
% matrix include multi-set matched elements
%
% assumptions:
% -input matrix is the outcome of serial pair-wise comparisons. i.e. each
% row has exactly two non-zero elements
% -input matrix is a matrix of indices, i.e. all elements are either zero
% or non-negative integers
%
% example
% >> mat = [ 2 1 0 0 0
%         2 0 3 0 0
%         0 1 0 3 0
%         0 0 3 3 0
%         4 4 0 0 0
%         4 0 5 0 0
%         0 6 8 0 0
%         5 0 0 4 0 ];
% >> pmat = prunemat( mat )
% pmat =
%      2     1     3     3     0
%      4     4     5     0     0
%      5     0     0     4     0
%      0     6     8     0     0
%
% calls: nothing

% 17-jan-13 ES

function y = prunemat( x )

[ m, n ]            = size( x );
y                   = [];
i                   = 1;

while m > 0
    
    % handling of a single element:
    val             = x( i );
    if val == 0
        i           = find( x( : ) > 0, 1 );
        continue;
    end
    col             = ceil( i / size( x, 1 ) );
    rows            = find( x( :, col ) == val );
    cols            = [ 1 : col - 1 col + 1 : n ];
    mat             = x( rows, cols );
    
    % now expand the matrix from this point:
    oldrows         = [];
    while ~isequal( oldrows, rows )
        oldrows     = rows;
        for j       = 1 : ( n - 1 )
            v       = unique( mat( :, j ) );
            v       = v( v > 0 );
            if isempty( v )
                continue
            elseif length( v ) > 1
                v   = max( v );
            end
            arows   = find( x( :, cols( j ) ) == v );
            rows    = unique( [ rows; arows ] );
        end
    end
    
    % update the connected graph
    newrow          = zeros( 1, n );
    for k           = 1 : n
        v           = unique( x( rows, k ) );
        v           = v( v > 0 );
        if isempty( v )
            continue
        elseif length( v ) > 1
            v       = max( v );
        end
        newrow( k ) = v;
    end
    y               = [ y; newrow ];
    % remove the exhausted rows and advance index
    x( rows, : )    = [];
    m               = size( x, 1 );
    if i > numel( x )
        i           = 1; 
    end
    if m > 0 && x( i ) == 0
        i           = i + 1;
    end
    if i > numel( x )
        i = 1; 
    end
    
end

return

% EOF