% getdatainranges        according to a 2-column matrix of ranges
%
% call                  [ out, idx, idx2, idx3, omat ] = getdatainranges( in, mat, pad )
%
% input: 
% in        data matrix with continuous data
% mat       range-specifying 2 column matrix, each row is 2
%               elements indicating a range. the matrix does not have to be sorted
% pad       inter-segment padder, must have same number of columns as in
%
% output: 
% out:      single matrix; the data in the elements of in corresponding to the
%               ranges specified in mat (concatenated)
%               size of out will be sum( diff( mat, [], 2 ) + 1 ) x size( in, 2 )
% idx:      indices of in corresponding to out s.t. out = in( idx, : );
% idx2:     indices for the inverted operation
% idx3:     indices for the corresponding range
% omat:     matrix indicating ranges (same format as mat)
%
% example:
% >> [ out idx idx2 idx3 ] = getdatainranges( 10 : 10 : 100, [ 2 4; 6 8; 20 40 ] ); [ out idx idx2 idx3 ]
% ans =
%     20     2     1     1
%     30     3     2     1
%     40     4     3     1
%     60     6     4     2
%     70     7     5     2
%     80     8     6     2
%
% see also: dilutesegments, geteventsinranges, inranges, intersectranges, isoverlap, plotranges, setdiffranges, sortranges, uniteranges

% 17-jul-12

% revisions
% 24-nov-12 bug fix <= n
% 28-jan-13 idx optional output added
% 03-mar-13 mat rows sorted
% 19-feb-14 1. support inter-segment padding
%           2. mirror back indices, s.t. 
%              >> isequal( in( idx, : ), out( idx2, : ) )
%              is true
% 10-dec-19 added idx3 that indicates to which range each value corresponds
% 12-dec-19 1. added omat that indicates the ranges concisely
%           2. cleaned up

function [ out, idx, idx2, idx3, omat ] = getdatainranges( in, mat, pad )

nargs                   = nargin;
nout                    = nargout;
if nargs < 2 || isempty( in ) || isempty( mat )
    out                 = []; 
    return
end
if nargs < 3 || isempty( pad )
    pad                 = []; 
end
if size( mat, 2 ) ~= 2
    error( '%s: not a matrix of ranges\n', upper( mfilename ) )
end
in                      = colvec( in );
n                       = size( in, 1 );
mat( sum( mat > n, 2 ) > 0, : ) = [];
if isempty( mat )
    out                 = [];
    idx                 = [];
    return
end
m                       = size( mat, 1 );
if sum( mat( :, 1 ) <= mat( :, 2 ) ) ~= m
    sidx                = mat( :, 1 ) > mat( :, 2 );
    t                   = mat( sidx, : );
    mat( sidx, : )      = t( :, [ 2 1 ] );
end
if ~isempty( pad ) && size( pad, 2 ) == size( in, 2 )
    npad                = size( pad, 1 );
else
    npad                = 0;
end
durs                    = diff( mat, [], 2 ) + 1;
eidx                    = cumsum( durs );
sidx                    = [ 1; eidx( 1 : m - 1 ) + 1 ];
out                     = zeros( min( eidx( m ) + npad * ( m - 1 ), n ), size( in, 2 ) );
if nout > 1
    idx                 = zeros( eidx( m ), 1 );
end
if nout > 2
    idx2                = idx;
end
if nout > 3
    idx3                = idx;
end
if nout > 4
    omat                = [ sidx eidx ];
end

for i                   = 1 : m
    if mat( i, 2 ) > n
        continue
    end
    tidx                = mat( i, 1 ) : mat( i, 2 );
    % padding:
    if npad
        idxI            = ( sidx( i ) : eidx( i ) ) + npad * ( i - 1 );
        out( eidx( i ) + npad * ( i - 1 ) + ( 1 : npad ), : ) = pad;
    else
        % no padding:
        idxI            = sidx( i ) : eidx( i );
    end
    out( idxI, : )      = in( tidx, : );
    if nout > 1
        % no padding:
        idx( sidx( i ) : eidx( i ) ) = tidx( : );
        % padding:
        %idx( idxI ) = tidx( : );
    end
    if nout > 2
        idx2( sidx( i ) : eidx( i ) ) = idxI;
    end
    if nout > 3
        idx3( sidx( i ) : eidx( i ) ) = i * ones( length( tidx ), 1 );
    end
end

if issparse( in )
    out                 = sparse( out );
end

return

% EOF