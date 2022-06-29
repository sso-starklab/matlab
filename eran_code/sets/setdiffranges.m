% setdiffranges         of two 2-column matrices of ranges
%
% CALL                  [ dmat, idx ] = setdiffranges( mat1, mat2 )
%
% GETS                  two 2-column matrices, each row is set of two numbers
%                           indicating a range. the matrices can have different number of rows.
%                           arrays do not have to be sorted
%
% RETURNS               dmat:     single matrix; sorted. all ranges that are in mat1 but not in mat2
% 
% EXAMPLE:
%               mat1 = [ 5 10; 15 20; 25 30 ]; 
%               mat2 = [ 6 11; 19 26 ];
%               dmat = setdiffranges( mat1, mat2 )
% 
% CALLS                 sortranges
% 
% See also              dilutesegments, getdatainranges, geteventsinranges, inranges, intersectranges, isoverlap, plotranges, sortranges, uniteranges

% 09-oct-12

% revisions
% 03-mar-13 efficiency
% 29-jul-13 indiced returned (preliminary)
% 17-aug-19 cleaned up

function [ dmat, idx ] = setdiffranges( mat1, mat2 )

nargs                       = nargin;
if nargs < 1 || isempty( mat1 )
    dmat                    = [];
    return
end
mat1                        = sortranges( mat1 );
if nargs < 2 || isempty( mat2 )
    dmat                    = mat1; 
    return
end
mat2                        = sortranges( mat2 );

m1                          = size( mat1, 1 );
m2                          = size( mat2, 1 );
dmat                        = zeros( m1 + m2, 2 );
h                           = 1;
j                           = 1;
k                           = 0;
idx                         = zeros( m1, 1 );
lidx                        = 0;
for i                       = 1 : m1
    j                       = j - k;
    k                       = 0;
    g                       = [];
    while j <= m2 && mat2( j, 1 ) <= mat1( i, 2 )
        if mat1( i, 2 ) >= mat2( j, 1 ) && mat1( i, 1 ) <= mat2( j, 2 )
            k               = k + 1;
            g( k )          = j;
        end
        j                   = j + 1;
    end
    if ~isempty( g )
        segs                = mat2( g, : );
        ns                  = size( segs, 1 );
        if segs( 1, 1 ) > mat1( i, 1 )    
            b1              = [ mat1( i, 1 ) segs( 1, 1 ) - 1 ];
        else
            b1              = [];
        end
        b2                  = [ segs( 1 : ns - 1, 2 ) + 1 segs( 2 : ns, 1 ) - 1 ];
        if mat1( i, 2 ) > segs( ns, 2 )
            b3              = [ segs( ns, 2 ) + 1 mat1( i, 2 ) ];
        else
            b3 = [];
        end
        b                   = [ b1; b2; b3 ];
        ns                  = size( b, 1 );
        dmat( h : ( h + ns - 1 ), : ) = b;
        h                   = h + ns;
        if ~isempty( b )
            %keyboard
            %fprintf( '%d \t', i )
        end
    else
        dmat( h, : )        = mat1( i, : );
        h                   = h + 1;
        lidx                = lidx + 1;
        idx( lidx )         = i;
    end
end
dmat( h : ( m1 + m2 ), : )  = [];
idx( lidx + 1 : m1 )        = []; % only correct if mat1 was sorted

return

% EOF
