% intersectranges       of two 2-column matrices of ranges
%
% call                  [ imat i1 i2 ] = intersectranges( mat1, mat2, flag )
%
% input:                2 column matrices, each row is set of two numbers
%                           indicating a range. 
%                           the matrices can have different number of rows.
%                           arrays do not have to be sorted
%
% optional input:       flag, {1}, passed to sortranges
%
% imat:                 single matrix; sorted
%                       containing all the ranges in both mat1 and mat2
% i1:                   for non-self-overlaps, imat = mat1( i1, : );
% i2:                   :, imat = mat2( i2, : );
%
% example:
% mat1 = [ 5 10; 15 20; 25 30 ]; mat2 = [ 6 11; 19 26 ];
% mat1 = [ 1 20; 25 30; 40 50 ]; mat2 = [ 5 8; 12 15; 18 28 ];
% imat = intersectranges( mat1, mat2 )
%
% see also: dilutesegments, getdatainranges, geteventsinranges, inranges, isoverlap, plotranges, setdiffranges, sortranges, uniteranges

% 17-jul-12

% revisions
% 06-jan-13 modified to be empty intersection if any is empty
% 03-mar-13 efficiency (fast linear algorithm)
% 10-sep-19 added flag to pass to sortranges

function [ imat, i1, i2 ] = intersectranges( mat1, mat2, flag )

% initialize output
i1              = [];
i2              = [];
imat            = [];

% process arguments
nargs = nargin;
if nargs < 2 || isempty( mat1 ) || isempty( mat2 )
    return
end
if nargs < 3 || isempty( flag )
    flag        = 1;
end
if nargout > 1
    kidx        = 1;
else
    kidx        = 0;
end

% prepare
mat1            = sortranges( mat1, flag );
mat2            = sortranges( mat2, flag );
m1              = size( mat1, 1 );
m2              = size( mat2, 1 );

% intersect
imat            = zeros( max( size( mat2, 1 ), m1 ), 2 );
h               = 0;
j               = 1;
k               = 0;
for i           = 1 : m1
    j           = j - k;
    k           = 0;
    g           = [];
    while j <= m2 && mat2( j, 1 ) <= mat1( i, 2 )
        if mat1( i, 2 ) >= mat2( j, 1 ) && mat1( i, 1 ) <= mat2( j, 2 )
            k   = k + 1;
            g( k ) = j;
        end
        j       = j + 1;
    end
    if ~isempty( g )
        if kidx
            i1  = [ i1; i ];
            i2  = [ i2 g ];
        end
        segs    = mat2( g, : );
        for si  = 1 : size( segs, 1 )
            h   = h + 1;
            imat( h, : ) = [ max( mat1( i, 1 ), segs( si, 1 ) ) min( mat1( i, 2 ), segs( si, 2 ) ) ];
        end
    end
end

if h == 0
    imat        = [];
else
    imat( h + 1 : end, : ) = [];
    if kidx
        i1      = unique( i1 );
        i2      = unique( i2( : ) );
    end
end

return

% EOF


