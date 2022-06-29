% isoverlap              between two 2-column matrices of ranges
%
% [ tf, loc ] = isoverlap( mat1, mat2 )
% 
% input:	2 column matrices, each row is set of two numbers
%           indicating a range. 
%           the matrices can have different number of rows.
%           arrays do not have to be sorted
%
% tf:       logical vector, same length as the number of mat1 rows. 1 if the
%               range in that row overlaps some range in mat2, 0 otherwise. 
% loc:      same dimensions as tf; the highest index of mat2 for which the 
%               overlap occurs, 0 if no overlap.
% 
% note: isoverlap( mat1, mat2 ) is NOT the same as isoverlap( mat2, mat1 )
% 
% example:
% mat1 = [ 1 22; 20 30; 40 50 ]; 
% mat2 = [ 4 6; 25 35; 60 70 ];
% [ tf1 loc1 ] = ismember( fliplr( mat1 ), mat2, 'rows' )
% [ tf2 loc2 ] = isoverlap( mat1, mat2 )
%
% see also: sortranges, dilutesegments, getdatainranges, geteventsinranges, inranges, intersectranges, plotranges, setdiffranges, uniteranges

% 26-jan-12 ES

% revisions
% 03-mar-13 efficiency (fast linear algorithm)

function [ tf, loc ] = isoverlap( mat1, mat2, flag )

nargs = nargin;
if nargs < 1 || isempty( mat1 )
    tf = [];
    loc = [];
    return;
end
if nargs < 3 || isempty( flag )
    flag = 1;
end
mat1 = sortranges( mat1, flag );
m1 = size( mat1, 1 );
tf = false( m1, 1 );
loc = zeros( m1, 1 );
if nargs < 2 || isempty( mat2 )
    return
end
mat2 = sortranges( mat2 );
m2 = size( mat2, 1 );

j = 1;
k = 0;
for i = 1 : m1
    j = j - k;
    k = 0;
    g = [];
    while j <= m2 && mat2( j, 1 ) <= mat1( i, 2 )
        if mat1( i, 2 ) >= mat2( j, 1 ) && mat1( i, 1 ) <= mat2( j, 2 )
            k = k + 1;
            g( k ) = j;
        end
        j = j + 1;
    end
    if ~isempty( g )
        tf( i ) = 1;
        loc( i ) = max( g );
    end
end

return

% EOF

