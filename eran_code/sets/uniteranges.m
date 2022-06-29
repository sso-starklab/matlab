% uniteranges           of two 2-column matrices of ranges
%
% call                  [ umat ] = uniteranges( mat1, mat2 )
%
% input: 2 column matrices, each row is set of two numbers
% indicating a range. the matrices can have different number of rows.
% arrays do not have to be sorted
%
% output: 
% umat: single matrix; sorted
% size will be size( mat1, 1 ) + size( mat2, 1 ) - sum( isoverlap( mat1, mat2 ) )
% 
% examples:
% mat1 = [ 5 10; 15 20; 25 30 ]; mat2 = [ 6 11; 19 26 ];
% mat1 = [ 5 10; 15 20; 25 30 ]; mat2 = [ 2 3; 6 11; 19 26 ];
% mat1 = [ 1 22; 20 30; 40 50 ]; mat2 = [ 4 6; 25 35; 60 70 ];
% mat1 = [ 10 20; 50 60; 55 60; 70 80 ]; mat2 = [ 10 20; 55 60; 70 80 ];
% umat = uniteranges( mat1, mat2 )
% 
% see also: dilutesegments, getdatainranges, geteventsinranges, inranges, intersectranges, isoverlap, plotranges, setdiffranges, sortranges

% 27-jan-12 ES

% revisions
% 03-mar-03 efficiency, bug fixes
% for large datasets, this routine is hundreds of times faster than the naive implementation

function [ umat ] = uniteranges( mat1, mat2 )

nargs = nargin;
if nargs < 1 || isempty( mat1 )
    umat = [];
    return
end
mat1 = sortranges( mat1 );
if nargs < 2 || isempty( mat2 )
    umat = mat1; 
    return
end
umat = sortranges( [ mat1; mat2 ] );

return

% EOF
