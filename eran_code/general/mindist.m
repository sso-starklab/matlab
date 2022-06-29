% mindist           find indices of minimal pair-wise Euclidean distances between two points
%
% call              idx1 = mindist( v1, v2 )
% 
% gets              vector v1
%                   vector v2 
%
% return            idx of v1 that correspond to minimal distance to the
%                   points in v2
% 
% note              for unique idx, must have: length( v2 ) <= length( v1 )
%
% calls             nothing

% 10-sep-19 ES

function idx1 = mindist( v1, v2 )

v1              = v1( : );
v2              = v2( : )';
nv1             = length( v1 );
nv2             = length( v2 );
mat1            = v1 * ones( 1, nv2 );
mat2            = ones( nv1, 1 ) * v2;
dmat            = mat1 - mat2;
[ ~, idx1 ]     = min( abs( dmat ) );

return

% EOF

