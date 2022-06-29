% mergeranges           by a count vector
%
% CALL                  [ mat, vec ] = mergeranges( mat, vec, mincount )
% 
% GETS                  mat
%                       vec
%                       mincount
%
% RETURNS               mat
%                       vec
% 
% DOES: 
%       mat and vec must have same number of rows/elements
%       for each count < minCount, the corresponding row of mat is merged with
%           its nearest neighbour (in terms of edges, not centers)
%
% CALLS                 parse, sortranges

% 28-may-13 ES

% revisions
% 17-aug-19 cleaned up

function [ mat, vec ] = mergeranges( mat, vec, mincount )

DEBUG                   = 0;

% argument check
if nargin < 3 || isempty( mat ) || isempty( vec ) || isempty( mincount )
    error( 'all arguments must be specified' )
end
if ~isequal( mat, sortranges( mat ) )
    error( 'mat must be a matrix of ranges' )
end
vec = colvec( vec );
if size( mat, 1 ) ~= size( vec, 1 )
    error( 'mat and vec must match' )
end
if length( mincount ) ~= 1
    error( 'mincount must be a scalar' )
end

if DEBUG
    vec0                = vec;
    mat0                = mat;
end

% algorithm
% (1) combine adjacent low-count ranges
t                       = vec < mincount;
if sum( t ) == 0
    return
end
if sum( t ) == length( t )
    mat                 = [ min( mat( : ) ) max( mat( : ) ) ];
    vec                 = sum( vec );
    if DEBUG
        fprintf( '\n%s: Pre:\n', upper( mfilename ) )
        disp( [ mat0 vec0 ] )
        fprintf( '%s: Post:\n', upper( mfilename ) )
        disp( [ mat vec ] )
    end
    return
end
x                       = parse( find( t ) );
idx                     = find( diff( x, 1, 2 ) > 0 );
for i                   = 1 : length( idx )
    ridx                = x( idx( i ), 1 ) : x( idx( i ), 2 );
    vec( ridx( 1 ) )    = sum( vec( ridx ) );
    vec( ridx( 2 : end ) ) = NaN;
    mat( ridx( 1 ), : ) = [ mat( ridx( 1 ), 1 ) mat( ridx( end ), 2 ) ];
    mat( ridx( 2 : end ), : ) = NaN;
end
if ~isempty( idx )
    rmv                 = setdiff( find( t ), x( :, 1 ) );
    mat( rmv, : )       = [];
    vec( rmv, : )       = [];
end

% (2) merge independent low-count ranges with their neighbours
t                       = vec < mincount;
if sum( t ) == 0 || sum( ~t ) == 0
    if DEBUG
        fprintf( '\nPre:\n' )
        disp( [ mat0 vec0 ] )
        fprintf( 'Post:\n' )
        disp( [ mat vec ] )
    end
    return
end
n                       = length( vec );
d                       = mat( 2 : n, 1 ) - mat( 1 : n - 1, 2 );
dpre                    = [ inf; d ];
dpost                   = [ d; inf ];
mpre                    = dpre( t ) <= dpost( t );
mpost                   = dpre( t ) > dpost( t );
ft                      = find( t );
x                       = ft( mpre );
for i                   = 1 : length( x )
    ridx                = x( i ) + [ -1 0 ];
    vec( ridx( 1 ) )    = sum( vec( ridx ) );
    vec( ridx( 2 ) )    = NaN;
    mat( ridx( 1 ), : ) = [ mat( ridx( 1 ), 1 ) mat( ridx( end ), 2 ) ];
    mat( ridx( 2 ), : ) = NaN;
end

x                   = ft( mpost );
for i = 1 : length( x )
    ridx                = x( i ) + [ 0 1 ];
    vec( ridx( 2 ) )    = sum( vec( ridx ) );
    vec( ridx( 1 ) )    = NaN;
    mat( ridx( 2 ), : ) = [ mat( ridx( 1 ), 1 ) mat( ridx( end ), 2 ) ];
    mat( ridx( 1 ), : ) = NaN;
end
mat( t, : )             = [];
vec( t, : )             = [];

if DEBUG
    fprintf( '\nPre:\n' )
    disp( [ mat0 vec0 ] )
    fprintf( 'Post:\n' )
    disp( [ mat vec ] )
end

return

% EOF
