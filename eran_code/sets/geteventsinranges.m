% geteventsinranges         according to a 2-column matrix of ranges
%
% call                      [ out, idx, midx ] = geteventsinranges( in, mat )
%
% input: 
% in        vector of discrete event times. 
% mat       range-specifying 2 column matrix, each row is 2
%               elements indicating a range. 
% flag      optional flag subtracts mat( j, 1 ) from each matching element
%
% output: 
% out:      vector, the data in the elements of in corresponding to the
%               ranges specified in mat (concatenated)
%               size of out will be sum( diff( mat, [], 2 ) + 1 ) x size( in, 2 )
% idx:      indices of in corresponding to out s.t. out = in( idx, : );
% midx:     indices of in corresponding to mat s.t. (for flag == 0)
%               out( midx == j ) >= mat( j, 1 ) && out( midx == j ) <= mat( j, 2 )
%               for any j.
%
% example:
% >> [ out idx ] = geteventsinranges( 10 : 10 : 100, [ 2 4; 6 8; 20 40 ], 0 ); [ out idx ]
% ans =
%     20     3
%     30     3
%     40     3
%
% this is a dummy routine, wrapper for inranges.m
%
% see also: dilutesegments, getdatainranges, inranges, intersectranges, isoverlap, plotranges, setdiffranges, sortranges, uniteranges

% 04-mar-13 ES

% revisions
% 23-mar-21 cleaned up

function [ out, idx, midx ] = geteventsinranges( in, mat, flag )

[ idx, midx, out ] = inranges( in, mat, flag );

return

% EOF