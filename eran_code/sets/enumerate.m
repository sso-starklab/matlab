% enumerate             the inverse of parse
%
% call                  vec = enumerate( mat )
%
% call example
% 
% >> mat = [ 1 4; 7 9 ];
% >> vec = enumerate( mat )
% 
% yields
% >> [ 1 2 3 4 7 8 9 ]
%
% in general, should always have
% >> isequal( mat, parse( vec ) )
%
% calls                 sortranges
% 
% see also              dilutesegments, getdatainranges, geteventsinranges, inranges, intersectranges, isoverlap, plotranges, setdiffranges, uniteranges

% 14-mar-13 ES

% revisions
% 18-jun-20 cleaned up

function vec = enumerate( mat )

% initialize output
vec                         = [];

% arguments
if isempty( mat )
    return
end

% code
mat                                 = sortranges( mat );
m                                   = size( mat, 1 );
durs                                = diff( mat, [], 2 ) + 1;
cdurs                               = cumsum( durs );
vec                                 = zeros( cdurs( m ), 1 );
idx                                 = [ 1; cdurs( 1 : m - 1 ) + 1 ];
for i                               = 1 : m
    vec( idx( i ) : cdurs( i ) )    = mat( i, 1 ) : mat( i, 2 );
end

return

% EOF
