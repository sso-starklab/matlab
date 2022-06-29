% get_rasters           for multiple units and events 
% 
% call                  [ r, bins ] = get_rasters( clu, res, trig, tim )
%
% gets                  clu, res      must be vectors of the same length
%                       trig, tim     must be vectors of the same length
%                       clu, trig     [labels], may overlap
%                       res, tim      [samples], same Fs
% 
% additional arguments, specified as parameter name/value pairs:
%                       binsize       {20}; [samples]
%                       halfwin       {50}; [bins]; use a 2-element vector to specify an assymetric range
%                       clunums       {[]}; labels to consider; defaults to unique( clu )
% 
% output:
%                       r         this generates 2D cell array of sparse matrices, each is: nbins x nreps
%                                 the array itself is nclu x ntrig. thus the cell array is m x n, and each
%                                 sparse matrix is ( p x qi ), where
%                                   m = length( unique( clu ) )
%                                   n = length( unique( trig ) )
%                                   p = 2 * halfwin + 1;
%                                   qi = sum( trig == utrig( i ) )
%                       bins      bin centers [samples]. to convert to seconds:
%                                     t = bins / Fs
% 
% does                  this routine cuts out spike times and deals them out to sparse arrays
%
% note:
% The default behavior is to cut data symmetrically around each trigger, for binsize * halfwin samples. 
% To cut a different arrangement, call with halfwin a 2-element vector of integers (may be positive and or negative)
% 
% calls                 inranges, ParseArgPairs, sumrows
%
% see also              get_spikes, multipeth

% 04-mar-13 ES

% revisions
% 14-oct-19 cleaned up and documented

function [ r, bins ] = get_rasters( clu, res, trig, tim, varargin )

r                       = [];
bins                    = [];

% arguments
nargs                   = nargin;
if nargs < 4 || isempty( clu ) || isempty( res ) 
    return
end
if isempty( trig ) && ~isempty( tim )
    trig                = ones( size( tim ) );
end
if length( clu ) ~= length( res ) || length( trig ) ~= length( tim )
    fprintf( '%s: input size mismatch\n', mfname );
    return
end
[ binsize, halfwin, clunums ] = ParseArgPairs(...
    { 'binsize', 'halfwin', 'clunums' }...
    , { 20, 50, [] }...
    , varargin{ : } );
if isempty( clunums )
    clunums         = unique( clu );
else
    kidx            = ismember( clu, clunums );
    res             = res( kidx );
    clu             = clu( kidx );
end
halfwin             = round( sort( halfwin ) );
if length( halfwin ) >= 2
    win             = halfwin( 1 : 2 );
end
if length( halfwin ) == 1
    win             = halfwin * [ -1 1 ];
end
bins                = ( win( 1 ) : win( end ) )' * binsize;

% remove irrelevant data
periods             = tim * [ 1 1 ] + ones( length( tim ), 1 ) * win * binsize;
kidx                = inranges( res, periods );
res                 = res( kidx );
clu                 = clu( kidx );
uclu                = unique( clu );

% initialize
[ ntrig, utrig ]    = uhist( trig );
n                   = length( utrig );
m                   = length( uclu );
r                   = cell( m, n );
p                   = diff( win ) + 1;
bidx                = round( ( win( 1 ) * binsize : win( 2 ) * binsize ) / binsize );
bidx                = bidx( : ) - bidx( 1 ) + 1;

% populate
pfull               = diff( win * binsize ) + 1;
for j               = 1 : n
    trignum         = utrig( j );
    per             = periods( trig == trignum, : );
    [ idx, midx, out ] = inranges( res, per, 1 );
    cluj            = clu( idx );
    for i           = 1 : m
        cidx        = cluj == uclu( i );
        mat         = sparse( out( cidx ), midx( cidx ), 1, pfull, ntrig( j ) ); % assume max 1 spikes/sample/clu
        bmat        = sparse( sumrows( mat, bidx ) );
        r{ i, j }   = bmat;
    end
end

% expand back to number of irrelevant units
eidx                = ismember( clunums, uclu );
if sum( eidx ) ~= length( eidx )
    z               = cell( 1, n ); 
    for j = 1 : n
        z{ j }      = sparse( p, ntrig( j ) ); 
    end
    r0              = r;
    r               = cell( length( clunums ), n );
    r( eidx, : )    = r0;
    r( ~eidx, : )   = repmat( z, [ sum( ~eidx ) 1 ] );
end

return

% EOF
