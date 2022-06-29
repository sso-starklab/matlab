% st2sparse             cell array of sparse to pointer + sparse
%
% call                  [ st tmat ] = st2sparse( r, uclu )
%
% gets                  r               cell array, nUnits x 1; each is nSamples x nTrials
%                       uclu            vector, nUnits x 1
%
% optional arguments (name/value pairs):
%                       map             matrix, nUnits x 3 (combined clunum, shank, clu)
%                       shankclu        matrix, nUnits x 3 (shank, clu, celltype)
%                       Detailed        {1}; 1 requires that shankclu and map will be supplied; 0 requires only uclu 
%                       doMUA           {0}; 0 returns SU matrix; 1 agglomerates into same-shank i/pMUA
%
% returns               st              nSamples x nInstances
%                       tmat            Detailed:   nInstances x 4 [ shanknum, clunum, celltype, trial index ]
%                                       otherwise:  nInstances x 2 [ combined clunum, trial index ]
%
% calls                 ParseArgPairs
%
% see also              gwnAnalysis1
%
% Notes: 
% 1. MUA requires Detailed information (i.e. shankclu and map must be supplied)
% 2. MUA can in principle result in st with entries other than 0/1, if two units spiked at the same time
%       (not likely for standard detection algorithms d.t. dead time, but possible for template-matching)

% 15-jun-14 ES

% revisions
% 16-jun-14 added nsources: 
%           uShanks x nTypes, counting the number of source clusters in each MUA
% 14-oct-19 cleaned up and documented

function [ st, tmat, nsources ] = st2sparse( r, uclu, varargin )

st                          = [];
tmat                        = [];
nsources                    = [];

%----------------------------------------------------------------------
% arguments
%----------------------------------------------------------------------
nargs                       = nargin;
if nargs < 2 || isempty( r ) || isempty( uclu )
    fprintf( '%s: missing arguments (r or uclu)!!\n', upper( mfilename ) )
    return
end
[ map, shankclu ...
    , Detailed, doMUA ]     = ParseArgPairs(...
    { 'map', 'shankclu' ...
    , 'Detailed', 'doMUA' }...
    , { [], [] ...
    , 1, 1 }...
    , varargin{ : } );

r                           = r( : );
sr1                         = size( r{ 1 } );
for i                       = 2 : length( r )
    if ~isequal( sr1, size( r{ i } ) )
        fprintf( '%s: input size (r array) mismatch\n', upper( mfilename ) )
        return
    end
end
nsamples                    = sr1( 1 );
ntrials                     = sr1( 2 );
nrelements                  = length( r );
if ~isempty( uclu ) && length( uclu ) ~= nrelements
    fprintf( '%s: input (clu vs. r) size mismatch\n', upper( mfilename ) )
    return
end

if Detailed && ~isempty( map ) && ~isempty( shankclu ) && size( shankclu, 1 ) == size( map, 1 )
    umap                    = map( ismember( map( :, 1 ), uclu ), 2 : 3 );
    umap                    = [ umap shankclu( ismember( shankclu( :, 1 : 2 ), umap, 'rows' ), 3 ) ];
else
    Detailed                = 0;
end
if Detailed == 0
    doMUA                   = 0;
end

%----------------------------------------------------------------------
% arrange
%----------------------------------------------------------------------

if doMUA == 0                                                               % requires: uclu, umap, r
    
    % preparations for SU spike trains:
    nclu                    = length( uclu );
    tidx                    = ( 1 : ntrials )' * ones( 1, nclu );
    tidx                    = tidx( : );
    uidx                    = ones( ntrials, 1 ) * uclu( : ).';
    uidx                    = uidx( : );
    
    if Detailed
        tmp1                = ones( ntrials, 1 ) * umap( :, 1 ).';
        tmp2                = ones( ntrials, 1 ) * umap( :, 2 ).';
        tmp3                = ones( ntrials, 1 ) * umap( :, 3 ).';
        eshankclu           = [ tmp1( : ) tmp2( : ) tmp3( : ) ];
        tmat                = single( [ eshankclu tidx ] );                 % [ shank; clu; type; trial index ]
    else
        tmat                = single( [ uidx tidx ] );                      % [ cluster number - NOT index; trial index ]
    end
    
    % the SU spike trains:
    st                      = sparse( nsamples, nclu * ntrials );
    for i                   = 1 : nclu
        idx                 = uidx == uclu( i );
        st( :, idx )        = r{ i };
    end
    nsources                = ones( nclu, 1 );
    
else                                                                        % agglomerate; requires: umap, r
    
    % make an identical structure (sparse matrix) in which all same-type/same-shank spikes are pooled:
    shanks                  = umap( :, 1 );
    celltypes               = umap( :, 3 );
    ushanks                 = unique( shanks );
    nshanks                 = length( ushanks );
    ucelltypes              = unique( celltypes );
    ncelltypes              = length( ucelltypes );
    
    % here the pointer includes: shank, cell type, trial
    tidx                    = ( 1 : ntrials )' * ones( 1, nshanks * ncelltypes );
    tidx                    = tidx( : );
    sidx                    = ones( ntrials * ncelltypes, 1 ) * ushanks( : ).';
    sidx                    = sidx( : );
    cidx                    = ones( ntrials, 1 ) * ucelltypes( : ).';
    cidx                    = cidx( : );
    cidx                    = repmat( cidx, [ nshanks 1 ] );
    nans                    = NaN * ones( size( cidx ) );
    tmat                    = single( [ sidx nans cidx tidx ] );            % [ shank; type; trial index ]
    % agglomerate:
    st                      = sparse( nsamples, nshanks * ncelltypes * ntrials );
    nsources                = zeros( nshanks, ncelltypes );
    for i                   = 1 : nshanks
        shanknum            = ushanks( i );
        for j               = 1 : ncelltypes
            ct              = ucelltypes( j );
            sidx            = find( umap( :, 1 ) == shanknum & umap( :, 3 ) == ct );    % source
            tidx            = tmat( :, 1 ) == shanknum & tmat( :, 3 ) == ct;            % target
            nsources( i, j ) = length( sidx );
            for k           = 1 : length( sidx )
                st( :, tidx ) = st( :, tidx ) + r{ sidx( k ) };
            end
        end
    end
    
end

return

% EOF
