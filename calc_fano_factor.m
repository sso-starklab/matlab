% calc_fano_factor          from multi-unit spike trains and window sizes
%
% call                      [ fanomat, myu, sd2, nwindows, nspikes ] = calc_fano_factor( res, clu )
%
% gets                      res         n-element vector list of spike times [samples]
%                           clu         n-element vector of spike identity
%
% optional arguments (given as name/value pairs)
%
%                           winsizes    {2^(4,5,...,14)};   m-element vector of window sizes [samples]
%                           drange      {[1 max( res ) ]};  time range of data
%                           verbose     { 1 };              screen output of progress
%                   
% returns                   fanomat     nclu x m matrix of Fano factors (:= variance/mean)
%                           myu         matrix of mean spike counts
%                           sd2         matrix of variance of spike counts
%                           nwindows    m-element vector, number of whole windows in data (same for all units)
%                           nspikes     nclu-element vector, number of spikes per unit (same for all windows)
%
% calls                     calc_com, uhist

% 22-dec-19 SSo + ES

% revisions
% 24-dec-19 completed

% to do:
% implement with diff; should improve from O(nlogn) to O(n)

function [ fanomat, myu, sd2, nwindows, nspikes ] = calc_fano_factor( res, clu, varargin )

% handle input
nargs                   = nargin;
if nargs < 2 || isempty( res ) || isempty( clu )
    return
end
[ winsizes, drange, verbose ] = ParseArgPairs(...
    { 'winsizes', 'drange', 'verbose' }...
    , { 2 .^ ( 4 : 14 ), [], 1 } ...
    , varargin{ : } );

res                     = res( : );
clu                     = clu( : );
if ~isequal( size( res ), size( clu ) )
    error( 'input size mismatch' )
end
if ~isequal( res, round( res ) )
    error( 'res must be in whole integers (samples)' )
end
if ~isequal( clu, round( clu ) )
    error( 'clu must be in whole integers (categories)' )
end
if ~isequal( winsizes, round( winsizes ) )
    error( 'winsizes must be in whole integers (samples)' )
end
if isempty( drange )
    drange              = [ 1 max( res ) ]; % this is an assumption
end

% initialize output
nwins                   = length( winsizes );
uclu                    = unique( clu );
nclu                    = length( uclu );
nspikes                 = NaN( nclu, 1 );
nwindows                = floor( diff( drange ) ./ winsizes );
fanomat                 = NaN( nclu, nwins );
myu                     = fanomat;
sd2                     = fanomat;

% go over units
for i                   = 1 : nclu
    
    % get the spikes for the unit
    cluI                = uclu( i );
    idx                 = clu == cluI;
    resI                = res( idx );
    nspikes( i )        = sum( idx );
    if verbose
        t0 = clock;
        fprintf( 1, 'clu %d', cluI )
    end
    
    % go over windows
    for wi              = 1 : nwins
        ws              = winsizes( wi );
        
        if verbose
            fprintf( 1, '.' )
        end
        % determine during which window each spike occurred
        nwin            = floor( ( diff( drange, [], 2 ) + 1 ) / ws );
        winI            = ceil( ( resI - drange( 1 ) + 1 )  / ws );
        
        % remove data from incomplete windows
        ridx            = winI > nwin;
        winI( ridx )    = [];
        nwin            = nwin - sum( ridx );
        
        % count number of spikes in each populated window
        cpwin           = uhist( winI );                                    % this is OK but may be improved with diff
        zwin            = nwin - length( cpwin );                           % this is the number of unpopulated windows
        [ weights, values ] = uhist( cpwin );                               % histogram of spike counts
        [ m, s ]        = calc_com( [ 0 values ]', [ zwin weights ]', 1 );  % mean and SD of the histogram
        myu( i, wi )    = m;
        sd2( i, wi )    = s.^2;
        
    end
    
    if verbose
        fprintf( 1, '%0.3g s\n', etime( clock, t0 ) )
    end
    
end

% compute the FF
fanomat                 = sd2 ./ myu;

return

% EOF

