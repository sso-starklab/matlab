% cch_mat           make CCH matrix (one row for each relevant spike )
%
% call              [ cchmat, ures ] = cch_mat( clu, res, BinSize, nBins, pairs, uClu )
%
% gets              clu, res, BinSize, nBins    input to CCG
%                   pairs                       output of CCG
%                   uClu                        2-element vector 
%
% returns           cchmat                      n x ( 2 * nBins + 1 )
%                   ures                        unique (used) spike times (n) of uClu(1)
%
% calls             nothing
%
% note              if the BinSize is 1, cchmat can be logical
%                   if BinSize is up to 255, can be uint8
%                   if BinSize is up to 65535, can be uint16
% 
% see also          calc_asg

% 05-oct-20 ES

% revisions
% 08-oct-20 fix instead of ceil in while loop
% 17-oct-20 error checking added

function [ cchmat, ures1 ] = cch_mat( clu, res, BinSize, nBins, pairs, uClu )

% argument handling
nargs                           = nargin;
if nargs < 6
    error( 'missing arguments' )
end
if ~isa( clu, 'numeric' ) || ~isa( res, 'numeric' ) || ~isvector( clu ) || ~isvector( res )
    error( 'clu, res: must be numeric vector arrays' )
end
if length( clu ) ~= length( res )
    error( 'clu/res: input size mismatch' )
end
if any( clu < 0 ) || sum( clu ~= round( clu ) ) 
    error( 'clu: must be an array of non-negative integers' )
end
if any( res < 0 ) || sum( res ~= round( res ) ) 
    error( 'res: must be an array of non-negative integers' )
end
if ~isa( BinSize, 'numeric' ) || length( BinSize ) ~= 1 || BinSize < 0 || BinSize ~= round( BinSize )
    error( 'BinSize: must be a non-negative integer (samples)' )
end
if ~isa( nBins, 'numeric' ) || length( nBins ) ~= 1 || nBins < 0 || nBins ~= round( nBins )
    error( 'nBins: must be a non-negative integer (bins, half-width of CCH)' )
end
if ~isa( pairs, 'numeric' ) || ~ismatrix( pairs ) || size( pairs, 2 ) ~= 2
    error( 'pairs must be a 2-column matrix' )
end
if ~isa( uClu, 'numeric' ) || length( uClu ) ~= 2
    error( 'uClu: input size mismatch' )
end

% get the relevant information
pidx                            = clu( pairs( :, 1 ) ) == uClu( 1 ) & clu( pairs( :, 2 ) ) == uClu( 2 );
npairs                          = sum( pidx );
res1                            = res( pairs( pidx, 1 ) );
res2                            = res( pairs( pidx, 2 ) );
ures1                           = unique( res1 );
nures1                          = length( ures1 );

% initialize and populate the matrix
cchmat                          = zeros( nures1, 2 * nBins + 1, 'single' );
r                               = 0;
for i                           = 1 : nures1
    resI                        = ures1( i );
    r                           = r + 1;
    while r < npairs
        c                       = fix( ( res2( r ) - res1( r ) ) / BinSize ) + nBins + 1;
        cchmat( i, c )          = cchmat( i, c ) + 1;
        if res1( r + 1 ) == resI
            r                   = r + 1;
        else
            break
        end
    end
end

return

% EOF
