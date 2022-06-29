% stmix                 shuffle/jitter spikes in a sparse array
%
% call                  xm = stmix( x, halfwin, nreps )
%
% gets                  x               (m x n) sparse array
%                       halfwin         {[]}; [samples]. empty shuffles
%                       nreps           {1}; if an integer > 1, replicates each column nreps times; then 
%                                           xm is m x nreps*n in blocks of nreps columns
%
% returns               xm              shuffled/jittered arrays
%
% calls                 nothing
%
% see also              streconstruct

% 17-jun-14 ES

% revisions
% 14-oct-19 (1) cleaned up and documented

% to do: interval jitter, randomize seed

function xm = stmix( x, h, p )

xm                      = [];

% arguments
nargs                   = nargin;
if nargs < 1 || isempty( x )
    return
end
if nargs < 2 || isempty( h )
    h                   = [];
end
if nargs < 3 || isempty( p )
    p                   = 1;
end

% intialize
if ~issparse( x )
    return
end
[ m, n ]                = size( x );
nspks                   = full( sum( x ) );
cspks                   = [ 0 cumsum( p * nspks ) ];
trN                     = zeros( sum( nspks ) * p, 1 );

% randomize
if isempty( h )         % shuffle (re-distribute over duration)
    [ ~, rnd ]          = sort( rand( m, n * p ), 1 );
    for i               = 1 : n
        cols            = ( 1 : p ) + ( i - 1 ) * p;
        idx             = ( cspks( i ) + 1  ) : cspks( i + 1 );
        d1              = ones( nspks( i ), 1 ) * ( cols - 1 ) * m;
        tr              = rnd( 1 : nspks( i ), cols ) + d1;
        trN( idx )      = tr( : );
    end
else                    % jitter (shift each spike)
    ntot                = sum( nspks );
    jit                 = round( 2 * h * rand( ntot * p, 1 ) - h );
    for i               = 1 : n
        cols            = ( 1 : p ) + ( i - 1 ) * p;
        idx             = ( cspks( i ) + 1  ) : cspks( i + 1 );
        d1              = ones( nspks( i ), 1 ) * ( cols - 1 ) * m;
        st              = find( x( :, i ) ) * ones( 1, p );
        d               = jit( idx );
        trN( idx )      = st( : ) + d1( : ) + d;
    end
    trN( trN < 1 )      = [];
    trN( trN > m * n * p ) = [];
end

% build arrays
[ ii, jj ]              = ind2sub( [ m n * p ], trN );
xm                      = sparse( ii, jj, 1, m, n * p );

return

% EOF