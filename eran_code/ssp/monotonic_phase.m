% monotonic_phase      force a phase vector to be monotonically increasing when unwrapped
%
% CALL                 [ phs, periods ]   = monotonic_phase( phs )
%
%
% GETS                 phs      phs vector (matrices are vectorized)
%
%
% RETURNS
%                      phs      the same vector, with fixed phase in regions of decrease in phase
%                      periods  ranges of indices in which there was a decrease in phs (and were changed)
%
% NOTE                  cycle allocation is maintained on average
%                       phase is monotonically increasing

% written by           ES+HS  29-Jun-21
% modified             HS     30-Jun-21  (1) added input checking
%                                        (2) wrote help
% 01-jul-21             removed tol

% to do:
% - if phs was originally between -pi:pi, return to that range at the end
%   of the function (simply wrapToPi makes the vector non-monotonically
%   increasing)

function [ phs, periods ]   = monotonic_phase( phs )

if nargin < 1 || isempty(phs)
    error('No input');
end

if ( min( phs ) < -pi ) || ( max( phs ) > ( 2 * pi ) )
    error( 'phs must be a circular input' );
end

% unwrap the phase
phs                         = phs( : ); 
%phs                         = wrapToPi( phs );
phs                         = wrapTo2Pi( phs );
cphs                        = unwrap( phs );
lvec                        = length( phs );

phs0                        = phs;
cphs0                       = cphs;

% detect segments of decreasing phase
periods                     = parse( find( diff( cphs ) < 0 ) );
nperiods                    = size( periods, 1 );
if nperiods == 0
    return
end

% go over segments and enforce each case to be flat
for i                       = 1 : nperiods
    st                      = periods( i, 1 );
    val                     = cphs( st );
    et                      = st + 1;
    while ( cphs( et ) < val ) && ( et < lvec )
        et                  = et + 1;
    end
    if et == lvec
        idx                 = st : et;
    else
        idx                 = st : ( et - 1 );
    end
    cphs( idx )            	= val;
end

% convert back to phases
phs                         = mod( cphs, 2 * pi );
%phs                         = wrapToPi( phs );

% problem:
% x = rand( 100, 1 ) * 10; isequal( unwrap( mod( x, 2 * pi ) ), x )
% isequal( unwrap( mod( 0.1 : 10, 2 * pi ) ), 0.1 : 10 ) 
% isequal( unwrap( mod( cphs, 2 * pi ) ), cphs )

% cyc0                        = ceil( unwrap( wrapTo2Pi( phs0 ) ) / ( 2 * pi ) );
% cyc                         = ceil( unwrap( wrapTo2Pi( phs ) ) / ( 2 * pi ) );
% sum( cyc0 ~= cyc )
% sum( ceil( cphs0 / ( 2 * pi ) ) ~= ceil( cphs / ( 2 * pi ) ) )
% sum( ceil( cphs0 / ( 2 * pi ) ) ~= cyc0 )
% sum( ceil( cphs / ( 2 * pi ) ) ~= cyc )

% make sure the function worked properly 
test                        = sum( diff( unwrap( phs ) ) < 0 );
if test > 0
    fprintf('Monotonic phase did not fix %d instences of non-monotony \n', test);
    keyboard
end

return

% EOF

figure
subplot( 2, 1, 1 ), plot( cphs0( 1 : 1250 ), 'b' ), hold on, plot( cphs( 1 : 1250 ), 'r' ), title( sum( diff( cphs ) < 0 ) )
subplot( 2, 1, 2 ), plot( phs0( 1 : 1250 ), 'b' ), hold on, plot( phs( 1 : 1250 ), 'r' ), title( sum( diff( unwrap( phs ) ) < 0 ) )

% call example:
sum( diff( unwrap( phs ) ) < 0 )
[ mphs, mperiods ] = monotonic_phase( phs );
sum( diff( mperiods, [], 2 ) + 1 )
sum( diff( unwrap( mphs ) ) < 0 )




mphs       = monotonic_phase( y );

    
block_size = 1e5;
lphs       = length(y);
blocks     = makeblocks( lphs, block_size, 0 );
nblocks    = size(blocks,1);

figure
for i = 1 : nblocks
    
    idx    = blocks(i,1):blocks(i,2);
    clf
    plot(y(idx),'--r'); hold on; plot(mphs(idx),'--b');
    title(sprintf('block %d', i));
    pause
    
end



