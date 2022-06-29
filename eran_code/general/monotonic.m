% monotonic         make a signal monotonically-increasing
%
% call              y = monotonic( x, step )
% 
% gets              x           vector (if matrix, vectorizes)
%                   step        {1}
%
% returns           y
% 
% calls             nothing

% 28-jul-13 ES

% revisions
% 23-mar-21 cleaned up

function y = monotonic( x, step )

y                               = [];
nargs                           = nargin;
if nargs < 1 || isempty( x )
    return;
end
if nargs < 2 || isempty( step )
    step                        = 1; 
end

y                               = x( : );
dx                              = [ 0; diff( y ) ];
neg                             = find( dx < 0 );
for i                           = 1 : length( neg )
    xidx                        = neg( i ) : length( y ); 
    y( xidx )                   = y( xidx ) + step;
end

return

% EOF