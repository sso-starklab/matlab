% calc_index            (a-b)./(a+b)
%
% call                  idx = calc_index( x, y )
%
% gets                  x,y     matrices of same size
%                       x,[]    x can be a two-column matrix
%
% returns               (x-y)./(x+y), limited to the range -1 to 1 if all
%                                     values are positive
%
% note:
% tanh( log10( a ./ b ) ) ~= (a-b)/(a+b)

% 22-dec-13 ES

function idx = calc_index( x, y )

if nargin < 2 || isempty( y )
    if size( x, 2 ) == 2
        y = x( :, 2 );
        x = x( :, 1 );
    else
        idx = [];
        return
    end
end
idx = ( x - y ) ./ ( x + y );

return

% EOF


