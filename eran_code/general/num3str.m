% NUM3STR       clips to range, converts to str, and formats
%
% STR = NUM3STR( NUM, NCHARS )
%
% examples:
% num3str( 3.4, 3 ) -> 003
% num3str( 3.4, 4 ) -> 0003
% num3str( 34, 3 ) -> 034
% num3str( 34, 2 ) -> 34
% num3str( 34, 1 ) -> 9

% 07-apr-11 ES

function str = num3str( num, nchars )

if nargin < 2, nchars = 3; end

num = num( : );
num = min( num, 10^nchars - 1 );
num = max( num, 0 );
num = round( num );
for i = 1 : length( num )
    pad = repmat( '0', [ 1 nchars - ceil( log10( num( i ) + 1 ) ) ] );
    str( i, : ) = sprintf( '%s%d', pad, num( i ) );
end

return