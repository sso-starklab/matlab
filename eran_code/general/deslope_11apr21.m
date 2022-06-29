% DESLOPE           remove 'slope' from detrended waveform.
%
% call              y = deslope( x, dflag )
%
% gets              x           may be a matrix ( works on columns )
%                   deslope     {0}; 1 also removes the 'dc' s.t. first and
%                                   last samples are at zero
%
% returns           y       desloped columns (first and last samples at same level)
%
% note              this function does not detrend
%
% see also          mydetrend 

% 22-may-03 ES

% revisions
% 11-apr-21 cleaned up

function y = deslope( x, dflag )

nargs                           = nargin;
if nargs < 1 || isempty( x ) || ~ismatrix( x )
    y                           = [];
    return
end
if nargin < 2 || isempty( dflag )
    dflag                       = 0;
end

[ n, m ]                        = size( x );
ridx                            = ( n / 2 : - 1 : - n / 2 + 1 )';
dx                              = x( n, : ) - x( 1, : );
ss                              = dx / ( n - 1 );
idx                             = ridx( :, ones( m, 1 ) );
mult                            = ss( ones( n, 1 ), : );
lt                              = idx .* mult - mult / 2;
y                               = x + lt;
if dflag
    dc                          = ones( n, 1 ) * y( 1, : );
    y                           = y - dc;
end

return

% EOF
