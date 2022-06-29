% BINCOLS           sum matrix columns in bins.
%
% call              B = BINCOLS( M, BS, BO, MODE )
%
% gets              M       input matrix
%                   BS      bin size
%                   BO      bin overlap (optional)
%                   MODE    {'sum'}, 'mean', or 'rms'
%
% returns           B       binned matrix

% 11-apr-03 ES

% revisions
% 11-may-04 BO handled
% 21-jan-06 mode added

% if binning is begun on first ms, some structure may be lost.
% assume for example binsize of 2:
%       0 1 1 0 -> 1 1
%       1 1 0 0 -> 2 0
% etc.

function B = bincols( M, BS, BO, mode )

nargs = nargin;
if nargs < 2, error( '2 arguments' ), end                           % check
if nargs < 3 || isempty( BO ), BO = 0; end                          % default
if nargs < 4 || isempty( mode ), mode = 'sum'; end                  % what to do
[ r c ] = size( M );                                                % size
d = BS - BO;                                                        % actual bin size
nr = floor( r / d );                                                % prepare trim
M( nr * d + 1 : r, : ) = [];                                        % trim
m = reshape( M, d, nr * c );                                        % prepare bin
switch mode
    case 'sum'
        b = sum( m, 1 );                                                    % bin
    case 'mean'
        b = mean( m, 1 );
    case 'rms'
        b = sqrt( mean( m.^2, 1 ) );
end
B = reshape( b, nr, c );                                            % back
if BO, 
    win = [ ones( round( BS / d ), 1 ); 0 ];
    win = win / sum( win );                                         
    B = firfilt( full( B ), win );                                  % approximate
end

return