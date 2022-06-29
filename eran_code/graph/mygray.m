% MYGRAY            modified gray colormap (logarithmically spaced).
%
% call              MAP = MYGRAY( M )
%
% gets              M       colormap length
%
% returns           MAP     the colormap ( M x 3 )
%
% calls             nothing
%
% see also          MYJET

% 08-may-04 ES

function map = mygray( m )

if nargin < 1 | isempty( m ), m = 64; end

map = flipud( [ 1 - ( logspace( 1, 0, m ) - 1 ) / 10 ]' * ones( 1, 3 ) );

return