% RANKCOLS          rank matrix column-wise.
%
% call              [ Y T ] = RANKCOLS( X )
%
% gets              X       data matrix
%
% returns           Y       matrix of ranks (of each column, ties handled)
%                   T       correction for tied variates
%
% calls             nothing.

% reference: Sokal & Rohlf 2001, p.426

% 12-jul-04 ES

% revisions
% 08-may-13 NaNs
% 18-aug-19 cleaned up

function [ y, T ] = rankcols( x )

[ m, n ]            = size( x );
if m == 1
    x               = x';
    m               = n;
    n               = 1;
end
nans                = isnan( x );
ridx                = m : -1 : 1;
cidx                = 1 : n;
[ x, idx ]          = sort( [ x x( ridx, : ) ] );
[ ~, ranks ]        = sort( idx );
y                   = ( ranks( :, cidx ) + ranks( ridx, cidx + n ) ) / 2;
y( nans )           = NaN;

if nargout > 1
    T               = zeros( 1, n );
    for i           = 1 : n
        t           = diff( parse( find( diff( x( :, i ) ) == 0 ) ), [], 2 ) + 2; 
        T( i )      = sum( ( t - 1 ) .* t .* ( t + 1 ) );
    end
end

return

% test:

n = 100;
m = 5;
x = round( rand( n, m ) * 20 ); 
for i = 1 : m, r1( :, i ) = tiedrank( x( :, i ) )';  end
r2 = rankcols( x );
for i = 1 : m, if ~isequal( r1( :, i ), r2( :, i ) ), disp( 'mismatch' ), end, end

% ex. from p. 428:
y1 = [ 104 109 112 114 116 118 118 119 121 123 125 126 126 128 128 128 ];
y2 = [ 100 105 107 107 108 110 116 120 121 123 ];

% example from p. 424:
y = [ 56 * ones( 3, 1 )
    57 * ones( 7, 1 )
    58 * ones( 7, 1 )
    59 * ones( 4, 1 )
    60 * ones( 4, 1 )
    61 * ones( 4, 1 )
    62 * ones( 4, 1 )
    63 * ones( 1, 1 )
    64 * ones( 1, 1 )
    65 * ones( 4, 1 )    
    66 * ones( 1, 1 )
    67 * ones( 4, 1 )    
    68 * ones( 1, 1 )
    70 * ones( 1, 1 )    
    71 * ones( 1, 1 )    
    75 * ones( 2, 1 )    
    76 * ones( 1, 1 )    
];