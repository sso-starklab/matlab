% GTEST             tests of goodness-of-fit and of independence.
%
% call              [ P, G ] = GTEST( OBS, TEST, CFLAG, MODE )
%
% gets              OBS             matrix
%                   TEST            {'lr'} (G test) or 'chi2' (Chi2 test)
%                   CFLAG           correction flag {1}
%                   MODE            independence test {'ind'}; if a
%                                       goodness-of-fit ('gof') test is desired, OBS
%                                       should be a 2 column matrix: [ observed expected ]
%
% returns           P               p-value of H0: independence
%                   G               G/Chi2 statistic
%
% notes
%                   1. William's correction for goodness-of-fit LR test is for an extrinsic hypothesis only
%                   2. Do not use corrections for additive (anova-like) uses of LR test
%                   3. Chi2 includes Yate's continuity correction (only for 2 x 2, n <= 200)
%                   4. n < 25 is done exactly only for 2 x 2
%                           goodness-of-fit (binomial; otherwise - should be multi-nomial)
%                   5. Model design is type I or II (all marginals free or
%                           either rows or columns fixed). If Model III
%                           (both marginals fixed), Fisher's test should be used (not included here)

% 15-jul-04 ES

% reference: Sokal & Rohlf 2001 ch. 17.

function [ p, g ] = gtest( obs, test, cflag, mode )

nargs = nargin;
if nargs < 2 | isempty( test ), test = 'g'; end
if nargs < 3 | isempty( cflag ), cflag = 1; end
if nargs < 4 | isempty( mode ), mode = 'ind'; end

p = NaN;
g = NaN;

% expected frequencies

[ r c ] = size( obs );
cols = sum( obs, 1 );
rows = sum( obs, 2 );
n = sum( rows );
switch lower( mode )
    case 'ind'
        expc = rows * cols / n;
    case 'gof'
        if c ~= 2% | cols( 1 ) ~= cols( 2 )
            error( 'goodness-of-fit option requires two columns (equal counts)' )
        end
        P = obs( :, 2 ) / cols( 2 );
        expc = P * cols( 1 );
        obs = obs( :, 1 );
        n = cols( 1 );
    otherwise
        error( 'unknown MODE requested' )        
end

% special cases

if n <= 25 & strcmp( lower( mode ), 'gof' ) & r == 2            % exact (n<=25, 2 class gof test)
    p=1-max( binomial( obs( 1 ), n, P( 1 ) ), binomial( obs( 2 ), n, P( 2 ) ) );
%    p = binomial( obs( 1 ), n, P( 1 ) ) + 1 - binomial( obs( 2 ) - 1, n, P( 1 ) );
    return
end

if sum( sum( expc < 5 ) ) / r / c > 0.2
    return
end

% minf = min( min( expc ) );                                      % Biometry recommendation: should be grouped
% if minf < 3 | ( r * c < 5 & minf < 5 )
%     p = NaN;
%     return
% end

% compute the statistic

dof = ( r - 1 ) * ( c - 1 );

switch lower( test )
    case { 'chi2', 'chi2_test' }
        if cflag & r == 2 & c == 2 & n <= 200                   % Yate's corrections to 2 x 2 tables
            switch lower( mode )
                case 'ind'                                      % p. 737
                    s = ( abs( det( obs ) ) - n / 2 ) .^ 2 * n ./ prod( rows ) ./ prod( cols );
                case 'gof'                                      % p. 704
                    s = ( abs( obs - expc ) - 0.5  ) .^ 2 ./ expc;
            end
        else                                                    % p. 697
            s = ( ( expc - obs ) .^ 2 ) ./ expc;
        end
        g = sum( sum( s ) );
    case { 'g', 'lr' }                                          % p. 692
        idx = obs > 0;
        s = obs( idx ) .* log( obs( idx ) ./ expc( idx ) );
        g = 2 * sum( s );
        if cflag                                                % William's corrections
            switch lower( mode )
                case 'ind'                                      % p. 738
                    q = 1 + ( n * sum( 1 ./ rows ) - 1 ) * ( n * sum( 1 ./ cols ) - 1 ) / ( 6 * n * dof );
                case 'gof'                                      % p.698
                    q = 1 + ( r + 1 ) / ( 6 * n );              % extrinsic hypothesis only
            end
            g = g / q;
        end
    otherwise
        error( 'unknown TEST requested' )
end

% evaluate probability

p = 1 - chi2cdf( g, dof );

return

% % example of p. 731
% mat = [ 12 22; 16 50 ];
% [ p g ] = gtest( mat, 'g', 0, 'ind' )
% [ p g ] = gtest( mat, 'g', 1, 'ind' )
% call_multinomial = inline( 'multinomial( mat( : ), [], mat( : ) / sum( sum( mat ) )  )', 'mat' );

% example of p. 717
mat = [ 83 47; 77 43; 110 96; 92 58; 51 31; 48 61; 70 42; 85 66];
[ r c ] = size( mat );
for i = 1 : r                                       % go over data
    [ p g( i ) ] = gtest( [ mat( i, : ); [ 9 7 ] ]', 'g', 0, 'gof' ); 
end
totalG = sum( g );
[ p pooledG ] = gtest( [ sum( mat, 1 ); [ 9 7 ] ]', 'g', 0, 'gof' ); 
heteroG = totalG - pooledG;                         % [ ign heteroG ] = gtest( mat, 'g', 0, 'ind' ); % mathematically the same
1 - chi2cdf( g, c - 1 )                             % each individual (here, pair-wise) comparison
1 - chi2cdf( pooledG, c - 1 )                       % all data together
1 - chi2cdf( heteroG, ( r - 1 ) * ( c - 1 ) )       % interactions
1 - chi2cdf( totalG, r * ( c - 1 ) )                % due to pooled / hetero