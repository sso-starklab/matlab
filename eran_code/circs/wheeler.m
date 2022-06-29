% WHEELER       Wheeler's 2-sample test for directional data.
%
% call          [ H0, p_value ] = wheeler( theta1, theta2, alpha, graph )
%
% gets          theta1          vector of angles, sampled from one distribution [rad]
%               theta2          angles from another distribution [rad]
%               alpha           {0.05} significance level
%               graph           {''}, optional presentation: 'compass',
%                                   compass plot; 'rose', rose plot;
%                                   'plot', scatter gram
%
% does          Test the null hypothesis (HO):
%                                               f(theta1) = f(theta2)
%               against the alternative (H1):
%                                               f(theta1) ~= f(theta2)
%
%               effectively, this tests the circular medians for equality,
%               akin to Mann-Whiney's U-test for linear data
%
% procedure     Observations are linearly rankled and then replaced
%               by uniformly distributed scores.
%
% returns       H0              1 if accepted, 0 if rejected
%               p_value         the p value (-1 if cannot compare, this
%                               will happen if any of the distributions
%                               are sampled sparsely with <4 samples)
%
% calls         nothing
% 
% reference     Mardia, 1972
% 
% see also      circperm, utest

% directional statistics package
% Dec-2001 ES

% note: need to untie matches

% revisions
% 12-dec-19 cleaned up

function [ H0, p_value ] = wheeler( a, b, alpha, graph )

%---------------------------------------------------------------
% handle arguments
%---------------------------------------------------------------
% initialize output
H0                      = 0;
p_value                 = NaN;

nargs = nargin;
if nargs < 2 || isempty( a ) || isempty( b )
    return
end

% check input size and make sure they are row vectors
[ s1, s2 ]              = size( a );
[ s3, s4 ]              = size( b );
if ( ( s1 ~= 1 ) && ( s2 ~= 1 ) ) || ( ( s3 ~= 1 ) && ( s4 ~= 1 ) )
    error('inputs should be vectorial')
end
if s1 > s2
    a                   = a';
end
if s3 > s4
    b                   = b';
end

% check requested alpha
if nargs < 3 || isempty( alpha )
    alpha               = 0.05;
end
if ( alpha <= 0) || ( alpha >= 1 )
    error('alpha should be between 0 and 1')
end

% check graphs
if nargs < 4 || isempty( graph )
    graph               = '';
end
if ~ischar( graph ) || ~ismember( graph, { 'compass', 'plot', 'rose' } )
    graph               = '';
end

%---------------------------------------------------------------
% estimate p-value
%---------------------------------------------------------------
% prepare input
a                       = sort( a );
b                       = sort( b );
na                      = length( a );
nb                      = length( b );

% cannot estimate, likely to be non-uniform
if na < 4 || nb < 4
    p_value             = -1;
    return
end

% arrange the signed data points
comb                    = sortrows( [ a b; ones( 1, na ) 2 * ones( 1, nb ) ]' );
beta                    = find( comb( :, 2 ) == 1 ) * 2 * pi / ( na + nb );

% compute the test statistic R
xa                      = cos( beta );
ya                      = sin( beta );
C                       = sum( xa );
S                       = sum( ya );
R                       = 2 * ( na + nb - 1 ) * ( C^2 + S^2 ) / ( na * nb );

% it is distributed as chi2, so
p_value                 = 1 - chi2cdf( R, 2 );

% compare
if p_value > alpha
    H0                  = 1;
else
    H0                  = 0;
end

%---------------------------------------------------------------
% graphics
%---------------------------------------------------------------

if ~isempty( graph )
    beta_b              = find( comb( :, 2 ) == 2 ) * 2 * pi / ( na + nb );
    xb                  = cos( beta_b );
    yb                  = sin( beta_b );
    tbuf                = sprintf( 'Wheeler uniform score test; p value = %g\nn1 = %g, n2 = %g'...
        , p_value, na, nb);
    
    figure
    switch graph
        case 'compass'
            compass( xa, ya, 'b' )
            hold on
            compass( xb, yb, 'r' )
        case 'plot'
            plot( xa, ya, '.b', xb, yb, '.r' )
            xlim( [ -1 1 ] )
            ylim( [ -1 1 ] )
        case 'rose'
            nbins = 20;
            rh = rose( beta, nbins );
            set( rh, 'color', [ 0 0 1 ] );
            hold on
            rh = rose( beta_b, nbins );
            set( rh, 'color', [ 1 0 0 ] );
    end
    hold off
    title( tbuf )
    
end

return

% EOF


% a better approximation - f distribution
% n = na + nb
% fR = R / (n - 1 - R)
% v1 = 1 + ( (n*(n+1)-6*na*nb) / (n*(na-1)*(nb-1)) )
% v2 = (n-3)/round(v1)
% fcdf( fR, round(v2), round(v1) )

% R2 = (C^2+S^2)
% R2 = R*n1*(n-n1)/2/(n-1)