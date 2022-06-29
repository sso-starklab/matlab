%   RAY_TEST        Rayleigh's LR test of uniformity for directional data.
%
%       call:   p = RAY_TEST( theta, f )
%       does:   Test the null hypothesis (HO):
%                                               f(theta) = 1/(2*pi)
%               against the alternative (H1):
%                                               f(theta) = g(theta),
%               where g corresponds to von Mises' pdf.
%               Based on the Pearson approximation for the nR^2 statistic
%       output: the p value
%
%       See also RAY_CDF

% directional statistics package
% Dec-2001 ES

% revisions
% 09-apr-13 modified for matrix use
% 31-dec-20 (1) support of NaNs
%           (2) cleaned up

function p = ray_test( theta, f )

nargs                           = nargin;
if nargs < 2 || isempty( f )
    f                           = ones( size( theta ) ); 
end

% check input size
if length( theta ) < 5
    p                           = 1;    %   cannot estimate, likely to be non-uniform
    return
end

[ tm, tn ]                      = size( theta );
[ fm, fn ]                      = size( f );
if tm == fm
    if tn == fn
        % ok
    elseif tn == 1
        theta                   = theta * ones( 1, fn );
    elseif fn == 1
        f                       = f * ones( 1, tn );
    end
elseif tn == 1 && tm == fn
    theta                       = theta * ones( 1, fm );
    f                           = f';
elseif fn == 1 && fm == tn
    f                           = f * ones( 1, tm );
    theta                       = theta';
end
if ~isequal( size( theta ), size( f ) )
    error( 'input size mismatch' )
end
        
% calculate the statistic R
n                               = nansum(f,1);
x                               = f.*cos(theta);
y                               = f.*sin(theta);

C                               = nansum(x,1)./n;
S                               = nansum(y,1)./n;   
R                               = (C.^2 + S.^2).^0.5;

% and the p value (~chi2 w/ 2 dof)
% p = 1 - chi2cdf( 2 * n .* R.^2, 2 );
% based on chi2cdf - valid for large n only

% and the p
p                               = ray_cdf( R, n );
p( p < 0 )                      = 0;
p( p > 1 )                      = 1;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RAY_CDF         find the cumulative distribution for the Rayleigh LR test.
%
%       call:   PR = RAY_CDF(R,N).
%       does:   gives the Pearson approximation for the nR^2 statistic,
%               namely, Pr(nR^2 >= K).
%
%               the return value indicates the probability of
%               R to be larger than or equal than K under the
%               null hypothesis that the data is uniformly distributed.
%               a large return value means the data is likely to
%               originate in a uniform distribution.

% directional statistics package
% Nov-2001 ES

function pr = ray_cdf( R, n )

k                               = n.*R.^2;
pr                              = exp(-k).*( 1 + (2*k - k.^2)./(4*n) - (24*k - 132*k.^2 +  76*k.^3 - 9*k.^4)./(288*n.^2) );

return

% EOF
