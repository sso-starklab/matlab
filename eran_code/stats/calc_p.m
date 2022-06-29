% CALC_P            calculate p value from empirical sample.
%
% call              P = CALC_P( X, Y, RANKF )
%
% gets              X       sample; does not have to be sorted or a cdf.
%                   Y       scalar
%                   RANKF   {0} - computes probabilities
%                           set to  1   calculate rank instead
%                                   2   rank of 2-sided test
%
% returns           P       fraction of sample larger than X
%                           if RANKF is on, fraction of sample smaller or equal
%                           to X.
%
% calls             nothing

% 25-dec-03 ES

% revisions
% 22-mar-04 NaNs in X ignored
% 12-jan-05 2-tailed p-value added
% 12-sep-19 cleaned up
% 19-dec-19 modified 2-tailed p-value to account for equality

function p      = calc_p( x, y, rankf )

nargs           = nargin;
if nargs < 2
    error( '2 arguments' )
end
if nargs < 3
    rankf       = 0;
end
if  numel( y )  ~= 1
    error( '2nd argument must be a scalar' )
end

x               = x( : );
xi              = ~isnan( x );

if rankf == 1       % rank
    p           = ( sum( y > x( xi ) ) + 1 ) / ( sum( xi ) + 1 );
elseif rankf == 2   % 2-tailed test, rank
    p           = sum( y >= x( xi ) ) / sum( xi );
else                % p-value
    p           = ( sum( y <= x( xi ) ) + 1 ) / ( sum( xi ) + 1 );
end

return

% EOF

