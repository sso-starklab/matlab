% poissonTest       probability of more/less events given an expected
% 
% [ pInc pDec surprise ] = poissonTest( baseRate, inCount, inTime )
%
% baseRate          rate, e.g. [spikes/s]
% inCount           count, e.g. [spikes]
% inTime            time, e.g. [s]
%
% pInc              probability to see the observed (or more), given the expected
% pDec              probability to see the observed (or less), given the expected
% suprise           composite measure (positive if the observed is higher
%                       than expected, negative if less)
% 
% note:
% the baseRate is typically the total spike count during a baseline period
% (could be many non-overlapping periods pooled together), divided by the
% total time in those periods
%
% assumptions:
% that everything is stationary, i.e. that the statistics during baseline
% and the tested periods are the same. Actually, what this routine does is
% test this assumption (can be thought of as the null hypothesis to be tested)

% 20-may-14 ES

% revisions
% 11-sep-19 cleaned up
% 10-dec-19 prevented zero values for probabilities

function [ pInc, pDec, surp ] = poissonTest( baseRate, inCount, inTime )

if nargin < 3 || ~isequal( size( baseRate ), size( inCount ), size( inTime ) )
    return
end
siz                 = size( baseRate );
baseRate            = baseRate( : );
inCount             = inCount( : );
inTime              = inTime( : );

lambdas             = baseRate .* inTime;
pInc                = 1 - poisscdf( inCount - 1, lambdas );
pDec                = poisscdf( inCount, lambdas );
pInc( pInc < eps )  = eps;
pDec( pDec < eps )  = eps;
surp                = log10( pDec ./ pInc );

pInc                = reshape( pInc, siz );
pDec                = reshape( pDec, siz );
surp                = reshape( surp, siz );

return

% EOF