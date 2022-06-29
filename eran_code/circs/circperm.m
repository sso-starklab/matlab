% circperm          two-sample circular permutation test for equality of means/dispersion/distributions
%
% call              [ p, r0, rm ] = circperm( p1, p2, method, nreps )
% 
% gets              p1, p2          samples [rad]
%                   method          {0}; 0/1/2 correspond to all/mean/dispersion
%                   nreps           {1000}, number of permutation
% 
% returns           p               p value
%                   r0              test statistic
%                   rm              resampled statistic
%
% calls             car2pol, circ_mean, mixmat, pol2car
%
% reference         Stark and Abeles 2005
%
% see also          utest, wheeler

% 18-jun-13 ES

% revisions
% 12-dec-19 cleaned up

function [ p, r0, rm ] = circperm( p1, p2, method, nreps )

% intialize output
p                       = [];
r0                      = [];
rm                      = [];

% handle arguments
nargs                   = nargin;
if nargs < 1 || nargs < 2 
    return
end
if nargs < 3 || isempty( method )
    method              = 0;
end
if numel( method ) ~= 1 || ~ismember( method, 0  : 2 ) 
    method              = 0;
end
if nargs < 4 || isempty( nreps )
    nreps               = 1000;
end

% prepare data
p1                      = p1( : );
p1( isnan( p1 ) )       = [];
p2                      = p2( : );
p2( isnan( p2 ) )       = [];
n1                      = length( p1 ); 
n2                      = length( p2 );

% compute test statistic
[ ph1, r1 ]             = circ_mean( p1 );
[ ph2, r2 ]             = circ_mean( p2 );
switch method
    case 0 % any difference
        [ x1, y1 ]      = pol2car( r1, ph1 ); 
        [ x2, y2 ]      = pol2car( r2, ph2 ); 
        r0              = car2pol( x1 + x2, y1 + y2 );
    case 1 % mean directions
        r0              = pi - abs( pi - abs( ph1 - ph2 ) );
    case 2 % spread
        r0              = abs( r1 - r2 );
end

% resample and recompute test statistic + p-value
p12                     = [ p1; p2 ];
p12hat                  = mixmat( p12 * ones( 1, nreps ), 1, 1 );
p1hat                   = p12hat( 1 : n1, : );
p2hat                   = p12hat( n1 + 2 : n1 + n2, : );

[ ph1m, r1m ]           = circ_mean( p1hat );
[ ph2m, r2m ]           = circ_mean( p2hat );
switch method
    case 0 % any difference
        [ x1, y1 ]      = pol2car( r1m, ph1m ); 
        [ x2, y2 ]      = pol2car( r2m, ph2m ); 
        rm              = car2pol( x1 + x2, y1 + y2 );
        p               = ( sum( r0 > rm ) + 1 ) / ( nreps + 1 );
    case 1 % mean directions
        rm              = pi - abs( pi - abs( ph1m - ph2m ) );
        p               = ( sum( r0 <= rm ) + 1 ) / ( nreps + 1 );
    case 2 % spread
        rm              = abs( r1m - r2m );
        p               = ( sum( r0 > rm ) + 1 ) / ( nreps + 1 );

end

return

% EOF
