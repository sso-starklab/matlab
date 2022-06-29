% circ_med              circular median - an angle on the diameter that divides a set of points into two ~equal sets
%
% call                  [ med, counts ] = circ_med( x )
%
% gets                  x           vector of n points (or matrix of multiple columns)
%                                   may include NaNs (ignored)
%
% returns               med         angle
%                       counts      number of points on each side of the circular median
%
% notes                 (1) there may be multiple valid circular medians. 
%                           to see this, simply consider an even number of
%                           equally-spaced samples on the unit circle
%                       (2) even for an odd number of samples, the circular
%                           median may be distinct from all samples. to see
%                           this, simply consider four samples - one in
%                           every quarter, all close to the horizontal
%
% calls                 nothing

% 03-apr-20 ES

%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% For a set x of n samples on the line, we define the linear median as the
% "half-point", the sample which is larger than half of the points and
% smaller than half of the points (for an odd n). 
% 
% For a set of n samples on the unit circle (angles), we define the 
% circular median as the "half-point", the sample (angle) for which half of
% the samples are larger (mod 2*pi) and half of the samples are smaller. 
%
% Due to the symmetry of the circle, this definition has two solution, and
% hence the circular median is in fact a diameter of the unit circle that
% divides the set into two equal (each sized (n-1)/2, or (n-2)/2) sets. To
% derive a single angle (radius) from this line (diameter), we can define
% the radius which is closer to the center-of-mass of samples (i.e. the 
% circular mean).

%------------------------------------------------------------------------
% Algorithm
%------------------------------------------------------------------------
% Input: sample of n points on the unit circle. 
% 
% (1) Start with any one sample. Compute the signed (e.g. clockwise) 
% circular distances between that point and each of the other points. 
% Count the number of distances that are below half-circle (0<d<pi) and the
% number above (pi<d<2*pi). Repeat for all points. 
% 
% Note that a point exactly opposite a point (or exactly overlapping a
% point) should not contribute to the count at either half-circle (hence
% the "<" in step 1, as opposed to "<=")
% 
% (2) For each sample, compute the assymetry - the difference between the
% number of points that are in each of the half-circles. 
% 
% (3) Select the "half-points" - those that have the minimal assymetery. 
%
% Note that the number of diameters can be as small as one, and as 
% large as n/2 (imagine n points equally spaced along the unit circle). 
% 
% (4) If multiple half-points are found, reflect them all to the same
% half-circle and then compute the mean of all (reflected) half-points.
% This yields a single half-point (an angle on the circle, or radius).
%
% Reflection is done by computing the mean of all half-points, computing
% the distances from the mean, and then adding a half-circle (pi) to those
% farther from the mean (then mod 2*pi). 
% 
% Note that in the extreme case of (an even) n equally-spaced points, the 
% circular mean is undefined (and indeed, in this case there are multiple
% equally-likely diameters) and thus the circular median is ill-defined
% (there are multiple solutions). 
% In such a case, the resultant vector will be close to zero. However, 
% there will always be a numerical solution and thus on a digital computer
% with finite resolution a a single solution will be selected (if the
% resolution is infinite, we can alternatively choose the angle closer to
% e.g. zero, or pi). 
% 
% (5) Now consider the two complementary radii. For the first (the
% half-point from step 4), count how many samples are up to a
% quarter-circle away (absolute circular distance of up to a
% quarter-circle, or pi/2). Do the same for the reflection of that
% half-point, and choose the half-point closer to more samples; this 
% is the median.
% 
% (6) Verify that this is indeed a valid solution by counting the number of
% points on each side of the median (see step 1). 
% 
% Output: 
% The median (an angle, possibly one of the samples in the input set)
% The number of points at each side of the median


%------------------------------------------------------------------------
% Code
%------------------------------------------------------------------------

function [ med, nc ] = circ_med( x )

% initialize output
med                 = [];
nc                  = NaN( 2, 1 );

% check input
nargs               = nargin;
if nargs < 1 || isempty( x )
    return
end
[ m, n ]            = size( x );
if m * n == 1 
    med             = x;
    return
end

% apply recursion to an n-column matrix
if n > 1
    med             = NaN( 1, n );
    nc              = NaN( 2, n );
    for i           = 1 : n
        [ med( i ), nc( :, i ) ] = circ_med( x( :, i ) );
    end
    return
end

% ensure x is an m-element column vector without NaNs
x                   = mod( x, 2 * pi );
x                   = x( ~isnan( x ) );
m                   = length( x );

% (1) signed distances between each sample and all others
n1                  = NaN( m, 1 );
n2                  = NaN( m, 1 );
for i               = 1 : m
    d               = mod( x - x( i ) + 2 * pi, 2 * pi );
    n1( i )         = sum( d < pi );
    n2( i )         = sum( d > pi & d < 2 * pi );
end

% (2) assymetry - difference between number of samples in each half-circle
a                   = abs( n1 - n2 );

% (3) half-points - those samples that have minimal assymetry
minval              = min( a );
hps                 = x( a == minval );

% (4.1) reflect all half-points to lie in the same half-circle
if length( hps ) == 2
    dh              = pi - abs( pi - abs( hps( 1 ) - hps( 2 ) ) );
    if dh > pi / 2                                              % check if more than quarter-circle apart
        hps( 1 )    = hps( 1 ) + pi;                            % reflect
    end
elseif length( hps ) > 2
    mhp             = circ_mean_local( hps );                   % compute circular mean
    dh              = pi - abs( pi - abs( hps - mhp ) );        % compute distance from mean
    idx             = dh > pi / 2;                              % check if in same half-circle as mean
    hps( idx )      = hps( idx ) + pi;                          % reflect to the other half-circle
end
hps                 = mod( hps, 2 * pi );

% (4.2) reduce to a single half-point (putative median)
hp                  = circ_mean_local( hps );

% (5) choose among the two angles on the diameter defined by the half-point
meds                = mod( [ hp; hp + pi ], 2 * pi );
d1                  = pi - abs( pi - abs( x - meds( 1 ) ) );
d2                  = pi - abs( pi - abs( x - meds( 2 ) ) );
nd1                 = sum( d1 < pi/2 );
nd2                 = sum( d2 < pi/2 );
if nd1 >= nd2
    med             = meds( 1 );
else
    med             = meds( 2 );
end

% (6) verify by counting samples on each half circle
if nargout > 1
    d               = mod( x - med + 2 * pi, 2 * pi );
    nc( 1 )         = sum( d < pi );
    nc( 2 )         = sum( d > pi & d < 2 * pi );
end

return

%------------------------------------------------------------------------
% helper function: circ_mean_local (see also circ_mean.m)
%------------------------------------------------------------------------

function phi = circ_mean_local( t )

% argument handling and trigonometric functions
nans                = isnan( t );
n                   = sum( ~nans, 1 );
x                   = cos( t );
y                   = sin( t );
% compute direction
sumx                = nansum( x, 1 );
sumy                = nansum( y, 1 );
C                   = sumx ./ n;
S                   = sumy ./ n;
phi                 = mod( atan2( S, C ), 2 * pi );

return

% EOF

