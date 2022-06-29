% CIRCLIN       regress a circular variable on a linear variable
%
% call          [ SLOPE, INTERCEPT, R, PVAL, Rs ] = CIRCLIN( THETA, X, EXTREMA, TOL, N, GRAPHICS, BLANK )
%
% does          fits the model: THETA = INTERCEPT + 2 * pi * SLOPE * X
%
% algorithm     given a circular dependent variable THETA, 
%               and a linear independent variable X:
%               (1) performs a global search (between EXTREMA(1) and EXTREMA(2)), 
%               at resolution 1000*TOL), ignoring data in BLANK (returned as Rs)
%               (2) refines using local optimization up to tolerance TOL
%               (3) computes significance by permuting the theta-x pairs N times
%               (4) if desired, plots
%
% defaults      EXTREMA     {[-1 1 ]}; can be overloaded by As
%               TOL         {1e-5}
%               N           {0}
%               GRAPHICS    {0}
%               BLANK       {[-0.03 0.03]}
%
% by default, the global search grid is on equally-spaced slopes:
%       As          := exterma( 1 ) : stepsize : extrema( 2 ) 
%       stepsize    := tol * 1000
%
% if extrema is overloaded (any vector of 3 or more elements), the grid
% search will be on those elements. an example is equally-scaled angles: 
%       As          = tan( 0 : pi/2/100 : pi/2 - pi/2/100 );
%
% calls         mixmat
%
% reference     Schmidt et al. 2009 JNS

% 27-may-11 ES

% revisions
% 29-dec-19 blanking of zero-slope during global search added
% 29-dec-19 (1) values of global search (Rs) returned
%           (2) permuation test adapted to blanking
% 21-jan-20 (1) extrema overloaded with As

function [ a, phi0, R, pval, Rs ] = circlin( theta, x, extrema, tol, nrands, graphics, blank )

mult                        = 1e3;
blanker                     = [ -3 3 ];

% arguments
nargs = nargin;
if nargs < 2 || isempty( theta ) || isempty( x )
    error( 'missing inputs' )
end
theta                       = theta( : );
x                           = x( : );
if length( x ) ~= length( theta )
    error( 'inpute size mismatch: theta/x' )
end
if nargs < 3 || isempty( extrema )
    extrema                 = [ -1 1 ];
end
As                          = [];
if length( extrema ) < 2
    error( 'input size mismatch: extrema' )
else
    As                      = sort( extrema( : ) )';
    extrema                 = [ min( As ) max( As ) ];
end
extrema                     = sort( extrema );
if nargs < 4 || isempty( tol )
    tol                     = 1e-5; % for fine tuning local search
end
tol                         = abs( tol( 1 ) );
if tol > diff( extrema )
    error( 'input size mismatch: tol/extrema' )
end
if nargs < 5 || isempty( nrands )
    nrands                  = 0; 
end 
if nargs < 6 || isempty( graphics )
    graphics                = 0; 
end
if nargs < 7 || isempty( blank )
    blank                   = blanker * mult * tol;
end

% step 1: screen the entire range for a global maximum
stepsize                    = tol * mult; % for initial search
if isempty( As )
    As                      = extrema( 1 ) : stepsize : extrema( 2 );
end
if ~isempty( blank )
    ridx                    = As >= blank( 1 ) & As <= blank( 2 );
    As( ridx )              = [];
end
n                           = length( As );
Rs                          = zeros( n, 1 );
for i                       = 1 : n
    t                       = theta - 2 * pi * As( i ) * x;
    Rs( i )                 = ( nanmean( cos( t ) ) .^ 2 + nanmean( sin( t ) ) .^ 2  ) .^ 0.5;
end
[ R0, idx ]                 = max( Rs );
a0                          = As( idx );

% step 2: iterate around the maximum
R                           = 0;
a                           = a0;
dir                         = 1;
step                        = stepsize;
maxiter                     = round( sqrt( 1/tol ) );
iter                        = 0;
while abs( R - R0 ) > tol && iter < maxiter
    R0                      = R;
    t                       = theta - 2 * pi * a * x;
    R                       = ( nanmean( cos( t ) ) .^ 2 + nanmean( sin( t ) ) .^ 2  ) .^ 0.5;
    if R > R0 % converging
        if dir == 1 % increased at previous step
            a               = a + step;
            dir             = 1;
        else
            a               = a - step;
            dir             = -1;
        end
        step                = step / 2;
    else
        if dir == 1
            a               = a - step;
            dir             = -1;
        else
            a               = a + step;
            dir             = 1;
        end
    end
    iter                    = iter + 1;
end

% step 3: return the values
t                           = theta - 2 * pi * a * x;
phi0                        = mod( atan2( sum( sin( t ) ), sum( cos( t ) ) ), 2 * pi );

% compute a p-value (recursion)
if nrands > 0
    mix                     = mixmat( ( 1 : size( theta, 1 ) )' * ones( 1, nrands ), 1, 1 );
    ros                     = rand( nrands ) * 2 * pi; % randomize the offset
    rmix                    = zeros( nrands, 1 );
    amix                    = zeros( nrands, 1 );
    for i                   = 1 : nrands
        tmix                = theta( mix( :, i ) ) + ros( i );
        [ amix( i ), ~, rmix( i ) ]    = circlin( tmix, x, extrema, tol, 0, 0, blank );
    end
    if ~isempty( blank )
        ridx                = amix >= blank( 1 ) & amix <= blank( 2 );
        amix( ridx )        = [];
        rmix( ridx )        = [];
        nrands              = sum( ~ridx );
    end    
    if nrands <= 100
        pval                = 1 - normcdf( atanh( R ), mean( atanh( rmix ) ), std( atanh( rmix ) ) );
    else
        pval                = ( sum( R < rmix ) + 1 ) / ( nrands + 1 );
    end
else
    pval                    = 1;
end

% graphical report
if graphics
    
    newplot
    clf
    
    subplot( 2, 1, 1 )
    plot( [ x; x ], [ theta; theta + 2 * pi ], 'ob' )
    xlabel( 'x' )
    ylabel( '{\theta}' )
    thetahat                = 2 * pi * a * x + phi0;
    title( sprintf( 'slope=%0.3g; intercept=%0.3g', a, phi0 ) )
    axis tight
    
    subplot( 2, 2, 4 )
    plot( As, Rs, '.-' )
    title( 'screening' )
    xlabel( 'slope' )
    ylabel( 'R' )
    
    subplot( 2, 2, 3 )
    plot( mod( theta, 2 * pi ), mod( thetahat, 2 * pi ), '.' )
    xlabel( '{\theta}' )
    ylabel( '{\theta}_{hat}' )
    title( sprintf( 'R=%0.3g; p=%0.3g', R, pval ) )
    xlim( [ 0 2 * pi ] )
    ylim( [ 0 2 * pi ] )
    
end

return

% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% example - generate some data
a = 1 * ( 2 * rand( 1 ) - 1 );
ph0 = rand( 1 ) * 2 * pi;
n = 10000;
x = sort( rand( n, 1 ) * 50 ) + 50;
phi = mod( 2 * pi * a * x + ph0, 2 * pi );
noise = rand( n, 1 ) * sqrt( 2 ) / 2 * range( phi );
% fit it with the same model
[ slope intercept R ] = circlin( phi + noise, x, [ -1 1 ], [], 10, 1 );
[ a, ph0 ]
[ slope intercept R ]


