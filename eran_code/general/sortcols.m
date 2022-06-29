% sortcols              of a matrix by the max/com value
%
% call                  [ xs, sidx ] = sortcols( x, scalemax, sortmode )
%
% gets                  x               matrix
%                       scalemax        {1};    1/2:    scale each column by its max value before sorting 
%                                               2:      also return max-scaled matrix
%                                               0:      do not scale columns
%                       sortmode        {1};    1:      by max
%                                               2:      by center-of-mass
%
% returns               xs              sorted matrix
%                       sidx            column indices such that xs = x( :, sidx );
% 
% calls                 calc_com

% 25-dec-12 ES

% revisions
% 17-apr-13 fixed erratic behavior of NaN columns
% 23-mar-21 cleaned up

function [ xs, sidx ] = sortcols( x, scalemax, sortmode )

nargs                           = nargin;
if nargs < 2 || isempty( scalemax )
    scalemax                    = 1;                                        % 1: scale by max; 0: do not
end 
if nargs < 3 || isempty( sortmode )
    sortmode                    = 1;                                        % 1: max; 2: center-of-mass
end 

m                               = size( x, 1 );
nans                            = sum( isnan( x ) ) == m;
x( :, nans )                    = [];

if scalemax
    mm                          = max( x );
    xn                          = bsxfun( @rdivide, x, mm );
else
    xn                          = x;
end

switch sortmode
    case -1 % min
        [ ~, maxidx ]           = min( xn );
    case 1 % max
        [ ~, maxidx ]           = max( xn );
    case 2 % center of mass
        maxidx                  = round( calc_com( ( 1 : size( x, 1 ) )' * ones( 1, size( x, 2 ) ), x ) );
    otherwise
        maxidx                  = 1 : size( xn, 1 );
end

[ ~, sidx ]                     = sort( maxidx );
switch scalemax
    case 2
        xs                      = xn( :, sidx );
    otherwise
        xs                      = x( :, sidx );
end

if sum( nans )
    xs                          = [ xs NaN * ones( m, sum( nans ) ) ];
    sidx0                       = sidx;
    fidx                        = find( ~nans ); 
    sidx                        = [ fidx( sidx0 ) find( nans ) ];
end

return

% EOF
