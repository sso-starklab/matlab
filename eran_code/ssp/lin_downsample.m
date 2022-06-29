% lin_downsample        downsample with a lowpass
% 
% CALL          y = lin_downsample( x, dsf, islin )
% 
% GETS          x             matrix (columns downsampled)
%               dsf           downsampling factor (integer); negative number indices number of
%                               samples in the output; vector indicates specific samples
%               islin         vector (same number of columns as in x)
%
% DOES:         low-pass filters, then downsamples. 
%               the dsf/2's sample of the input is the 1st output
%               for instance, if dsf is 32, then the first 3 samples of the
%               output are samples 16, 48, and 80 of the filtered input
%
% CALLS         makefir, firfilt, circfilt

% 03-feb-13 ES

% revisions
% 25-feb-18 (1) converted to a targeted lowpass filter
%           (2) support of circular variables
%           (3) improved argument handling

function y = lin_downsample( x, dsf, islin )

% initialize output
y               = [];

% argument check
nargs = nargin;
if nargs < 2 || isempty( x ) || isempty( dsf )
    return
end
if nargs < 3 || isempty( islin )
    islin               =  true( 1, size( x, 2 ) );
end
if ~isnumeric( x )
    error( 'argument size/type mismatch' )
end
if isvector( x )
    x                   = x( : );
end
n                       = size( x, 1 );
islin                   = islin( : ).';
if size( x, 2 ) ~= length( islin ) || any( ~islogical( islin ) )
    error( 'argument size/type mismatch' )
end

% determine downsampling indices
if length( dsf ) > 1
    didx                = dsf;
else
    if dsf > 1 % downsample
        didx            = round( dsf / 2 : dsf : n );
        MA              = dsf;
    elseif dsf < -1
        ns              = abs( dsf );
        didx            = round( [ fliplr( n / 2 : -n/ns : 1 ) n / 2 + n / ns : n/ns : n ] );
        MA              = ceil( mean( diff( didx ) ) );
    end
end
if ~exist( 'didx', 'var' )
    return
end
didx( 1 )               = max( didx( 1 ), 1 );
didx( end )             = min( didx( end ), n );

% filter
win                     = makefir( [ 0 ceil( 1000 / MA ) ], 1000, [], 'low' );
win                     = win( : ) / sum( win );
if all( islin == 1 )
    xf                  = firfilt( x, win );
elseif all( islin == 0 )
    xf                  = circfilt( x, win );
else
    xf                  = zeros( size( x ) );
    xf( :, islin )      = firfilt( x( :, islin ), win );
    xf( :, ~islin )     = circfilt( x( :, ~islin ), win );
end

% downsample    
y                       = xf( didx, : );

return

% EOF
