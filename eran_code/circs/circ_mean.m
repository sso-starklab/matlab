% CIRC_MEAN         Compute mean direction - matrix columns
%
% call              [ phi, R, s ] = circ_mean( t, f, sf )
%
% gets              t           matrix of angles
%                   f           corresponding matrix of weights
%                   sf          {0} optional mixing factor that determines
%                                   rescaling and reordering of output arguments; 
%                                   0: [ phi R SD ]
%                                   1: [ sem R phi ]
%                                   2: [ R phi s ]
%
% returns           phi         mean direction, computed for each column separately
%                   R           resultant (")
%                   s           circular SD
%
% calls             nothing
%
% note: NaNs and weights are supported

% 26-mar-04 ES

% revisions
% 04-apr-13 if any of the input are vectors, they are expanded
% 12-dec-19 cleaned up
% 27-jul-21 R limited to the 0-1 range

function [ phi, R, s ] = circ_mean( t, f, sf )

% argument handling and trigonometric functions
nargs                   = nargin;
if nargs < 2 || isempty( f )
    nans                = isnan( t );
    n                   = sum( ~nans, 1 );
    x                   = cos( t );
    y                   = sin( t );
else
    if size( t, 2 ) == 1
        t               = t * ones( 1, size( f, 2 ) );
    elseif size( f, 2 ) == 1
        f               = f * ones( 1, size( t, 2 ) );
    end
    if ~isequal( size( f ), size( t ) )
        error( 'input size mismatch' )
    end
    n                   = nansum( f );
    x                   = f .* cos( t );
    y                   = f .* sin( t );
end

% compute direction
sumx                    = nansum( x, 1 );
sumy                    = nansum( y, 1 );
C                       = sumx ./ n;
S                       = sumy ./ n;
phi                     = mod( atan2( S, C ), 2 * pi );

% compute dispersion
if nargout > 1 || nargs >= 3
    R                   = ( C.^2 + S.^2) .^ 0.5;
    R( R > 1 )          = 1;
    R( R < 0 )          = 0;
end
if nargout > 2 || nargs >= 3
    s                   = ( -2 * log( R ) ) .^ 0.5;
end

% optionally reorder output
if nargs >= 3
    if sf == 1
        % output is: [ sem R phi ]
        sem             = s / sqrt( n );
        phi0            = phi;
        phi             = sem;
        s               = phi0;
    elseif sf == 2
        % output is [ R phi s ]
        phi0            = phi;
        phi             = R;
        R               = phi0;
    end
end

return

% EOF
