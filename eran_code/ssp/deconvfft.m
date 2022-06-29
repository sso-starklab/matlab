% deconvfft             recover impulse response of an LTI system using frequency-domain deconvolution
%
% signal model:         y = x * h + n, 
%                           where * indicates time-domain convolution 
%                           we observe y and x, 
%                           want to determine h in the presence of noise n
%
% call                  [ h, hx, x ] = deconvfft( y, x, m, nfft, method )
%
% gets                  y           the output to the system
%                       x           the input to the system
%                       m           number of a-causal samples in h
%                                   {0}, defaults to zero, corresponding to a causal filter
%                       nfft        number of frequency domain bins
%                       method      {'fft'} or 'wiener'
%
% does                  h = ifft ( Y( f ) * R )
%                       
%                       where, for FFT filtering: 
% 
%                               R( f ) = 1 ./ X( f )
% 
%                       and for Wiener filtering:
%
%                               R( f ) = 1 / X * Y^2 / ( Y^2 + N^2 )
%                       
%                       for matrix inputs, works on columns (computes a 
%                           filter for each pair of columns) so the two
%                           matrices must have the same number of columns
%
% returns               h           impulse response of the LTI system
%                                       for a matrix input, matrix output
%                                       number of rows will be 2*m+1 (if m is not 0) or nfft (if m is 0)
%                       hx          lag vector for h
%                       x           clipped input. if m is 0, same as input
%
% note                  if m is not zero, then 
%
%                               y = conv( x, h )
%
%                       if m is zero, then 
%
%                               fft( y ) = fft( h ) .* fft( x )
%
% calls                 nothing
%
%-------------------------------------------------------------------------
% theory
%-------------------------------------------------------------------------
% The standard signal, LTI system, and Gaussian noise model is:
% 
%               y( t ) = x(t) * h(t) + n(t)                         (1)
% 
% where x(t) is some random signal, h(t) is the impulse response of the LTI
% system, and n(t) is independent Gaussian noise. Given this model, the 
% input x(t) can be reconstructed from h(t) and y(t) using deconvolution:
%
%               x( t ) = y( t ) /* h( t )                           (2)
%
% which is most easily implemented in the frequency domain using
%
%               X( f ) = Y( f ) / H( f )                            (3) 
%
% where X( f ) is the power spectrum of x( t ). If we define
%
%               R( f ) := 1 / H( f )                                (4) 
%
% then eq. (3) can be rewritten as 
%
%               X( f ) = Y( f ) * R( f )                            (5) 
%
% In the case of Wiener deconvolution, R( f ) can be defined differently
% to reduce the amplification of high-frequency noise:
% 
%               R := 1 / H * H^2 * X^2 / ( H^2 * X^2 + N^2 )        (6) 
%
% Plugging eq. (3) into eq. (6) yields the Wiener filter:
%
%               R := 1 / H * Y^2 / ( Y^2 + N^2 )                    (7)
%
% Regardless of whether an FFT filter (eq. 4) or a Wiener filter
% (eq. 7) is used in the deconvolution process (eq. 5), the construction of
% the system (eq. 1) is symmetric. Therefore, the roles of h( t ) and 
% x( t ) can be exchanged. 
%
% Thus, if x( t ) and y( t ) are known, we can use 'fft' filtering to find
% h( t ) using the inverse Fourier transform of
%
%               H( f ) = Y( f ) / X( f )                            (8)
%
% Or 'wiener' filtering to find h( t )
%
%               G( f ) = Y( f )^2 / ( Y( f )^2 + N( f )^2 )         (9.1)
%               H( f ) = Y( f ) / X( f ) * G( f )                   (9.2)
%
% The noise spectrum N( f ) is assumed to be Gaussian (white), and can be
% estimated from the tails of y.
%-------------------------------------------------------------------------

% 27-jun-19 ES

% revisions
% 14-oct-20 (1) matrix input supported
%           (2) real part of the signal taken
% 17-oct-20 (1) wiener filtering added
%           (2) documentation expanded

function [ h, hx, x ] = deconvfft( y, x, m, nfft, method )

% check input size, make into column vectors or matrices with same number of columns
nargs                           = nargin;
if nargs < 2 || isempty( y ) || isempty( x )
    return
end
if ndims( x ) == 1
    y                           = y( : );
end
if ndims( y ) == 1
    x                           = x( : );
end
nx                              = size( x, 1 );
ny                              = size( y, 1 );
if size( x, 2 ) ~= size( y, 2 )
    return
end
if nargs < 3 || isempty( m )
    m                           = 0; 
end
if nargs < 4 || isempty( nfft )
    nfft                        = max( nx, ny ); 
end
if nargs < 5 || isempty( method )
    method                      = 'fft';
end
method                          = lower( method );
if ~ismember( method, { 'fft', 'wiener' } )
    error( 'unsupported deconvolution method' )
end

% clip input to preserve phase
if m ~= 0
    nh                          = 2 * m + 1;                                % size of a-causal h with symmetric support
    d                           = nx - ny + nh - 1;                         % number of samples to remove from input vector
    d                           = max( d, 0 );                              % cannot remove negative samples (may occur if x already shifted)
    kidx                        = ( d / 2 + 1 ) : ( nx - d / 2 );
    x                           = x( kidx, : );
    nx                          = length( kidx );
end

% compute FFT
Sx                              = fft( x, nfft );
Sy                              = fft( y, nfft );
if strcmp( method, 'fft' ) || nx >= ny
    F                           = 1 ./ Sx;                                  % phase-preserving inverse filter
else
    k                           = floor( ( ny - nx ) / 2 );
    tidx                        = [ 1 : k ( ny - k + 1 : ny ) ];
    vn                          = ones( nfft, 1 ) * var( y( tidx, : ), [], 1 );
    G                           = abs( Sy ).^2 ./ ( abs( Sy ).^2 + vn );    % denoising filter
    F                           = G ./ Sx;                                  % Wiener filter
end
h                               = real( ifft( Sy .* F, nfft ) );

% shift back and determine lags
if m == 0
    hx                          = 1 : nfft;
else
    h                           = h( 1 : nh, : );
    hx                          = -m : m;
end

return

% EOF
