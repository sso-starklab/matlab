% cch_deconv        deconvolve ACH from CCH to esimate the STC
%
% call              [ ht, t ] = cch_deconv( cch, ach1, nspks1 )
%                   [ ..., kidx ] = cch_deconv( ..., ach2, nspks2, method, nfft )
% 
% gets              cch                 CCH (vector or matrix), [counts]
%                   ach1                ACH of the trigger spike train
%                   nspks1              number of spike in the trigger train
%                   ach2                ACH of the referred train
%                   nspks2              number of spike in the referred train
%                   method              {'fft'}         FFT deconvolution (inverse filtering)
%                                       'wiener'        Wiener deconvolution (inverse and lowpass filtering)
%                   nfft                argument 4 of deconvfft
%                   dcmode              {1|2};  1   deconvolve ach1 from cch
%                                               2   deconvolve ach1 then ach2 from the cch
%                                               3   convolve ach1 and ach2 before deconvolving from the cch
%
% returns           ht                  the STC (spike transmission curve), [counts]
%                   t                   time vector, [samples]
%                   kidx                CCH bins used for ht, [indices]
%
% does              (1) removes extra samples from ACH to preserve phase
%                   (2) scale the ACH to integral of zero
%                   (3) sets zero-lag bin of ACH to number of spikes to
%                       conserve impulse-response behavior
%                   (4) call a computational routine
%
% calls             deconvfft, firfilt
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
% Thus, if x( t ) and y( t ) are known, we can find h( t ) using
%
%               H( f ) = Y( f ) / X( f )                            (8)
%
% In the case of spike trains and second order statistics ACH and CCH, the
% formulation equivalent to eq. (1) is
%
%               cch( t ) = ach( t ) * h( t ) + e( t )               (9)
%
% Since we are typically interested in h( t ) but measure cch( t ), the
% approach used in eq. (8) can be used to recover h( t ) by an inverse
% Fourier transform of:
%
%               H( f ) = CCH( f ) / ACH( f )                        (10)
%
% The Wiender deconvolution version (in the frequency domain) would be
%
%               H( f ) = CCH / ACH * CCH^2 / ( CCH^2 + N^2 )        (11)
%
% Which attenuates the high-frequency noise. Regardless of the specific 
% filter (eq. 10 or eq. 11), the convolution of two m-element vectors is 2*m-1. 
%
% Thus, to obtain an m-element h( t ), the CCH must be at least 2*m-1 long. 
% 
% Second, the zero-lag bin of the ACH must be set to the spike count (if
% the ACH is flat, this would be equivalent to deconvolving a delta
% function, keeping h=cch. 
%-------------------------------------------------------------------------

% 14-oct-20 ES

% revisions
% 17-oct-20 methods supported: 'fft' and 'wiener'
% 24-nov-20 nfft added as input argument
% 03-dec-20 added option for sequential deconvolution
% 06-dec-20 scale ach to integral of zero before making a delta function
% 10-dec-20 corrected bug in line 215
% 11-apr-21 modified ach scaling to intergral of exactly 1 (previously,
%               was intergral of approximately 1, and center of exactly 1)

function [ ht, t, kidx ] = cch_deconv( cch, ach1, nspks1, ach2, nspks2, method, nfft, dcmode )

% argument handling
nargs                       = nargin;
if nargs < 3 || isempty( cch ) || isempty( ach1 ) || isempty( nspks1 )
    error( 'missing arguments' )
end
if isvector( cch )
    cch                     = cch( : );
end
if isvector( ach1 )
    ach1                    = ach1( : );
end
[ m,  n0 ]                  = size( cch );
[ m1, n1 ]                  = size( ach1 );
if m ~= m1
    error( 'input size mismatch: cch and ach1 must have same number of rows' )
end
if n0 > 1 && n1 > 1 && n0 ~= n1 
    error( 'input size mismatch: cch and ach1 must have same number of columns' )
end
if n0 > 1 && n1 == 1
    ach1                        = ach1 * ones( 1, n0 );
end
if numel( nspks1 ) == 1 && n0 > 1
    nspks1                      = nspks1 * ones( 1, n0 );
end

% argument handling - dual deconvolution 
if nargs < 4 || isempty( ach2 )
    ach2                        = [];
end
if nargs < 5 || isempty( nspks2 )
    nspks2                      = [];
end
if ~isempty( ach2 ) && ~isempty( nspks2 )
    if isvector( ach2 )
        ach2                    = ach2( : );
    end
    [ m2, n2 ]                  = size( ach2 );
    if m ~= m2
        error( 'input size mismatch: cch and ach2 must have same number of rows' )
    end
    if n0 > 1 && n2 > 1 && n0 ~= n2
        error( 'input size mismatch: cch and ach2 must have same number of columns' )
    end
    if n0 > 1 && n2 == 1
        ach2                    = ach2 * ones( 1, n0 );
    end
    if numel( nspks2 ) == 1 && n0 > 1
        nspks2                  = nspks2 * ones( 1, n0 );
    end
    if ~exist( 'dcmode', 'var' ) || isempty( dcmode )
        dcmode                  = 2;
    end
else
    dcmode                      = 1;
end

% deconvolution method
if nargs < 6 || isempty( method )
    method                      = 'fft';
end
if ~ismember( method, { 'fft', 'wiener' } )
    error( 'unsupported deconvolution method' )
end

% nfft
if nargs < 7 || isempty( nfft )
    nfft                        = [];
end
if isempty( dcmode )
    dcmode                      = 1;
end

% prepare
nBins                           = ( m - 1 ) / 2;
hw                              = floor( nBins / 2 );
nw                              = 2 * hw + 1;
na                              = m - nw + 1;
kidx                            = ( ( m - na ) / 2 + 1 ) : m - ( m - na ) / 2;
ach1k                           = ach1( kidx, : );
ach1k                           = ach1k - sum( ach1k ) / nw;
ach1k                           = ach1k ./ nspks1;

% deconvolve
switch dcmode
    case 1
        
        % deconvolve ach1 from cch
        hidx                    = [ 1 : hw ( hw + 2 ) : nw ];
        ach1k( hw + 1, : )      = 1 - sum( ach1k( hidx, : ) );
        ht                      = deconvfft( cch, ach1k, hw, nfft, method );  % will also preserve phase if ach1k is same dimension as cch
        % organize time vector
        t                      	= ( -hw : hw )';

    case 2
        
        % deconvolve ach1 from cch, then ach2 from the remainder
        hidx                    = [ 1 : hw ( hw + 2 ) : nw ];
        ach1k( hw + 1, : )      = 1 - sum( ach1k( hidx, : ) );
        
        ht1                     = deconvfft( cch, ach1k, hw, nfft, method );  % will also preserve phase if ach1k is same dimension as cch
        
        kidx1                   = kidx;
        ach2k                   = ach2( kidx1, : );
        m2                      = length( kidx1 );
        nBins2               	= ( m2 - 1 ) / 2;
        hw2                  	= floor( nBins2 / 2 );
        nw2                   	= 2 * hw2 + 1;
        na2                   	= m2 - nw2 + 1;
        kidx2                	= ( ( m2 - na2 ) / 2 + 1 ) : m2 - ( m2 - na2 ) / 2;
        ach2k                   = ach2k( kidx2, : );
        ach2k                   = ach2k - sum( ach2k ) / nw2;
        ach2k                   = ach2k ./ nspks2;
        
        hidx                    = [ 1 : hw2 ( hw2 + 2 ) : nw2 ];
        ach2k( hw2 + 1, : )     = 1 - sum( ach2k( hidx, : ) );
        
        ht2                     = deconvfft( flipud( ht1 ), ach2k, hw2, nfft, method );
        ht                      = flipud( ht2 );
        kidx                    = kidx1( kidx2 );
        
        % organize time vector
        t                    	= ( -hw2 : hw2 )';
        
    case 3
        
        % convolve the ACHs before deconvolution
        ach2k                   = ach2( kidx, : );
        ach2k                   = ach2k - sum( ach2k ) / nw; 
        ach2k                   = ach2k ./ nspks2;
        ach12k                  = firfilt( ach1k, ach2k );
        hidx                    = [ 1 : hw ( hw + 2 ) : nw ];
        ach12k( hw + 1, : )     = 1 - sum( ach12k( hidx, : ) );
        
        % deconvolve the ACH product
        ht                      = deconvfft( cch, ach12k, hw, nfft, method );
        % organize time vector
        t                      	= ( -hw : hw )';

end

% organize
ht( ht < 0 )                    = 0;                                        % clip negatives (may occur numerically) to zero

return

% EOF

