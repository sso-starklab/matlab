% CALC_SPK_SNRS         compute single-spike SNR
%
% CALL                  [ snr, chan, snrs ] = calc_spk_snrs( spks, mat, ndim )
%
% GETS                  SPKS        nchans x nsamples x nspikes (or nsamples x nspikes)
%                       MAT         optional: basis to project spikes onto (library PCs)
%                       NDIM        force 3-dim structure (relevant for single spikes)
%
% DOES: 
% for each spike on each channel:
%
% if a basis (MAT) is given, projects spikes on it and computes the ratio between the
% signal (energy of the reconstruction at the best offset) and the noise (residuals)
% in this case the asymptotic value of SNR in case of noise depends on the basis
%
% otherwise, computes the ratio between the signal (max-min of central third of the waveform) 
% and the noise (the two-thirds of the waveform, i.e. the tails), with a normalization factor that
% ensures snr => 1 in the case of gaussian noise
% 
% then snr = max( snrs )
%
% CALLS                 optproj

% 22-apr-10 ES

% revisions
% 01-nov-10 2D input supported as well
% 24-mar-11 protect from case of a dead offset channel
% 02-apr-10 added option for projection on a basis (e.g. library PCs)
% 07-jun-12 added argument ndim to force structure of nchans x nsamples x
%               nspikes (relevant for the 1-spike case only)
% 17-aug-19 cleaned up

function [ snr, chan, snrs ] = calc_spk_snrs( spks, mat, ndim )

% arguments
nargs               = nargin;
if nargs < 2 || isempty( mat )
    f               = 1;
else
    f                   = 0;
end
if nargs < 3 || isempty( ndim )
    ndim            = ndims( spks );
end
if ndims( spks ) > ndim
    ndim            = ndims( spks );
end
switch ndim
    case 2
        [ nsamples, nspikes ]               = size( spks );
    case 3
        [ nchans, nsamples, nspikes ]       = size( spks );
end

% compute
si                                          = ceil( nsamples / 3 ) + 1 : nsamples - ceil( nsamples / 3 );
ni                                          = [ 1 : ceil( nsamples / 3 ) nsamples - ceil( nsamples / 3 ) + 1 : nsamples ];
switch ndim
    case 2
        if f
            snrs                            = ( ( max( spks( si, : ), [], 1 ) - min( spks( si, : ), [], 1 ) ) ./ std( spks( ni, : ), [], 1 ) ) / log( length( ni ) );
        else
            [ ~, snrs ]                     = optproj( spks, mat );
        end
        snr                                 = snrs;
        chan                                = ones( size( snr ) );
    case 3
        if f
            snrs                            = squeeze( ( max( spks( :, si, : ), [], 2 ) - min( spks( :, si, : ), [], 2 ) ) ./ std( spks( :, ni, : ), [], 2 ) ) / log( length( ni ) );
        else
            snrs                            = zeros( nchans, nspikes );
            for i                           = 1 : nchans
                [ ~, snrs( i, : ) ]         = optproj( squeeze( spks( i, :, : ) ), mat );
            end
        end
        snrs( isinf( snrs ) )               = NaN;
        [ snr, chan ]                       = max( snrs );
end

return

% EOF

