% classify_waveform         from peak spectrum and trough-peak time
%
% CALL                      [ wide, pval, fsep ] = classify_waveform( xy )
%
% GETS                      xy      two column vector:
%                               first column is the trough-peak time (ms)
%                               second column is the spike width (1/freq*1000)
%                               classification is by GMM and a p-value is assigned
%
% RETURNS                   wide    a logical vector:
%                                   1 for a "pyramidal" cell (wide waveform, long trough-to-peak time), 
%                                   0 for an "interneuron" (narrow waveform, short peak-to-trough)
%
% CALLS                     nothing (hard-coded separatrix)
%
% NOTE: 
% (1) the separatrix was derived from mono-synaptic connectivity
%   patterns in CX and CA1 of mice and rats and from light activation of PV
%   cells in PV mice, but may not hold for other brain regions/species
% (2) xy should be derived from detrended waveforms (not high-pass filtered)
% (3) the tp and width parameters discriminate much better than half-width, asymmetry,
%   rate, first moment of the ACH or any combination thereof
%   see spike_stats.m for methods to compute these parameters from raw data

% 25-jan-12 ES

% revisions
% 25-may-13 GMM added
% 23-feb-18 GMM by default
% 30-aug-18 hpf temporary workaround
% 17-aug-19 cleaned up

function [ wide, pval, fsep ] = classify_waveform( xy, ft, hpf )

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( xy )
    wide                    = [];
    pval                    = [];
    return
end
if nargs < 2 || isempty( ft )
    ft                      = 0;    % 1: use a linear classified (obsolete)
end
ft                          = ft( 1 );
if size( xy, 2 ) ~= 2
    error( 'input size mismatch' )
end
if nargs < 3 || isempty( hpf )
    hpf                     = 0;    % 1: use median-filtered waveforms (suboptimal)
end
if hpf
    offset                  = -0.17;
else
    offset                  = 0;
end

% classify
if ~ft                              % use the GMM
    mixp                =  [   0.14815     0.85185 ];
    mu                  =  [   0.328       0.90036
                               0.6491      1.078 ] + offset;
    Sigma(:,:,1)        = [
                0.0085982   0.0052262
                0.0052262   0.0056747 ];
    Sigma(:,:,2)        = [
                0.0020751   0.0013943
                0.0013943   0.0049193 ];
    gm                  = gmdistribution( mu, Sigma, mixp );   
    p                   = posterior( gm, xy );
    [ ~, cls ]          = max( p, [], 2 );
    pval                = min( p, [], 2 );
    wide                = cls == 2;
else                                % use a linear separatrix - OBSOLETE
    % in this case, xy columns are: [ frequency, trough-to-peak ]
    xx                  = [ 600 1450 ];
    yy                  = [ 0.2 0.9 ];
    a                   = diff( yy ) / diff( xx ); 
    b                   = yy( 1 ) - a * xx( 1 );        % y = ax+b
    xint                = -b/a;                         % x-intercept
    wide                = atan( xy( :, 2 ) ./ ( xy( :, 1 ) - xint ) ) > atan( a );
    pval                = NaN * ones( size( wide ) );
end

if nargout > 2
    if ft
        fsep            = @(x,y) b + a*x - y;
    else
        K               = 48.905;
        L               = [   -335.3  117.01 ]';
        Q               = [   166.3   36.272
                               36.272  -72.892 ];
        fsep            = @(x,y) K + L(1)*x + L(2)*y ...
            + Q(1,1)*x.^2 + (Q(1,2)+Q(2,1))*x.*y + Q(2,2)*y.^2;
    end
end

return

% EOF

% to actually plot the separatrix: 
xx = [ 0.1 0.9 ]; yy = [ 0.6 1.4 ];
fh = fimplicit( fsep, [ xx yy ] ); 
set( fh, 'color', [ 0 0 0 ] );
title( '' )