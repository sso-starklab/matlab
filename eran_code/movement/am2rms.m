% AM2RMS            convert AM data to RMS values
%
% [ Y XF ] = AM2RMS( X, FS, GRAPHICS, ZFLAG )
%
% X             multi-channel accelerometer data (columns), e.g. 3-axis from ADXL330
% Fs            {1250}; sampling rate, Hz
% GRAPHICS      {0}; optional (2 figs)
% ZFLAG         {0}; standardize before computing RMS
%                       (does not enable comparison of diff. state files)
%
% XF            low-pass filtered (10 Hz) readings, postural effects (2 sec MA) removed 
% Y             RMS (50 ms) of the smoothed and baseline-corrected data, averaged across all channels
%
% thus, if X is in units of g, XF will be measured in changes from the
% temporally local baseline (that is, acceleration relative to posture) and
% y will measure the moving-average RMS 
%
% calls          makefir, firfilt, zs, ma_rms

% 02-jun-12 ES

% revisions
% 19-nov-12 (1) no standardization by default
%           (2) vector sum of axes before RMS
%           (3) optional parameters

function [ y, xf ] = am2rms( x, Fs, graphics, zflag, varargin )

% arguments
nargs = nargin;
if nargs < 1 || isempty( x )
    y = x;
    xf = x;
    return
end
if nargs < 2 || isempty( Fs ), Fs = 1250; end
if nargs < 3 || isempty( graphics )
    graphics = 0;
end
if nargs < 4 || isempty( zflag )
    zflag = 0; 
end

% "constants"
[ MA, LP, T ] = DefaultArgs( varargin, { 2, 10, 0.05 } ); % sec, Hz, sec

% preparations
mwin = ones( ceil( MA * Fs ), 1 ); 
mwin = mwin / sum( mwin );
h = makefir( [ 0 LP ], Fs, [], 'low' );
h = h / sum( h );

% compute
xm = firfilt( x, mwin );                % time-varying mean (posture)
x0 = x - xm;                            % acceleration free of postural changes
xf = firfilt( x0, h );                  % low-pass to remove 5 Hz noise and fast transients:
if zflag
    xf = zs( xf );                      % standardize axes to same weight (bad idea)
end
r = sqrt( sum( xf .^ 2, 2 ) );          % vector sum
y = ma_rms( r, ceil( T * Fs ) );        % short-time RMS

% plot
if graphics
    tidx = 1 : min( length( y ), 100 * Fs );
    %tidx = Fs * 450 : Fs * 700;
    t = tidx / Fs;
    figure,
    subplot( 2, 2, 1 ), plot( t, x( tidx, : ) ), xlim( t( [ 1 end ] ) ), title( 'raw' )
    subplot( 2, 2, 2 ), plot( t, xf( tidx, : ) ), xlim( t( [ 1 end ] ) ), title( sprintf( 'smoothed (%0.2g sec; %0.3g Hz)', MA, LP ) )
    subplot( 2, 2, 3 ), plot( t, r( tidx, : ) ), xlim( t( [ 1 end ] ) ), title( 'vector sum' ), ylabel( 'g' ), xlabel( 'Time (sec)' )
    subplot( 2, 2, 4 ), plot( t, y( tidx ) ), xlim( t( [ 1 end ] ) ), title( sprintf( 'RMS (%0.3gsec)', T ) ), ylabel( 'RMS' )
end

return

% EOF




