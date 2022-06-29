% plotProbeHFOsPeriods      rudimentary version based on plotProbeHFOs

% 01-dec-20 ES & SSo

function [ fig, foo, too, yoo ] = plotProbeHFOsPeriods( filebase, periods, specChannel, channels, varargin )

% constants
tempSD_DEFAULT            	= 0;
spatBin_DEFAULT           	= 1;
nT_calc_DEFAULT             = [ -1 1 ]; % number of whole periods (T) before/after center of time ranges

% arguments
nargs                       = nargin;
if nargs < 4 || isempty( periods ) || isempty( specChannel ) || isempty( channels )
    return
end
[ plotSpec, tempSD, spatBin, nT_calc ] = ParseArgPairs(...
    { 'plotSpec', 'tempSD', 'spatBin', 'nT_calc' }...
    , { 1, tempSD_DEFAULT, spatBin_DEFAULT, nT_calc_DEFAULT }...
    , varargin{ : } );

if plotSpec == 1
    avgMode                 = 'wlt';
else
    avgMode                 = 'eeg';
end

% par
par                         = LoadXml( filebase );
try
    ngroups                 = length( par.SpkGrps );
catch
    ngroups                 = 1;
end
probe                       = get_probe( par, 0 );
if isempty( probe )
    probe                   = 1;
end
nsites                      = size( probe, 1 );
scalefactor                 = par.VoltageRange / 2^par.nBits / par.Amplification * 1e3; % A2DU to mV (not uV)
spkFs                       = par.SampleRate;

% derive trigs and staWin from periods
T                           = mean( diff( periods, [], 2 ) + 1 );           % [samples]
staWin                      = nT_calc * T / spkFs;                          % [s]
trigs                    	= mean( periods, 2 ) / spkFs;                   % [s]

% compute
[ avgcsd, avglfp, tim, ~, ~, too, foo, yoo0, yoo1 ] = pt_avg( par...
    , channels, trigs, 'specChans', specChannel...
    , 'stawinSEC', staWin, 'tempSD', tempSD, 'spatBin', spatBin...
    , 'suffix', avgMode, 'graphics', 0 ...
    , 'scalefactor', scalefactor );
% colormap( myjet )
% subplot( 2, 2, 2 )
% alines( [ -0.5 0.5 ] * T / spkFs * 1000, 'x', 'linestyle', '--', 'color', [ 1 1 1 ] );

% plot the CSD
yoo                         = nanmean( yoo1 + eps, 3 );

% plotting parameters
xCalib                      = 10;                               % 10 ms 
yCalib                      = 50;                               % 50 Hz
calib                       = [ xCalib yCalib ];
CALIBCOLOR                  = [ 0 0.7 0 ];
CALIBWIDTH                  = 2;
USF                         = 10; 

fig                         = figure;
plotSpectrogram( too * 1000, foo, yoo', USF, 'contour', 'linear' );
alines( [ 80 160 240 ], 'y', 'color', [ 1 1 1 ], 'linestyle', '--' );
xlabel( 'Time [ms]' )
ylim( [ 40 260 ] )
%line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
calibration( calib, 1i* 0.05 * [ 1 1 ], gca, 0, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH );
axis square
ah0                         = gca;
set( ah0, 'tickdir', 'out', 'box', 'off' )
alines( [ -0.5 0.5 ] * T / spkFs * 1000, 'x', 'linestyle', '--', 'color', [ 1 1 1 ] );

return

% EOF

L0 = load( [ filebase '.sps' ], '-mat' )
L0.vote % these are the channels with the peak ripple power


L = load( [ fileparts( fileparts( fileparts( filebase ) ) ) '/mat/mC400_21_stims_42_0.3V' ], '-mat' )
periods = L.stims.times; % @ spkFs
specChannel = 18;
channels = 17 : 24;
plotProbeHFOsPeriods( filebase, periods, specChannel, channels )
fig1 = gcf;
clim1 = get( gca, 'CLim' );
title( '' )
% alines( [ 80 160 240 ], 'y', 'color', [ 1 1 1 ], 'linestyle', '--' );
% ylim( [ 40 260 ] )

L = load( [ fileparts( fileparts( fileparts( filebase ) ) ) '/mat/mC400_21_stims_42_0.1V' ], '-mat' )
periods = L.stims.times; % @ spkFs
specChannel = 18;
channels = 17 : 24;
plotProbeHFOsPeriods( filebase, periods, specChannel, channels )
fig2 = gcf;
clim2 = get( gca, 'Clim' );
set( gca, 'Clim', clim1 )
title( '' )
% alines( [ 80 160 240 ], 'y', 'color', [ 1 1 1 ], 'linestyle', '--' );
% ylim( [ 40 260 ] )

