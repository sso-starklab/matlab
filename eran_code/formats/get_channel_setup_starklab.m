% get_channel_setup    get channel parameters for recordings sessions
%
% intended as a source for updating an xml file
% all relevant info necessary to be packaged into a *prm.xml file is provided
% this file is partially redundant with the *xml file and readable by ndmanager
% but is not overwritten by neuroscope
%
% conventions:
%
% each channel has the following fields
% type              'neuronal', 'stim', 'sensor', 'ASD', 'trig', 'solenoid', 'sync', 'am'
% voltagerange      a number ('2', '8', '20')
%
% stim channels also have the following fields:
% target            number: shank number; negative number: above a given
%                   recording site; 0/NaN - not above any electrode (e.g.
%                   MS fiber)
% source            'LED', 'LD', 'DPSS', or 'SHUT'
% wavelength        [nm]
% power             [uW] actual max power, measured on probe (not intensity)
% voltage           [V]  measured voltage (equivalent to current for diodes) at max power
%
% other (behavioral) channels have additional fields:
% location          'left', 'right', 'center' (for 'sensor' or 'solenoid')
% reward            'water', 'sugar'
%
% see also          create_prm_file

% 07-jan-13 ES

% revisions
% 03-sep-19 modified for stark lab

function [ data, dat, animals ] = get_channel_setup( setup, filenumber, varargin )

% constants
daq                 = 'BAD';

% initialize output
data                = [];
dat                 = [];
animals             = [];

% arguments
nargs               = nargin;
if nargs < 2 || isempty( filenumber )
    filenumber      = [];
end

% in case this is a filebase:
[ ~, filename ]     = fileparts( setup );
if ~isempty( filename )
    setup           = filename;
end

% assign the values
switch setup
    
    case 'mV99' % A1x32 Edge 20 um probe
        % channel allocation
        
        % mapping
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        data{  2 } = {   	  33, 'type', 'stim' };         %
        data{  3 } = {   34 : 36, 'type', 'am' };
        data{  4 } = {        37, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 32, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = {      33, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 34 : 36, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      37, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      33, 'target', '1' };           % channel 32
        data{ 31 } = {      33, 'source', 'LED' };
        data{ 41 } = {      33, 'wavelength', '470' };  % [nm]
        data{ 51 } = {      33, 'power', '45' };        % [uW]
        data{ 61 } = {      33, 'voltage', '5' };       % [V]
        
    case 'mF84' % A1x32 Edge 20 um probe
        % channel allocation
        
        % mapping
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        data{  2 } = {        33, 'type', 'x' };            %
        data{  3 } = {        34, 'type', 'theta' };            %
        data{  4 } = {   	  35, 'type', 'stim' };         %
        data{  5 } = {   36 : 38, 'type', 'am' };
        data{  6 } = {        39, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 32, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = {      35, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 36 : 38, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      39, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      35, 'target', '1' };           % channel 32
        
        data{ 31 } = {      35, 'source', 'LED' };
        data{ 41 } = {      35, 'wavelength', '470' };  % [nm]
        
        data{ 51 } = {      35, 'power', '30' };        % [uW]
        
        data{ 61 } = {      35, 'voltage', '5' };       % [V]
        
    case 'mF91' % A1x32 Edge 100 um probe
        % channel allocation
        
        % mapping
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        data{  2 } = {   33 : 34, 'type', 'stim' };         % 2 CS channels (S#3; S#5)
        data{  3 } = {        35, 'type', 'x' };            %
        data{  4 } = {        36, 'type', 'theta' };            %
        data{  5 } = {   37 : 39, 'type', 'am' };
        data{  6 } = {        40, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 32, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = { 33 : 34, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 37 : 39, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      40, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      33, 'target', '1' };           % channel 32
        data{ 22 } = {      34, 'target', '1' };           % channel 25
        
        data{ 31 } = { 33 : 34, 'source', 'LED' };
        data{ 41 } = { 33 : 34, 'wavelength', '470' };  % [nm]
        
        data{ 51 } = {      33, 'power', '32' };        % [uW]
        data{ 52 } = {      34, 'power', '3' };        % [uW]
        
        data{ 61 } = { 33 : 34, 'voltage', '5' };       % [V]
        
    case 'mF93' % 6-shank Stark64 probe with 5 LEDs
        % channel allocation
        % shanks 1, 6:      10 sites
        % shanks 2,3,4,5:   11 sites
        
        % mapping
        data{  1 } = {    1 : 64, 'type', 'neuronal' };
        data{  2 } = {        65, 'type', 'x' };            %
        data{  3 } = {        66, 'type', 'theta' };            %
        data{  4 } = {   67 : 72, 'type', 'stim' };         % 2 CS channels (S#3; S#5)
        data{  5 } = {   73 : 75, 'type', 'am' };
        data{  6 } = {        76, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 64, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = { 67 : 72, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 73 : 75, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      76, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      67, 'target', '1' };
        data{ 22 } = {      68, 'target', '2' };
        data{ 23 } = {      69, 'target', '3' };
        data{ 24 } = {      70, 'target', '4' };
        data{ 25 } = {      71, 'target', '5' };
        data{ 26 } = {      72, 'target', '6' };
        
        data{ 31 } = { 67 : 72, 'source', 'LED' };
        data{ 41 } = { 67 : 72, 'wavelength', '470' };  % [nm]
        
        data{ 51 } = {      67, 'power', '34' };        % [uW]
        data{ 52 } = {      68, 'power', '46' };        % [uW]
        data{ 53 } = {      69, 'power', '44' };        % [uW]
        data{ 54 } = {      70, 'power', '15' };        % [uW]
        data{ 55 } = {      71, 'power', '18' };        % [uW]
        data{ 56 } = {      72, 'power', '0' };        % [uW]
        
        data{ 61 } = { 67 : 72, 'voltage', '5' };       % [V]
        
    case 'mF79' % 6-shank Stark64 probe with 2 LEDs
        % channel allocation
        % shanks 1, 6:      10 sites
        % shanks 2,3,4,5:   11 sites
        
        % mapping
        data{  1 } = {    1 : 64, 'type', 'neuronal' };
        data{  2 } = {        65, 'type', 'x' };            %
        data{  3 } = {        66, 'type', 'y' };            %
        data{  4 } = {   67 : 68, 'type', 'stim' };         % 2 CS channels (S#3; S#5)
        data{  5 } = {   69 : 71, 'type', 'am' };
        data{  6 } = {        72, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 64, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = { 65 : 68, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 69 : 71, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      72, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      67, 'target', '3' };
        data{ 22 } = {      68, 'target', '5' };
        
        data{ 31 } = { 67 : 68, 'source', 'LED' };
        data{ 41 } = { 67 : 68, 'wavelength', '470' };  % [nm]
        
        data{ 51 } = {      67, 'power', '35' };        % [uW]
        data{ 52 } = {      68, 'power', '34' };        % [uW]
        
        data{ 61 } = { 67 : 68, 'voltage', '5' };       % [V]
        
        
    case 'mP23' % 4-shank/4-LD, 2-LED probe
        layout = 540;
        if ~isempty( filenumber ) && filenumber >= 23
            layout = 549;
        end
        % channel allocation
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        switch layout
            case 540
                data{  2 } = {   33 : 37, 'type', 'stim' };
                data{  3 } = {        38, 'type', 'stimdead' };
                data{  5 } = {        40, 'type', 'temperature' };
            case 549
                data{  2 } = {   33 : 38, 'type', 'stim' };
                data{  5 } = {        40, 'type', 'theta' };
        end
        data{  4 } = {        39, 'type', 'x' };
        data{  6 } = {   41 : 43, 'type', 'am' };
        data{  7 } = {        44, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 32, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = { 33 : 38, 'voltageRange', '0.2' };    % gain (from Iout to Vout of the CS) = 50; range (of digitization of AIS in BAD)=10; range/gain=0.2 (note x50 gain at CS)
        data{ 13 } = { 39 : 40, 'voltageRange', '10' };     % amp=1; range=10
        data{ 14 } = { 41 : 43, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 15 } = {      44, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      33, 'target', '1' };
        data{ 22 } = {      34, 'target', '2' };
        data{ 23 } = {      35, 'target', '3' };
        data{ 24 } = {      36, 'target', '4' };
        data{ 25 } = {      37, 'target', '1' };
        data{ 26 } = {      38, 'target', '4' };
        
        data{ 31 } = { [ 33 36 : 38 ], 'source', 'LD' };
        data{ 32 } = { [ 34 35 ], 'source', 'LED' };
        data{ 41 } = { [ 33 36 ], 'wavelength', '450' };  % [nm]
        data{ 42 } = { [ 34 35 ], 'wavelength', '470' };  % [nm]
        data{ 43 } = { [ 37 38 ], 'wavelength', '638' };  % [nm]
        
        data{ 51 } = {      33, 'power', '202' };         % [uW]
        data{ 52 } = {      34, 'power', '41' };         % [uW]
        data{ 53 } = {      35, 'power', '41' };         % [uW]
        data{ 54 } = {      36, 'power', '103' };       % [uW]
        data{ 55 } = {      37, 'power', '1300' };      % [uW]
        data{ 56 } = {      38, 'power', '940' };       % [uW]
        data{ 61 } = { 33 : 38, 'voltage', '5' };       % [V]
        
    case 'mC41' % 6-shank Stark64 probe with 6 LEDs
        % channel allocation
        % shanks 1, 6:      10 sites
        % shanks 2,3,4,5:   11 sites
        % in practice:
        % shank 2 - broken
        % shanks 3 and 4 - top site high impedance
        % thus only 51 sites (10 per shank, shank "4" (physically 5) with 11
        % from now on ignore S#2 and channel 53 and treat S#3 as S#2 etc
        % thus, this is effectively a FIVE SHANK PROBE
        data{  1 } = {    1 : 51, 'type', 'neuronal' };
        data{  2 } = {   52 : 57, 'type', 'stim' };         % 6 CS channels (52: S#1; 54-57: S#2-5)
        data{  3 } = {        58, 'type', 'x' };            % MC DAC.0
        data{  4 } = {        59, 'type', 'theta' };        % MC DAC.3
        data{  5 } = {   60 : 62, 'type', 'am' };
        data{  6 } = {        63, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 51, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = { 52 : 57, 'voltageRange', '1' };      % amp=10 (gain x10 in CS); range=10 (of AIS inside BAD)
        data{ 13 } = { 60 : 62, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      63, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      52, 'target', '1' };
        data{ 22 } = {      53, 'target', '6' };
        data{ 23 } = {      54, 'target', '2' };
        data{ 24 } = {      55, 'target', '3' };
        data{ 25 } = {      56, 'target', '4' };
        data{ 26 } = {      57, 'target', '5' };
        
        data{ 31 } = { 52 : 57, 'source', 'LED' };
        
        data{ 41 } = { 52 : 55, 'wavelength', '470' };  % [nm]
        data{ 42 } = { 	    56, 'wavelength', '365' };  % [nm]
        data{ 43 } = { 	    57, 'wavelength', '470' };  % [nm]
        
        data{ 51 } = {      52, 'power', '39' };        % [uW]
        data{ 52 } = {      53, 'power', '46' };        % [uW]
        data{ 53 } = {      54, 'power', '37' };        % [uW]
        data{ 54 } = {      55, 'power', '43' };        % [uW]
        data{ 55 } = {      56, 'power', '78' };        % [uW]
        data{ 56 } = {      57, 'power', '11.6' };      % [uW]
        
        data{ 61 } = { 52 : 57, 'voltage', '5' };       % [V]
        
    case 'm588' % Edge32-100 linear/2-diode probe
    case 'm590' % 2-shank/6-uLED probe
    case 'm637' % 4-shank/4-LED probe
        % channel allocation
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        data{  2 } = {        33, 'type', 'x' };
        data{  3 } = {        34, 'type', 'y' };
        data{  4 } = {        35, 'type', 'theta' };
        data{  5 } = {        36, 'type', 'am' };
        data{  6 } = {   37 : 40, 'type', 'stim' };
        data{  7 } = {        41, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 32, 'voltageRange', '2.45' };
        data{ 12 } = { 33 : 36, 'voltageRange', '1' };      % amp=10 (gain x10 in CS); range=10 (of AIS inside BAD)
        data{ 13 } = { 37 : 40, 'voltageRange', '10' };
        data{ 14 } = { 	    41, 'voltageRange', '1' };
        
        % stimulation channels
        data{ 21 } = {      37, 'target', '1' };
        data{ 22 } = {      38, 'target', '2' };
        data{ 23 } = {      39, 'target', '3' };
        data{ 24 } = {      40, 'target', '4' };
        
        data{ 31 } = { 37 : 40, 'source', 'LED' };
        data{ 41 } = { 37 : 40, 'wavelength', '470' };  % [nm]
        data{ 51 } = {      37, 'power', '8' };         % [uW]
        data{ 52 } = {      38, 'power', '6.8' };       % [uW]
        data{ 53 } = {      39, 'power', '10.8' };      % [uW]
        data{ 54 } = {      40, 'power', '7.5' };       % [uW]
        data{ 61 } = { 37 : 40, 'voltage', '5' };       % [V]
        
    case 'm746' % 3-tetrode/1-LED probe
end


dat.daq = daq;
switch daq
    case { 'Intan', 'BAD' }
        dat.nBits = 16;
        dat.samplingRate = 20000;
        dat.voltageRange = 2.45;
        dat.amplification = 192;
end

return

% EOF