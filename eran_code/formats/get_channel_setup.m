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
    
        case 'mA335' % 2 red LD
         % channel allocation
       if any (filenumber == [47,48,49])

        % mapping
        data{  2 } = {      1: 2, 'type', 'stim' };         
        data{  3 } = {        3, 'type', 'x' };            
        data{  4 } = {        4, 'type', 'y' };                 
        data{  5 } = {        5, 'type', 'theta' };
        data{  6 } = {        6, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 12 } = {  1 : 2, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = {  3 : 5, 'voltageRange', '10' };     %  range of digitization of AIS in BAD =10
        data{ 14 } = {      6, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
       else 
        data{  2 } = {      1 : 2, 'type', 'stim' };         
        data{  3 } = {      3 : 5, 'type', 'am' };            
        data{  4 } = {        6, 'type', 'x' };                 
        data{  5 } = {        7, 'type', 'y' };
        data{  6 } = {        8, 'type', 'theta' };
        data{  7 } = {        9, 'type', 'digitalin' };

        % voltage ranges
        data{ 12 } = {  1 : 2, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = {  3 : 5, 'voltageRange', '2.45' };     % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {  6 : 8, 'voltageRange', '10' };      % range of digitization of AIS in BAD =10
        data{ 15 } = {      9, 'voltageRange', '1' };      % digital units (no amp, 16bit range)

       end
        % stimulation channels
        data{ 21 } = { 1, 'target', '1' };           % channel 1
        data{ 31 } = { 1, 'source', 'LD' };
        data{ 41 } = { 1, 'wavelength', '638' };  % [nm]
        data{ 51 } = { 1, 'power', '6500' };        % [uW]
        data{ 61 } = { 1, 'voltage', '5' };       % [V]
        
        data{ 71 } = { 2, 'target', '2' };           % channel 2
        data{ 81 } = { 2, 'source', 'LD' };
        data{ 91 } = { 2, 'wavelength', '638' };  % [nm]
        data{ 101 } = { 2, 'power', '6500' };        % [uW]
        data{ 111 } = { 2, 'voltage', '5' };       % [V]
        
    case 'mS303' % proccessing only intanL of 4-shank/12-uLED probe, 2 red LD on shanks 1,3, 2 red LED on shanks 2,4
                 % so the dat file oncludes 4-shank/4-uLED (8 analog channels, 3 aux (x,y,theta), 47 channels in total)
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        data{  2 } = {   33 : 36, 'type', 'stim' }; % high power channels (i.e. LDs and LEDs)
        data{  3 } = {   37 : 40, 'type', 'stim' }; % uLEDS recorded on intan L
        data{  5 } = {   44, 'type', 'x' };
        data{  6 } = {   45, 'type', 'y' };
        data{  7 } = {   46, 'type', 'theta' };
        data{  8 } = {   41 : 43, 'type', 'am' };
        data{  9 } = {        47, 'type', 'digitalin' };

        % voltage ranges
        data{ 10 } = {  1 : 32, 'voltageRange', '2.45' };     % amplification: 192
        data{ 11 } = { 33 : 36, 'voltageRange', '1' };        % gain (from Iout to Vout) = 10; range=10; range/gain=1 (assuming no x50 gain at CS_Rafi high power channels) 
        data{ 12 } = { 37 : 40, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
        data{ 14 } = { 41 : 43, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 15 } = {      47, 'voltageRange', '1' };      % digital units (no amp, 16bit range)

        % stimulation channels
        data{ 16 } = {      33, 'target', '1' };
        data{ 17 } = { 37 : 39, 'target', '1' };
        data{ 18 } = {      34, 'target', '2' };
        data{ 19 } = {      40, 'target', '2' };
        data{ 21 } = {      35, 'target', '3' };
        data{ 23 } = {      36, 'target', '4' };

        data{ 25 } = { 33, 'source', 'LD' };
        data{ 26 } = { 34, 'source', 'LED' };
        data{ 27 } = { 35, 'source', 'LD' };
        data{ 28 } = { 36, 'source', 'LED' };
        data{ 29 } = { 33, 'wavelength', '638' };
        data{ 30 } = { 34, 'wavelength', '617' };
        data{ 31 } = { 35, 'wavelength', '638' };
        data{ 32 } = { 36, 'wavelength', '617' };

        data{ 33 } = { 37 : 40, 'source', 'LED' };
        data{ 35 } = { 37 : 40, 'wavelength', '470' };  % [nm]

        data{ 40 } = {      37, 'power', '6.3' };       % [uW]
        data{ 41 } = {      38, 'power', '6.2' };       % [uW]
        data{ 42 } = {      39, 'power', '5.8' };       % [uW]
        data{ 43 } = {      40, 'power', '6' };       % [uW]           
        data{ 52 } = {      33, 'power', '6460' };       % [uW]    
        data{ 53 } = {      34, 'power', '24.9' };       % [uW]    
        data{ 54 } = {      35, 'power', '4450' };       % [uW]    
        data{ 55 } = {      36, 'power', '22.8' };       % [uW]    
        data{ 56 } = { 37 : 40, 'voltage', '4' };       % [V]
        data{ 58 } = { 33 : 36, 'voltage', '5' };       % [V] 
        
     case 'mA303' % 4-shank/12-uLED probe, 2 red LD on shanks 1,3, 2 red LED on shanks 2,4
        if filenumber == 1
        data{  1 } = {    1 : 31, 'type', 'neuronal' };
        data{  2 } = {   32 : 35, 'type', 'stim' }; % high power channels (i.e. LDs and LEDs)
        data{  3 } = {   36 : 39, 'type', 'stim' }; % uLEDS recorded on intan L
        data{  4 } = {   47 : 54, 'type', 'stim' }; % uLEDS recorded on intan R
        data{  5 } = {   43, 'type', 'x' };
        data{  6 } = {   44, 'type', 'y' };
        data{  7 } = {   45, 'type', 'theta' };
        data{  8 } = {   40 : 42, 'type', 'am' };
        data{  9 } = {        46, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 10 } = {  1 : 31, 'voltageRange', '2.45' };     % amplification: 192
        data{ 11 } = { 32 : 35, 'voltageRange', '1' };        % gain (from Iout to Vout) = 10; range=10; range/gain=1 (assuming no x50 gain at CS_Rafi high power channels) 
        data{ 12 } = { 36 : 39, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
        data{ 13 } = { 47 : 54, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)

        data{ 14 } = { 40 : 42, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 15 } = {      46, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 16 } = {      32, 'target', '1' };
        data{ 17 } = { 36 : 38, 'target', '1' };
        data{ 18 } = {      33, 'target', '2' };
        data{ 19 } = {      39, 'target', '2' };
        data{ 20 } = { 47 : 48, 'target', '2' };
        data{ 21 } = {      34, 'target', '3' };
        data{ 22 } = { 49 : 51, 'target', '3' };
        data{ 23 } = {      35, 'target', '4' };
        data{ 24 } = { 52 : 54, 'target', '4' };
        
        data{ 25 } = { 32, 'source', 'LD' };
        data{ 26 } = { 33, 'source', 'LED' };
        data{ 27 } = { 34, 'source', 'LD' };
        data{ 28 } = { 35, 'source', 'LED' };
        data{ 29 } = { 32, 'wavelength', '638' };
        data{ 30 } = { 33, 'wavelength', '617' };
        data{ 31 } = { 34, 'wavelength', '638' };
        data{ 32 } = { 35, 'wavelength', '617' };
        
        data{ 33 } = { 36 : 39, 'source', 'LED' };
        data{ 34 } = { 47 : 54, 'source', 'LED' };
        data{ 35 } = { 36 : 39, 'wavelength', '470' };  % [nm]
        data{ 36 } = { 47 : 54, 'wavelength', '470' };  % [nm]

        data{ 40 } = {      36, 'power', '6.3' };       % [uW]
        data{ 41 } = {      37, 'power', '6.2' };       % [uW]
        data{ 42 } = {      38, 'power', '5.8' };       % [uW]
        data{ 43 } = {      39, 'power', '6' };       % [uW]           
        data{ 44 } = {      47, 'power', '6' };       % [uW]
        data{ 45 } = {      48, 'power', '5.5' };       % [uW]
        data{ 46 } = {      49, 'power', '6.1' };       % [uW]
        data{ 47 } = {      50, 'power', '6.1' };       % [uW]
        data{ 48 } = {      51, 'power', '6.1' };       % [uW]
        data{ 49 } = {      52, 'power', '6.4' };       % [uW]
        data{ 50 } = {      53, 'power', '6.4' };       % [uW]
        data{ 51 } = {      54, 'power', '5.5' };       % [uW]    
        data{ 52 } = {      32, 'power', '6460' };       % [uW]    
        data{ 53 } = {      33, 'power', '24.9' };       % [uW]    
        data{ 54 } = {      34, 'power', '4450' };       % [uW]    
        data{ 55 } = {      35, 'power', '22.8' };       % [uW]    
        data{ 56 } = { 36 : 39, 'voltage', '4' };       % [V]
        data{ 57 } = { 47 : 54, 'voltage', '4' };       % [V] 
        data{ 58 } = { 32 : 35, 'voltage', '5' };       % [V] 
        else 
            data{  1 } = {    1 : 32, 'type', 'neuronal' };
            data{  2 } = {   33 : 36, 'type', 'stim' }; % high power channels (i.e. LDs and LEDs)
            data{  3 } = {   37 : 40, 'type', 'stim' }; % uLEDS recorded on intan L
            data{  4 } = {   48 : 55, 'type', 'stim' }; % uLEDS recorded on intan R
            data{  5 } = {   44, 'type', 'x' };
            data{  6 } = {   45, 'type', 'y' };
            data{  7 } = {   46, 'type', 'theta' };
            data{  8 } = {   41 : 43, 'type', 'am' };
            data{  9 } = {        47, 'type', 'digitalin' };

            % voltage ranges
            data{ 10 } = {  1 : 32, 'voltageRange', '2.45' };     % amplification: 192
            data{ 11 } = { 33 : 36, 'voltageRange', '1' };        % gain (from Iout to Vout) = 10; range=10; range/gain=1 (assuming no x50 gain at CS_Rafi high power channels) 
            data{ 12 } = { 37 : 40, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
            data{ 13 } = { 48 : 55, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)

            data{ 14 } = { 41 : 43, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
            data{ 15 } = {      47, 'voltageRange', '1' };      % digital units (no amp, 16bit range)

            % stimulation channels
            data{ 16 } = {      33, 'target', '1' };
            data{ 17 } = { 37 : 39, 'target', '1' };
            data{ 18 } = {      34, 'target', '2' };
            data{ 19 } = {      40, 'target', '2' };
            data{ 20 } = { 48 : 49, 'target', '2' };
            data{ 21 } = {      35, 'target', '3' };
            data{ 22 } = { 50 : 52, 'target', '3' };
            data{ 23 } = {      36, 'target', '4' };
            data{ 24 } = { 53 : 55, 'target', '4' };

            data{ 25 } = { 33, 'source', 'LD' };
            data{ 26 } = { 34, 'source', 'LED' };
            data{ 27 } = { 35, 'source', 'LD' };
            data{ 28 } = { 36, 'source', 'LED' };
            data{ 29 } = { 33, 'wavelength', '638' };
            data{ 30 } = { 34, 'wavelength', '617' };
            data{ 31 } = { 35, 'wavelength', '638' };
            data{ 32 } = { 36, 'wavelength', '617' };

            data{ 33 } = { 37 : 40, 'source', 'LED' };
            data{ 34 } = { 48 : 55, 'source', 'LED' };
            data{ 35 } = { 37 : 40, 'wavelength', '470' };  % [nm]
            data{ 36 } = { 48 : 55, 'wavelength', '470' };  % [nm]

            data{ 40 } = {      37, 'power', '6.3' };       % [uW]
            data{ 41 } = {      38, 'power', '6.2' };       % [uW]
            data{ 42 } = {      39, 'power', '5.8' };       % [uW]
            data{ 43 } = {      40, 'power', '6' };       % [uW]           
            data{ 44 } = {      48, 'power', '6' };       % [uW]
            data{ 45 } = {      49, 'power', '5.5' };       % [uW]
            data{ 46 } = {      50, 'power', '6.1' };       % [uW]
            data{ 47 } = {      51, 'power', '6.1' };       % [uW]
            data{ 48 } = {      52, 'power', '6.1' };       % [uW]
            data{ 49 } = {      53, 'power', '6.4' };       % [uW]
            data{ 50 } = {      54, 'power', '6.4' };       % [uW]
            data{ 51 } = {      55, 'power', '5.5' };       % [uW]    
            data{ 52 } = {      33, 'power', '6460' };       % [uW]    
            data{ 53 } = {      34, 'power', '24.9' };       % [uW]    
            data{ 54 } = {      35, 'power', '4450' };       % [uW]    
            data{ 55 } = {      36, 'power', '22.8' };       % [uW]    
            data{ 56 } = { 37 : 40, 'voltage', '4' };       % [V]
            data{ 57 } = { 48 : 55, 'voltage', '4' };       % [V] 
            data{ 58 } = { 33 : 36, 'voltage', '5' };       % [V] 
        end
    case 'mS51' %128 channels, 1 LD, 1 LED
         % channel allocation
        
        % mapping
        data{  1 } = {    1 : 126, 'type', 'neuronal' };
        data{  2 } = {   127 : 128, 'type', 'stim' };         %
        data{  3 } = {        129, 'type', 'x' };            %
%        data{  4 } = {        34, 'type', 'theta' };            %     
        data{  5 } = {   130 : 132, 'type', 'am' };
        data{  6 } = {        133, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 126, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = {  127 : 128, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 130 : 132, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      133, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = { 127, 'target', '1' };           % channel 127
        data{ 31 } = { 127, 'source', 'LED' };
        data{ 41 } = { 127, 'wavelength', '470' };  % [nm]
        data{ 51 } = { 127, 'power', '15' };        % [uW]
        data{ 61 } = { 127, 'voltage', '5' };       % [V]
        
        data{ 71 } = { 128, 'target', '8' };           % channel 128
        data{ 81 } = { 128, 'source', 'LD' };
        data{ 91 } = { 128, 'wavelength', '450' };  % [nm]
        data{ 101 } = { 128, 'power', '493' };        % [uW]
        data{ 111 } = { 128, 'voltage', '5' };       % [V]
        
    case 'mK01' % A1x32 Edge 20 um probe
        % channel allocation
        
        % mapping
        data{  1 } = {    1 : 31, 'type', 'neuronal' };
        data{  2 } = {   	  32, 'type', 'stim' };         %
        data{  3 } = {        33, 'type', 'x' };            %
        data{  4 } = {        34, 'type', 'theta' };            %     
        data{  5 } = {   35 : 37, 'type', 'am' };
        data{  6 } = {        38, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 11 } = {  1 : 31, 'voltageRange', '2.45' };   % amplification: 192
        data{ 12 } = {      32, 'voltageRange', '0.2' };    % gain (from Iout to Vout) = 50; range=10; range/gain=0.2 (assuming x50 gain at CS)
        data{ 13 } = { 35 : 37, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 14 } = {      38, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 21 } = {      32, 'target', '1' };           % channel 31
        data{ 31 } = {      32, 'source', 'LED' };
        data{ 41 } = {      32, 'wavelength', '470' };  % [nm]
        data{ 51 } = {      32, 'power', '32.3' };        % [uW]
        data{ 61 } = {      32, 'voltage', '5' };       % [V]
        
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
 
      case 'mC400' % 4-shank/12-uLED probe
          
        outanON = [2:11,13,15,17,19,21,23,25,26,27];
        data{  1 } = {    1 : 31, 'type', 'neuronal' };
        data{  2 } = {   32 : 35, 'type', 'stim' };
        data{  3 } = {   40 : 47, 'type', 'stim' };
        data{  4 } = {   36 : 38, 'type', 'am' };
        data{  5 } = {        39, 'type', 'digitalin' };
        

        % voltage ranges
        data{ 6 } = {  1 : 31, 'voltageRange', '2.45' };   % amplification: 192
        data{ 7 } = { 32 : 35, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
        data{ 8 } = { 40 : 47, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
        data{ 10 } = { 36 : 38, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 11 } = {      39, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        

        % stimulation channels
        data{ 12 } = { 32 : 34, 'target', '1' };
        data{ 13 } = {      35, 'target', '2' };
        data{ 14 } = { 40 : 41, 'target', '2' };
        data{ 15 } = { 42 : 44, 'target', '3' };
        data{ 16 } = { 45 : 47, 'target', '4' };
        

        data{ 17 } = { 32 : 35, 'source', 'LED' };
        data{ 18 } = { 40 : 47, 'source', 'LED' };
        data{ 19 } = { 32 : 35, 'wavelength', '470' };  % [nm]
        data{ 20 } = { 40 : 47, 'wavelength', '470' };  % [nm]
        
        if ismember (filenumber, outanON) 

            data{ 24 } = {      32, 'power', '7.72' };       % [uW]
            data{ 25 } = {      33, 'power', '7.72' };       % [uW]
            data{ 26 } = {      34, 'power', '7.72' };       % [uW]
            data{ 27 } = {      35, 'power', '7.72' };       % [uW]           
            data{ 28 } = {      40, 'power', '7.72' };       % [uW]
            data{ 29 } = {      41, 'power', '0' };       % [uW]
            data{ 30 } = {      42, 'power', '7.72' };       % [uW]
            data{ 31 } = {      43, 'power', '7.72' };       % [uW]
            data{ 32 } = {      44, 'power', '7.72' };       % [uW]
            data{ 33 } = {      45, 'power', '7.72' };       % [uW]
            data{ 34 } = {      46, 'power', '7.72' };       % [uW]
            data{ 35 } = {      47, 'power', '7.72' };       % [uW]    
            data{ 40 } = { 32 : 35, 'voltage', '5' };       % [V]
            data{ 41 } = { 40 : 47, 'voltage', '5' };       % [V] 
                       
        else
            data{ 24 } = {      29, 'power', '5.73' };     % [uW]
            data{ 25 } = {      30, 'power', '5.51' };         % [uW]
            data{ 26 } = {      31, 'power', '6.15' };         % [uW]
            data{ 27 } = {      32, 'power', '5.7' };       % [uW]
            data{ 28 } = {      37, 'power', '6.36' };       % [uW]
            data{ 26 } = {      38, 'power', '0' };       % [uW]
            data{ 27 } = {      39, 'power', '6.25' };       % [uW]
            data{ 28 } = {      40, 'power', '6.34' };       % [uW]
            data{ 29 } = {      41, 'power', '6.38' };       % [uW]
            data{ 30 } = {      42, 'power', '6.59' };       % [uW]
            data{ 31 } = {      43, 'power', '6.32' };       % [uW]
            data{ 32 } = {      44, 'power', '6.3' };       % [uW]
            data{ 40 } = { 32 : 35, 'voltage', '4' };       % [V]
            data{ 41 } = { 40 : 47, 'voltage', '4' };       % [V] 
        end
        
        if (filenumber ==23)
             data{  1 } = {    1 : 30, 'type', 'neuronal' };
            data{  2 } = {   31 : 34, 'type', 'stim' };
            data{  3 } = {   39 : 46, 'type', 'stim' };
            data{  4 } = {   35 : 37, 'type', 'am' };
            data{  5 } = {        38, 'type', 'digitalin' };   
            % voltage ranges
            data{ 6 } = {  1 : 30, 'voltageRange', '2.45' };   % amplification: 192
            data{ 7 } = { 31 : 34, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
            data{ 8 } = { 39 : 46, 'voltageRange', '0.001' };    % gain (from Iout to Vout of the CS) = 10000; range (of digitization of AIS in BAD)=10; range/gain=0.001 (note x50 gain at CS)
            data{ 10 } = { 35 : 37, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
            data{ 11 } = {      38, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
            % stimulation channels
            data{ 12 } = { 31 : 33, 'target', '1' };
            data{ 13 } = {      34, 'target', '2' };
            data{ 14 } = { 39 : 49, 'target', '2' };
            data{ 15 } = { 41 : 43, 'target', '3' };
            data{ 16 } = { 45 : 47, 'target', '4' };
            data{ 17 } = { 31 : 34, 'source', 'LED' };
            data{ 18 } = { 39 : 46, 'source', 'LED' };
            data{ 19 } = { 31 : 34, 'wavelength', '470' };  % [nm]
            data{ 20 } = { 39 : 46, 'wavelength', '470' };  % [nm]
            data{ 24 } = {      31, 'power', '7.72' };       % [uW]
            data{ 25 } = {      32, 'power', '7.72' };       % [uW]
            data{ 26 } = {      33, 'power', '7.72' };       % [uW]
            data{ 27 } = {      34, 'power', '7.72' };       % [uW]           
            data{ 28 } = {      39, 'power', '7.72' };       % [uW]
            data{ 29 } = {      40, 'power', '0' };       % [uW]
            data{ 30 } = {      41, 'power', '7.72' };       % [uW]
            data{ 31 } = {      42, 'power', '7.72' };       % [uW]
            data{ 32 } = {      43, 'power', '7.72' };       % [uW]
            data{ 33 } = {      44, 'power', '7.72' };       % [uW]
            data{ 34 } = {      45, 'power', '7.72' };       % [uW]
            data{ 35 } = {      46, 'power', '7.72' };       % [uW]    
            data{ 40 } = { 31 : 34, 'voltage', '5' };       % [V]
            data{ 41 } = { 39 : 46, 'voltage', '5' };       % [V] 
        end 
         case 'mDL5' % 4-shank/4-LD, 2-LED probe
        
        data{  1 } = {    1 : 32, 'type', 'neuronal' };

        data{  8 } = {   37:38 'type', 'stim' };
        data{  2 } = {   33:35 'type', 'stim' };
        data{  3 } = {   36            , 'type', 'stimdead'};
        data{  5 } = {        40, 'type', 'theta' };
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
        
        data{ 51 } = {      33, 'power', '68.3' };         % [uW]
        data{ 52 } = {      34, 'power', '45.5' };         % [uW]
        data{ 53 } = {      35, 'power', '39.7' };         % [uW]
        data{ 54 } = {      36, 'power', '0' };       % [uW]
        data{ 55 } = {      37, 'power', '604' };      % [uW]
        data{ 56 } = {      38, 'power', '1000' };       % [uW]
        data{ 61 } = { 33 : 38, 'voltage', '5' };       % [V]
        
          case 'mP101' % 4-shank 4-LED probe
        
        data{  1 } = {    1 : 30, 'type', 'neuronal' };
        data{  2 } = {   31:34 'type', 'stim' };
        data{  3 } = {        36, 'type', 'theta' };
        data{  4 } = {        35, 'type', 'x' };
        data{  5 } = {   37 : 39, 'type', 'am' };
        data{  6 } = {        40, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 7 } = {  1 : 30, 'voltageRange', '2.45' };   % amplification: 192
        data{ 8 } = { 31:34, 'voltageRange', '0.2' };    % gain (from Iout to Vout of the CS) = 50; range (of digitization of AIS in BAD)=10; range/gain=0.2 (note x50 gain at CS)
        data{ 9 } = { 35 : 36, 'voltageRange', '10' };     % amp=1; range=10
        data{ 10 } = { 37:39, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 11 } = {      40, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 12 } = {      31, 'target', '1' };
        data{ 13 } = {      32, 'target', '2' };
        data{ 14 } = {      33, 'target', '3' };
        data{ 15 } = {      34, 'target', '4' };       
 
        data{ 16 } = { 31:34, 'source', 'LED' };
        data{ 17 } = { 31:34, 'wavelength', '470' };  % [nm]
        data{ 18 } = {      31, 'power', '40.1' };         % [uW]
        data{ 19 } = {      32, 'power', '45.9' };         % [uW]
        data{ 20 } = {      33, 'power', '29.9' };         % [uW]
        data{ 21 } = {      34, 'power', '23.3' };       % [uW]
        data{ 22 } = { 31 : 34, 'voltage', '5' };       % [V]
        
        case 'mS234' % 4-shank 3-LED probe
        
        data{  1 } = {    1 : 28, 'type', 'neuronal' };
        data{  2 } = {   29:31 'type', 'stim' };
        data{  3 } = {        32, 'type', 'x' };
        data{  4 } = {        33, 'type', 'y' };
        data{  5 } = {   34 : 36, 'type', 'am' };
        data{  6 } = {        37, 'type', 'digitalin' };
        
        % voltage ranges
        data{ 7 } = {  1 : 28, 'voltageRange', '2.45' };   % amplification: 192
        data{ 8 } = { 29:31, 'voltageRange', '0.2' };    % gain (from Iout to Vout of the CS) = 50; range (of digitization of AIS in BAD)=10; range/gain=0.2 (note x50 gain at CS)
        data{ 9 } = { 32 : 33, 'voltageRange', '10' };     % amp=1; range=10
        data{ 10 } = { 34:36, 'voltageRange', '2.45' };   % amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 11 } = {      37, 'voltageRange', '1' };      % digital units (no amp, 16bit range)
        
        % stimulation channels
        data{ 12 } = {      29, 'target', '1' };
        data{ 13 } = {      30, 'target', '2' };
        data{ 14 } = {      31, 'target', '4' };
 
        data{ 16 } = { 29:31, 'source', 'LED' };
        data{ 17 } = { 29:31, 'wavelength', '470' };  % [nm]
        data{ 18 } = {      29, 'power', '5.3' };         % [uW]
        data{ 19 } = {      30, 'power', '21.2' };         % [uW]
        data{ 20 } = {      31, 'power', '16.8' };         % [uW]
        data{ 22 } = { 29 : 31, 'voltage', '5' };       % [V]
        
        case 'mA234'
        % channel allocation
        data{  1 } = {    1 : 32, 'type', 'neuronal' };
        data{  2 } = {   33 : 35, 'type', 'stim' };
        data{  5 } = {        36, 'type', 'x' };
        data{  4 } = {        37, 'type', 'theta' };
        data{  6 } = {   38 : 40, 'type', 'am' };
        data{  7 } = {        41, 'type', 'digitalin' };
       
        % voltage ranges
        data{ 11 } = {  1 : 32, 'voltageRange', '2.45' };   % neuronal: amplification: 192
        data{ 12 } = { 33 : 35, 'voltageRange', '0.2' };    % stimulation: gain (from Iout to Vout of the CS) = 50; range (of digitization of AIS in BAD)=10; range/gain=0.2 (note x50 gain at CS)
        data{ 13 } = { 36 : 37, 'voltageRange', '10' };     % theta and x: amp=1; range=10
        data{ 14 } = { 38 : 40, 'voltageRange', '2.45' };   % accelomter: amp, ; actual values are only 0.1 to 2.45V (Aux; on-chip accelerometer)
        data{ 15 } = {      41, 'voltageRange', '1' };      % digital: digital units (no amp, 16bit range)
        data{ 16 } = { 52 : 57, 'voltageRange', '1' };      % amp=10 (gain x10 in CS); range=10 (of AIS inside BAD)

        % stimulation channels
        data{ 21 } = {      33, 'target', '1' };
        data{ 22 } = {      34, 'target', '2' };
        data{ 23 } = {      35, 'target', '4' };
       
        data{ 32 } = {  33:35 , 'source', 'LED' };
        data{ 42 } = {  33:35 , 'wavelength', '470' };  % [nm]
       
        data{ 51 } = {      33, 'power', '5.3' };         % [uW]
        data{ 52 } = {      34, 'power', '21.2' };         % [uW]
        data{ 53 } = {      35, 'power', '16.8' };         % [uW]
        data{ 61 } = {  33:35, 'voltage', '5' };    

        
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
        
          case 'mB142' % 4-shank/4-LD, 2-LED probe
       
        data{  1 } = {    1 : 32, 'type', 'neuronal' };

        data{  8 } = {   37:38 'type', 'stim' };
        data{  2 } = {   33:36 'type', 'stim' };
        data{  5 } = {        40, 'type', 'theta' };
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
       
        data{ 51 } = {      33, 'power', '248' };         % [uW]
        data{ 52 } = {      34, 'power', '32.5' };         % [uW]
        data{ 53 } = {      35, 'power', '34.9' };         % [uW]
        data{ 54 } = {      36, 'power', '158' };       % [uW]
        data{ 55 } = {      37, 'power', '1600' };      % [uW]
        data{ 56 } = {      38, 'power', '1500' };       % [uW]
        data{ 61 } = { 33 : 38, 'voltage', '5' };       % [V]
        
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