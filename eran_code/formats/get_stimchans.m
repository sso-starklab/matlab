% get_stimchans             get list of channels that correspond to some type from (*prm.xml file)
%
% CALL          [ chans, target, voltagerange, source, wavelength, power, maxvoltage, par ] = get_stimchans( par, chans, type )
%
% GETS          par           par structure from *prm.xml file
%               chans         optional subset of channels to search in
%               type          {'stim'}, can be anything else ('sensor', 'sync'...)
%
% RETURNS       chans         channels carrying the stimulation signal
%               target        shank numbers
%               voltagerange  of the recording system
%               source        description (text), e.g. 'LED', 'LD', ...
%               wavelength
%               power
%               maxvoltage
%
% CALLS         LoadXml, get_egroup

% 03-jan-13 ES

% revisions
% 21-jan-13 extended support of non-stim channels
% 28-mar-13 added a post-hoc to convert negative channel numbers (probe
%               geometry) to egroups (functional anatomy)
% 28-may-13 wavelength and power output BEFORE par
% 16-nov-13 maxvoltage added to enable intensity computations from par
% 17-aug-19 cleaned up
% 31-aug-19 modified negative case 

function [ chans, target, voltagerange, source, wavelength, power, maxvoltage, par ] = get_stimchans( par, dchans, type, cnvrt )

chans                       = [];
target                      = [];
source                      = [];
voltagerange                = [];
wavelength                  = [];
power                       = [];
maxvoltage                  = [];

nargs                       = nargin;
if nargs < 1 || isempty( par )
    return
end
if nargs < 2 || isempty( dchans )
    dchans                  = [];
end
if nargs < 3 || isempty( type )
    type                    = 'stim';
end
if nargs < 4 || isempty( cnvrt )
    cnvrt                   = 1;
end

if ~isa( par, 'struct' ) && isa( par, 'char' ) && exist( fileparts( par ), 'dir' )
    par                     = LoadXml( [ par '.prm.xml' ] );
end
if ~isa( par, 'struct' )
    return
end

for i                   = 1 : length( par.AnatGrps )
    idx                 = find( ismember( par.AnatGrps( i ).Type, type ) );
    if sum( idx )
        % channel
        ichans          = par.AnatGrps( i ).Channels( idx ) + 1;
        if ~isempty( dchans )
            idx         = idx( ismember( ichans, dchans ) );
        end
        achan           = par.AnatGrps( i ).Channels( idx ) + 1;
        chans           = [ chans achan ];
        % voltage range
        voltagerange    = [ voltagerange par.AnatGrps( i ).VoltageRange( idx ) ];
        % voltage range
        maxvoltage      = [ maxvoltage par.AnatGrps( i ).MaxVoltage( idx ) ];
        % target
        if isempty( par.AnatGrps( i ).Target ) ...
                || ~isempty( idx ) && length(  par.AnatGrps( i ).Target ) < max( idx ) ...
                || isempty( par.AnatGrps( i ).Target( idx ) )
            if isempty( dchans ) || ~isempty( idx )
                target  = [ target NaN ];
            end
        else
            target      = [ target par.AnatGrps( i ).Target( idx ) ];
        end
        % source
        if isempty( par.AnatGrps( i ).Source ) ...
                || ~isempty( idx ) && length( par.AnatGrps( i ).Source ) < max( idx ) ...
                || isempty( par.AnatGrps( i ).Source( idx ) )
            if isempty( dchans ) || ~isempty( idx )
                source  = [ source {''} ];
            end
        else
            source  = [ source par.AnatGrps( i ).Source( idx ) ];
        end
        % wavelength
        if isempty( par.AnatGrps( i ).Wavelength ) ...
                || ~isempty( idx ) && length(  par.AnatGrps( i ).Wavelength ) < max( idx ) ...
                || isempty( par.AnatGrps( i ).Wavelength( idx ) )
            if isempty( dchans ) || ~isempty( idx )
                wavelength = [ wavelength NaN ];
            end
        else
            wavelength = [ wavelength par.AnatGrps( i ).Wavelength( idx ) ];
        end
        % power
        if isempty( par.AnatGrps( i ).Power ) ...
                || ~isempty( idx ) && length(  par.AnatGrps( i ).Power ) < max( idx ) ...
                || isempty( par.AnatGrps( i ).Power( idx ) )
            if isempty( dchans ) || ~isempty( idx )
                power = [ target NaN ];
            end
        else
            power   = [ power par.AnatGrps( i ).Power( idx ) ];
        end
        
        if isequal( chans, dchans )
            break
        end
    end
end

if ~strcmpi( type, 'stim' )
    return
end

% remove partial entries 
n                       = length( source );
ridx                    = false( n, 1 );
for i                   = 1 : n
    if isempty( source{ i } )
        ridx( i )       = 1;
    end
end
if sum( ridx )
    fprintf( '%s: Removing channels %s (no source designated)\n'...
        , upper( mfilename ), num2str( chans( ridx ) ) )
    chans( ridx )           = [];
    target( ridx )          = [];
    voltagerange( ridx )    = [];
    source( ridx )          = [];
    wavelength( ridx )      = [];
    power( ridx )           = [];
end
    
% post-hoc: convert channel numbers( indicated as negative targets ) to electrode groups
nidx                        = target < 0;
if sum( nidx ) && cnvrt
    nidx                    = find( nidx );
    targetG                 = zeros( 1, sum( nidx ) );
    for i = nidx( : ).'
        targetG( i )        = get_egroup( par, -target( nidx( i ) ) );
    end
    target                  = targetG;
end

return

% EOF
