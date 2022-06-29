% get_states            get state ranges
%
% [ periods, msg ]= get_states( filebase, statename )
%
% filebase          full path + base
% statename         the state name. any 3-character string will be checked.
%                       The standard format is: filebase.sts.statename 
%                       for instance *sts.the is theta, *sts.sws is
%                       slow-wave-sleep (no-theta, immobile animal)
%
% periods           [samples] 2-column matrix, at eeg sampling rate
% msg               error message if file not found
%
% see also          segmentBehavior (the generating routine)

% 17-mar-13 ES

% revisions
% 04-apr-13 extended support for (1) state 'all' (2) cell array of states

function [ periods, msg ] = get_states( filebase, statename )

if nargin < 2 || ~isa( filebase, 'char' ) || ~isa( statename, 'char' ) && ~isa( statename, 'cell' )
    error( 'Call: get_states( filebase, statename );' )
end

if isa( statename, 'cell' )
    periods = [];
    msg = '';
    for i = 1 : length( statename )
        if ~isa( statename{ i }, 'char' )
            msg = 'all cell array elements must be strings';
            return
        end
        [ per, ms ] = get_states( filebase, statename{ i } );
        periods = [ periods; per ];
        msg = [ msg ms ];
    end
    return
end
           
periods = [];
if ~exist( fileparts( filebase ), 'dir' )
    msg = 'missing filebase';
    return
end
statebase = [ filebase '.sts.' ];
if isempty( dir( [ statebase '*' ] ) )
    msg = 'no *sts* files';
    return
end

statename = lower( statename );
switch statename
    case { 'immobile' }
        suffix = 'imb';
    case { 'mobile' }
        suffix = 'mob';
    case { 'moving', 'move' }
        suffix = 'mov';
    case { 'notheta', 'nonetheta' }
        suffix = 'not';
    case { 'theta' }
        suffix = 'the';
    case { 'all' }
        suffix = 'all';
    otherwise
        if length( statename ) == 3 % 'rem', 'run', 'sws'
            suffix = statename;
        else
            msg = sprintf( 'unrecognized state %s', statename );
            return
        end
end

if strcmp( suffix, 'all' )
    eegfname = [ filebase '.eeg' ];
    info = dir( eegfname ); 
    if isempty( info )
        msg = sprintf( 'missing file %s', eegfname );
        return
    end
    par = LoadXml( filebase ); 
    periods = [ 1 info.bytes / 2 / par.nChannels ];
    msg = '';
    return
end

filename = [ filebase '.sts.' suffix ];
if ~exist( filename, 'file' )
    msg = sprintf( 'missing file %s', filename );
    return
end

periods = load( filename );
msg = '';

return

% EOF
