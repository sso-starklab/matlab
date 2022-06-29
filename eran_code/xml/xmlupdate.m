% xmlupdate         add details for each channel
%
% rxml = xmlupdate( infile, data, outfile )
%
% where data is a cell array, with each element having the following
% structure:
% { channel_numbers, tagname, value }
%
% infile        full filename
% data          cell array, with each element having the following structure:
%                   { channel_numbers, tagname, value }
% outfile       {''}, then does not save
%
% does
% (1) gives a number to each anatomical group (0-based)
% (2) changes skip status of all channels not in data to 0
% (3) adds/modifies children in the anatomicalDescription child to match
% the requested data
% 
% note
% this routine is specific to the xml structure used in ndmanager. thus
% different source (filebase.xml) and target () files should be used, since
% neuroscope will otherwise overwrite all changes
%
% calls         xmltools, xmlfindchild
%
% see also      xmlcopy, xmlgetchans, xmlgetcolors, xmlupdatecolors

% example       update channel types for a buzsaki32sp diode-probe w/ linear track behavior
%
% data{ 1 } = {  1 : 32, 'type', 'neuronal' }; 
% data{ 2 } = { 33 : 36, 'type', 'stim' }; 
% data{ 3 } = { 41 : 43, 'type', 'sensor' };
% data{ 4 } = { 45 : 46, 'type', 'ASD' }; 
% data{ 5 } = { 47, 'type', 'trig' }; 
% data{ 6 } = { 49 : 52 , 'type', 'solenoid' }; 
% data{ 7 } = { 53, 'type', 'sync' }; 
% data{ 8 } = { 54 : 56, 'type', 'am' };
% 
% rxml = xmltools( infile )
% rxml = xmlupdate( rxml, data );
% xmltools( rxml )

% 02-jan-13 ES

% revisions
% 07-jan-13 make sure all channels have type field

function rxml = xmlupdate( infile, data, outfile )

vflag = 0;
rxml = [];

mfname = upper( mfilename );
nargs = nargin;
if nargs < 2 || isempty( infile ) || isempty( data )
    error( 'missing data' )
end
if nargs < 3 || isempty( outfile )
    outfile = '';
end 

% load
if isa( infile, 'struct' )
    rxml = infile;
else
    verb( sprintf( '%s: Loading %s...', mfname, infile ), vflag );
    rxml = xmltools( infile );
end

% get the channel group field and map channel numbers (0-based) to fields
verb( sprintf( '%s: Mapping file...', mfname ), -vflag );
i1 = xmlfindchild( rxml.child( 2 ), 'anatomicalDescription' );
if ~isempty( i1 )
    i2 = xmlfindchild( rxml.child( 2 ).child( i1 ), 'channelGroups' );
end
if isempty( i1 ) || isempty( i2 )
    error( 'xml file structure differs' )
end
k = 0;
mod = 0;
map = [];
ngroups = length( rxml.child( 2 ).child( i1 ).child( i2 ).child );
for gi = 1 : ngroups
    nchans = length( rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child ); % get the number of channels
    rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).value = sprintf( '%d', gi - 1 ); % update the group number
    for ci = 1 : nchans
        k = k + 1;
        str = rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).value;
        map( k, : ) = [ str2double( str ) + 1 gi ci ];
        tmpattribs = rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).attribs( 1 );
        if strcmp( tmpattribs.name, 'skip' ) && ~isequal( tmpattribs.value, '1' )
            mod = 1;
            rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).attribs( 1 ).value = '1';
        end
    end
end

% updte the requested fields
verb( sprintf( 'Updating fields...' ), -vflag );
skip = 0;
for ti = 1 : length( data ) % go over data fields
    
    if length( data{ ti } ) ~= 3
        continue
    end
    chans = data{ ti }{ 1 };
    newchild.tag = data{ ti }{ 2 };
    newchild.attribs.name = ''; 
    newchild.attribs.value = ''; 
    newchild.value = data{ ti }{ 3 }; 
    newchild.child = [];
    
    idx = map( ismember( map( :, 1 ), chans ), 2 : 3 );
    for i = 1 : size( idx, 1 ) % go over channels
        gi = idx( i, 1 );
        ci = idx( i, 2 );
        % add a child to each channel in the anatomical groups (same idea as used to give details about the units)
        nchilds = length( rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child );
        for j = 1 : nchilds
            tmpchild = rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child( j );
            if isequal( tmpchild, newchild )
                % an identical child already exists at that node
                skip = 1;
                break
            elseif isequal( tmpchild.tag, newchild.tag ) && isequal( tmpchild.attribs, newchild.attribs )
                % this is an update of the value of an existing field
                rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child( j ).value = newchild.value;
                skip = 1;
                mod = 1;
                break
            end
        end
        if skip
            skip = 0;
            break
        end
        % no such child exists, append
        mod = 1;
        if nchilds
            rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child( nchilds + 1 ) = newchild;
        else
            rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child = newchild;
        end
        if strcmp( rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).attribs( 1 ).name, 'skip' )
                   rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).attribs( 1 ).value = '0';
        end
    end
end

% make sure all channels have 'type' children
newchild.tag = 'type';
newchild.attribs.name = '';
newchild.attribs.value = '';
newchild.value = '';
newchild.child = [];
for ti = 1 : size( map, 1 )
    gi = map( ti, 2 );
    ci = map( ti, 3 );
    skip = 0;
    % add a child to each channel in the anatomical groups (same idea as used to give details about the units)
    nchilds = length( rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child );
    for j = 1 : nchilds
        tmpchild = rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child( j );
        if isequal( tmpchild.tag, newchild.tag )
            % an identical child already exists at that node
            skip = 1;
            break
        end
    end
    if skip
        continue
    end
    if nchilds
        rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child( nchilds + 1 ) = newchild;
    else
        rxml.child( 2 ).child( i1 ).child( i2 ).child( gi ).child( ci ).child = newchild;
    end
end

% summarize and save
if ~mod
    verb( sprintf( 'All fields already updated!\n' ), -vflag );
end
if ~isempty( outfile ) && exist( fileparts( outfile ), 'dir' )
    verb( sprintf( 'Saving modified file...' ), -vflag );
    xmltools( rxml, outfile )
end
verb( sprintf( 'Done!' ), vflag );

return

% EOF

% for instance, if using mice w/ 4-shank diode-probe w/ linear-track
% behavior setup, the update should be:

% copying without the unnecessary data:
filebase = '/Volumes/My Passport/phaser7/mouse552/m552r1/dat/m552r1.005/m552r1.005';
children = { 'generalInfo', 'acquisitionSystem', 'fieldPotentials', 'anatomicalDescription', 'spikeDetection' };
insuffix = 'xml';
outsuffix = 'prm.xml';

% adding (updating) channel fields:
data{ 1 } = {  1 : 32, 'type', 'neuronal' }; 
data{ 2 } = { 33 : 36, 'type', 'stim' }; 
data{ 3 } = { 41 : 43, 'type', 'sensor' };
data{ 4 } = { 45 : 46, 'type', 'ASD' }; 
data{ 5 } = { 47, 'type', 'trig' }; 
data{ 6 } = { 49 : 52 , 'type', 'solenoid' }; 
data{ 7 } = { 53, 'type', 'sync' }; 
data{ 8 } = { 54 : 56, 'type', 'am' };

data{ 9 } =  { 33, 'target', 'shank 1' }; 
data{ 10 } = { 34, 'target', 'shank 2' }; 
data{ 11 } = { 35, 'target', 'shank 3' }; 
data{ 12 } = { 36, 'target', 'shank 4' };
data{ 13 } = { 33, 'source', 'LED' }; 
data{ 14 } = { 34, 'source', 'LED' }; 
data{ 15 } = { 35, 'source', 'LED' }; 
data{ 16 } = { 36, 'source', 'LED' };
data{ 17 } = { 33, 'wavelength', '470' }; 
data{ 18 } = { 34, 'wavelength', '470' }; 
data{ 19 } = { 35, 'wavelength', '470' }; 
data{ 20 } = { 36, 'wavelength', '470' }; 
data{ 25 } = { [ 41 49 50 ], 'location', 'left' };
data{ 26 } = { [ 42 51 52 ], 'location', 'right' };
data{ 27 } = { [ 43 ], 'location', 'center' };
data{ 28 } = { [ 50 52 ], 'reward', 'water' };
data{ 29 } = { [ 49 51 ], 'reward', 'sugar' };

data{ 21 } = { 1 : 32, 'voltageRange', '8' };
data{ 22 } = { 33 : 40, 'voltageRange', '2' };
data{ 23 } = { [ 41 : 47 49 : 53 ], 'voltageRange', '20' };
data{ 24 } = { 54 : 56, 'voltageRange', '8' };

% modifying existing fields:
pathtag1 = { 'acquisitionSystem' };
pathval1 = { '' };
tag1 = 'voltageRange';
value1 = '8';

% pathtag2 = { 'generalInfo' };
% pathval2 = { '' };
% tag2 = 'date';
% value2 = datestr(date, 'yyyy-mm-dd' );

% actual work:
rxml0 = xmltools( [ filebase '.' insuffix ] );                  % load the existing strcutre
rxml1 = xmlcopy( rxml0, children );                             % copy just some children
rxml2 = xmlupdate( rxml1, data );                               % add channel information
rxml3 = xmlmodify( rxml2, pathtag1, pathval1, tag1, value1 );   % modify some additional fields
%rxml3 = xmlmodify( rxml3, pathtag2, pathval2, tag2, value2 );
if strcmp( rxml3.child( 2 ).attribs( 1 ).name, 'creator' )      % modify the creator field
    rxml3.child( 2 ).attribs( 1 ).value = [ rxml3.child( 2 ).attribs( 1 ).value ' ' computer ]; 
end             
xmltools( rxml3 )                                               % display the final structure
xmltools(rxml3, [ filebase '.' outsuffix ] );                   % save it

par = LoadXml( [ filebase '.prm.xml' ] );                       % will contain all the 
% this is completely compatible with the LoadXml, have to input the 
% [ filebase '.prm' ] (or the full name with the xml suffix)
% the only thing missing will be the high-pass filtering parameter
% 
% does not have the neuroscope fields, this is intentional so could
% not be loaded by neuroscope and written over.
