% create_prm_file       create a *.prm.xml from an *.xml file
%
% par = create_prm_file( filebase, data, Overwrite, voltageRange )
%
% copy data from xml file (without the programs info), add relevant data,
% and save. return the data structure loaded from the saved file.
%
% compatibility:
% the format is completely compatible with the LoadXml, 
% the only thing missing is the high-pass filtering parameter
% 
% note:
% does not have the neuroscope fields, this is intentional so could
% not be loaded/overwritten by neuroscope.
% thus, neuroscope works with the *.xml file only
%
% example:
% [ data, animals ] = get_dmax_channel_setup( 'mouse552' );
% filebase = '/Volumes/My Passport/phaser7/mouse552/m552r1/dat/m552r1.005/m552r1.005';
% par = create_prm_file( filebase, data );
%
% calls         xmltools, xmlcopy, xmlmodify    (xml, generic)
%               xmlupdate                       (xml, format-specific)
%               LoadXml                         (blab)
%               verb                            (general)
% 
% see also      get_dmax_channel_setup

% 03-jan-13 ES

% revisions
% 20-aug-13 voltageRange argument
% 02-sep-13 voltageRange argument overloaded to dat structure w/ possible
%               fields such as 
%                <acquisitionSystem>
%                 <nBits>16</nBits>
%                 <nChannels>72</nChannels>
%                 <samplingRate>20000</samplingRate>
%                 <voltageRange>10</voltageRange>
%                 <amplification>400</amplification>
%                 <offset>0</offset>
% 23-feb-18 clean up

function [ par, xml ] = create_prm_file( filebase, data, Overwrite, voltageRange )

par                 = [];
xml                 = [];
vflag               = 1;

%--------------------------------------------------------------------%
% constants
%--------------------------------------------------------------------%
% copying without the unnecessary data (neuroscope, programs, multiFileProcessing)
children            = { 'generalInfo', 'acquisitionSystem', 'fieldPotentials', 'anatomicalDescription', 'spikeDetection' };
insuffix            = 'xml';
outsuffix           = 'prm.xml';

% modifying existing fields
pathtag1            = { 'acquisitionSystem' };
pathval1            = { '' };
% tag1 = 'voltageRange';

%--------------------------------------------------------------------%
% arguments
%--------------------------------------------------------------------%
mfname              = upper( mfilename );
nargs               = nargin;
if nargs < 1 || isempty( filebase )
    return
end
infile              = [ filebase '.' insuffix ];
outfile             = [ filebase '.' outsuffix ];
if ~exist( infile, 'file' )
    verb( sprintf( '%s: Missing input file %s', mfname, infile ), vflag )
    return
end
if nargs < 2 || isempty( data )
    return
end
if nargs < 3 || isempty( Overwrite )
    Overwrite       = -2;
end
if nargs < 4 || isempty( voltageRange )
    voltageRange    = 8;                        % the default range for the DataMAX system (16 bits)
end
if ~isempty( voltageRange ) && isa( voltageRange, 'struct' )
    dat             = voltageRange;
    fields          = fieldnames( dat );
    ridx            = ismember( fields, 'daq' );
    if sum( ridx )
        dat         = rmfield( dat, 'daq' );
        fields( ridx ) = [];
    end
    nfields         = length( fields );
    vals            = cell( nfields, 1 );
    for i = 1 : nfields
        vals{ i }   = num2str( getfield( dat, fields{ i } ) );
    end
else
    nfields         = 1;
    fields{ 1 }     = 'voltageRange';
    vals{ 1 }       = num2str( voltageRange );
end
%value1 = num2str( voltageRange );

if Overwrite <= 0 && exist( outfile, 'file' )
    if nargout > 0
        par         = LoadXml( outfile );
    end
    return
end
  
%--------------------------------------------------------------------%
% actual work:
%--------------------------------------------------------------------%
verb( sprintf( '%s: Converting file %s...', mfname, infile ), -vflag );
rxml0               = xmltools( infile );                                     % load the existing strcutre
rxml1               = xmlcopy( rxml0, children );                             % copy just some children
rxml2               = xmlupdate( rxml1, data );                               % add channel information
% rxml3 = xmlmodify( rxml2, pathtag1, pathval1, tag1, value1 );   % modify some additional fields
rxml3               = rxml2;
for i = 1 : nfields
    rxml3           = xmlmodify( rxml3, pathtag1, pathval1, fields{ i }, vals{ i } );   % modify some additional fields
end
if strcmp( rxml3.child( 2 ).attribs( 1 ).name, 'creator' )      % modify the creator field
    rxml3.child( 2 ).attribs( end + 1 ).name    = 'modifier';
    rxml3.child( 2 ).attribs( end ).value       = computer;
    rxml3.child( 2 ).attribs( end + 1 ).name    = 'date';
    rxml3.child( 2 ).attribs( end ).value       = datestr(date, 'yyyy-mm-dd' );
end
xml = rxml3;

% output
if exist( fileparts( outfile ), 'dir' )
    if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( outfile, 'file' ) )
        %verb( sprintf( '%s: Saving %s...', mfname, outfile ), -vflag );
        xmltools( xml, outfile );                                      % save the structure
        verb( sprintf( 'Saved *prm.xml file!!' ), vflag );
    else
        verb( sprintf( 'Done, not saving.' ), vflag );        
    end
else
    verb( sprintf( 'Done, cannot save!!' ), vflag );
end
par = LoadXml( outfile );                                       % output

return

% EOF

% batch procesing - go over all xml files in a given directory
% NOTE should modify the data structure above for different animals
% (or for different files when there are e.g. less LC fibers than DPSS
% sousrce, less analog channels than diodes, etc)

% example 1: all files in the same directory
setup = 'mouse552';
setup = 'mouse660';
setup = 'mouse260';
data = get_dmax_channel_setup( setup );
datdir = '/Volumes/My Passport/phaser7/mouse552/m552r1/dat';
datdir = '/Volumes/My Passport/phaser7/mouse660/m660r1/dat'; OW = 1;
datdir = '/Volumes/My Passport/phaser7/mouse260/m260r1/dat'; OW = 1;
dirs = dir( datdir );
for i = 1 : length( dirs )
    if dirs( i ).isdir && dirs( i ).name( 1 ) ~= '.'
        dirname = dirs( i ).name; 
        filebase = sprintf( '%s/%s/%s', datdir, dirname, dirname );
        par = create_prm_file( filebase, data, OW );
    end 
end

% example 2: separate directory for each recording date
animal = 'mouse482';
data = get_dmax_channel_setup( animal );
basedir = '/Volumes/Data/phaser4/mouse482';
subdirs = { '19mar12', '21mar12', '22mar12', '26mar12' };
for i = 1 : length( subdirs )
    datdir = sprintf( '%s/%s/dat', basedir, subdirs{ i } );
    dirs = dir( datdir );
    for j = 1 : length( dirs )
        if dirs( j ).isdir && dirs( j ).name( 1 ) ~= '.'
            dirname = dirs( j ).name;
            filebase = sprintf( '%s/%s/%s', datdir, dirname, dirname );
            par = create_prm_file( filebase, data );
        end
    end
end

% example 3: unknown directory structure, use get_animal_sessions
animal = 'mouse365';
animal = 'mouse428';
data = get_dmax_channel_setup( animal );
basedir = [ '/Volumes/Data/phaser3/' animal ];
subdirs = get_animal_sessions( animal );
%subdirs = { '19mar12', '21mar12', '22mar12', '26mar12' };
for i = 1 : length( subdirs )
    datdir = sprintf( '%s/%s/dat', basedir, subdirs{ i } );
    dirs = dir( datdir );
    for j = 1 : length( dirs )
        if dirs( j ).isdir && dirs( j ).name( 1 ) ~= '.'
            dirname = dirs( j ).name;
            filebase = sprintf( '%s/%s/%s', datdir, dirname, dirname );
            par = create_prm_file( filebase, data );
        end
    end
end

