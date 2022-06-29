% xmlmodify         modify/add a child to an xml file
%
% rxml = xmlmodify( infile, pathtag, pathvalue, tag, value, outfile )
%
% path is a cell array of child tag names, starting from child( 2 ) as the
% root. at the end of the path, the children are checked for the tag TAG;
% if such a tag is found, the value is updated; otherwise the child is
% appended
% 
% infile            either a full path or an xmltools strucutre
% pathtag,value     cell arrays
% tag, value        strings
% outfile           {''}; optional argument to save the modified structure
%
% calls             xmltools, xmlfindchild
%
% see also          xmlcopy, xmlgetchans, xmlgetcolors, xmlupdatecolors
%
% example
%
% pathtag = { 'acquisitionSystem' };
% pathval = { '' };
% tag = 'voltageRange';
% value = '8';
% rxml = xmlmodify( infile, pathtag, pathval, tag, value )

% 02-jan-12 ES


function rxml = xmlmodify( infile, pathtag, pathvalue, tag, value, outfile )

vflag = 0;
rxml = [];

mfname = upper( mfilename );
nargs = nargin;
if nargs < 5 || isempty( infile ) || isempty( pathtag ) || isempty( pathvalue ) || isempty( tag ) || isempty( value )
    error( 'missing data' )
end
if nargs < 6 || isempty( outfile )
    outfile = '';
end 
if ~isa( value, 'char' ) || ~isa( tag, 'char' )
    value = char( value );
    tag = char( tag );
    verb( sprintf( '%s: Note tag: %s; value: %s...', mfname, tag, value ), vflag );
end

% load
if isa( infile, 'struct' )
    rxml = infile;
else
    verb( sprintf( '%s: Loading %s...', mfname, infile ), vflag );
    rxml = xmltools( infile );
end

% get indices to requested mode
s = rxml.child( 2 );
ok = 1;
for i = 1 : length( pathtag ),
    idx = xmlfindchild( s, pathtag{ i } ); 
    for j = 1 : length( idx )
        if ~strcmp( s.child( idx( j ) ).value, pathvalue{ i } ), 
            idx( j ) = NaN; 
        end
    end
    idx( isnan( idx ) ) = [];
    if isempty( idx )
        ok = 0; 
        break
    end
    out( i ) = idx;
    s = s.child( out( i ) );
end

% update the relevant field
depth = length( out );
str = 'rxml.child( 2 )'; 
for i = 1 : depth
    str = sprintf( '%s.child( %d )', str, out( i ) ); 
end
s = eval( str );
for i = 1 : length( s.child )
    if strcmp( s.child( i ).tag, tag ) 
        idx = i; % only the first one right now
        break
    else
        idx = 0;
    end
end
if idx % child with requested tag exists, just modify the value
    str = sprintf( '%s.child( %d ).value = value;', str, idx );
else % add child
    newchild.tag = tag;
    newchild.attribs.name = '';
    newchild.attribs.value = '';
    newchild.value = value;
    newchild.child = [];
    str = sprintf( '%s.child( 1 ) = newchild;', str );
end
eval( str );

% save the new structure
if ~isempty( outfile ) && exist( fileparts( outfile ), 'dir' )
    verb( sprintf( 'Saving modified file...' ), -vflag );
    xmltools( rxml, outfile )
end

verb( sprintf( '%s: Done!', mfname ), vflag );

return

% EOF

% for instance, if using mice w/ 4-shank diode-probe w/ linear-track
% behavior setup, the update should be:
%     
% path = { 'acquisitionSystem', 'voltageRange' };
pathtag = { 'anatomicalDescription', 'channelGroups', 'group', 'channel' }; % group 1, channel 4
pathval = { '', '', '1', '4' };
tag = 'voltageRange'; % channel 4
value = '8'; % 8 V

pathtag = { 'acquisitionSystem' };
pathval = { '' };
tag = 'voltageRange';
value = '8';
rxml = xmlmodify( infile, pathtag, pathval, tag, value )
% i1 = xmlfindchild( rxml.child( 2 ), 'acquisitionSystem' ); i2 = xmlfindchild( rxml.child( 2 ).child( i1 ), 'voltageRange' ); rxml.child( 2 ).child( i1 ).child( i2 ).value = 8;

