% xmlcopy           copies xml file with the selected children
%
% rxml = xmlcopy( infile, children, outfile )
%
% infile        full path or xml strucutre
% children      to copy
% outfile       if empty, nothing is saved
%
% calls         xmlfindchild
%
% 02-jan-12 ES

function rxml = xmlcopy( infile, children, outfile )

vflag = 0;
mfname = upper( mfilename );

nargs = nargin;
if nargs < 1 || isempty( infile )
    rxml = '';
    return;
end
if nargs < 2 || isempty( children )
    children = '';
end
if nargs < 3 || isempty( outfile )
    outfile = '';
end

if isa( infile, 'struct' )
    rxml = infile;
else
    % load the file
    verb( sprintf( '%s: Loading %s...', mfname, infile ), vflag );
    rxml = xmltools( infile );
end

% get the channel group field
verb( sprintf( '%s: Mapping file...', mfname ), -vflag );
tokeep = false( length( rxml.child( 2 ).child ), 1 );
for i = 1 : length( children )
    idx = xmlfindchild( rxml.child( 2 ), children{ i } );
    if isnan( idx )
        continue
    end
    tokeep( idx ) = 1;
end

% save the new xml file
if sum( ~tokeep )
    verb( sprintf( 'Removing extra fields...' ), -vflag );
    rxml.child( 2 ).child( ~tokeep ) = [];
    if ~isempty( outfile ) && exist( fileparts( outfile ), 'dir' )
        verb( sprintf( 'Saving modified file...' ), -vflag );
        xmltools( rxml, outfile )
    end
end

verb( sprintf( 'Done!' ), vflag );

return

% EOF

% this is a generic function

% example (with the ndmanager format): 
% to generate an xml file with ONLY the generalInfo, acquisitioSystem,
% fieldPotentials, and anatomicalDescription. 

children = { 'generalInfo', 'acquisitionSystem', 'fieldPotentials', 'anatomicalDescription' };
insuffix = 'xml';
outsuffix = 'prm.xml';
rxml = xmlcopy( [ filebase '.' insuffix ], children, [ filebase '.' insuoutsuffix ] );


