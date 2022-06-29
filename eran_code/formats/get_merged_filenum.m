% get_merged_filenum        and name
%
% CALL              [ mfilebase, mfname, mfnum ] = get_merged_filenum( filebase )
%
% GETS              filebase          full path + base or par structure
%
% RETURNS           mfnum             number (assuming an '_' notation)
%                   mfname            name
%                   mfilebase         full path + base
%
% CALLS             get_source_filenums

% 04-feb-13 ES

% revisions
% 25-feb-13 extended to support call with par structure (must have 'FileName' field)
% 20-mar-13 if no merged file, return the filebase
% 23-feb-18 clean up
% 17-aug-19 cleaned up

function [ mfilebase, mfname, mfnum ] = get_merged_filenum( filebase )

mfilebase                       = '';
mfname                          = '';
mfnum                           = NaN;

if isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
elseif isa( filebase, 'struct' )
    if isfield( filebase, 'FileName' )
        filebase                = filebase.FileName;
    else
        fprintf( '%s: missing filebase\n', upper( mfilename ) );
        return
    end
else
    fprintf( '%s: missing filebase\n', upper( mfilename ) );
    return
end

[ pathname, filename, extname ] = fileparts( fileparts( filebase ) );
filename                        = [ filename extname ];
dirs                            = dir( pathname );
for i = 1 : length( dirs )
    if ~dirs( i ).isdir || dirs( i ).name( 1 ) == '.'
        continue
    end
    srsbase                     = sprintf( '%s/%s/%s', pathname, dirs( i ).name, dirs( i ).name );
    [ ~, fnames, ismerged ]     = get_source_filenums( srsbase );
    if ismerged && ismember( filename, fnames )
        mfilebase               = srsbase;
        break
    end
end
if isempty( mfilebase )
    srsfname                    = dir( [ filebase '*srs' ] );
    if ~isempty( srsfname )
        [ ~, mfname ]           = fileparts( srsfname.name );
        mfilebase               = filebase;
    end
else
    [ ~, mfname ]               = fileparts( mfilebase );
end
if isempty( mfname )
    mfilebase                   = filebase;
    return
end
[ ~, rem ]                      = strtok( mfname, '_' );
if length( rem ) > 1
    mfnum                       = str2num( rem( 2 : end ) );
end

return

% EOF
