% get_source_filenums       from *srs file
%
% [ fnums fnames ismerged ] = get_source_filenums( filebase )
%
% filebase          full path w/ base or par structure
%
% fnums             vector
% fnames            always a cell array
% ismerged          logical
%
% see also makesrslen

% 09-dec-12 ES

% revisions
% 25-feb-13 extended to support call with par structure (must have 'FileName' field)
% 25-sep-13 flexible format support
% 15-may-14 exception added..
% 17-aug-19 clean up

function [ fnums, fnames, ismerged ] = get_source_filenums( filebase )

fnums                               = [];
fnames                              = [];
ismerged                            = 0;

[ pathname, filename ]              = fileparts( filebase );
if isa( filebase, 'char' ) && exist( pathname, 'dir' )
elseif isa( filebase, 'struct' )
    if isfield( filebase, 'FileName' )
        filebase                    = filebase.FileName;
    else
        fprintf( '%s: missing filebase\n', upper( mfilename ) );
        return
    end
else
    fprintf( '%s: missing filebase\n', upper( mfilename ) );
    return
end

srsfname                            = sprintf( '%s.srs', filebase );
if ~exist( srsfname, 'file' )
    [ pathname, filename, extname ] = fileparts( filebase );
    if exist( pathname, 'dir' ) && ~isempty( extname )
        fnums                       = str2num( extname( 2 : end ) );
        fnames                      = { [ filename extname ] };
    elseif exist( [ filebase '.eeg' ], 'file' )
        fnums                       = 1;
        fnames                      = { filename };
        ismerged                    = 1;
    end
    return
end

ismerged                            = 1;
fid                                 = fopen( srsfname, 'r' );
filenames                           = fgetl( fid );
fclose( fid );
i = 0;
while ~isempty( filenames )
    i = i + 1;
    [ tok, filenames ]              = strtok( filenames );
    fnames{ i }                     = tok;
    [ bs, num ]                     = strtok( tok, '.' );
    if isempty( num )
        % file numbers are given by *-##
        [ mun, sb ]                 = strtok( fliplr( tok ), '-' ); 
        if length( mun ) > 9 && isequal( fliplr( mun( 1 : 8 ) ), '_stopped' )
            % exception..
            mun                     = mun( 9 : end );
        end
        fnums( i )                  = str2num( fliplr( mun ) );
    else % file numbers are given by *.###
        fnums( i )                  = str2num( num( 2 : end ) );
    end
end

return

% EOF
