% create_prm_files          for a directory and all source directories
% 
% par = create_prm_files( filebase, setup, Overwrite )
%
% filebase      full path
% setup         argument to get_dmax_channel_setup.m
% Overwrite     {-2}; passed to create_prm_file.m
%
% par           see create_prm_file for details
%
% see also      get_dmax_channel_setup, create_prm_file, LoadXml
%
% example:
% filebase = '/Volumes/lab/phaser8/mouse258/m258r1/dat/m258r1_78/m258r1_78';
% par = create_prm_files( filebase, 'm258r1' )

% 06-feb-13 ES

% revisions
% 02-apr-13 data modified for each file
% 23-feb-18 clean up
% 04-jul-19 mirror the file number of the merged file to get_channel_setup

function [ par, pars ] = create_prm_files( filebase, setup, Overwrite )

par                 = [];
pars                = [];
nargs               = nargin;
if nargs < 2 || isempty( filebase ) || isempty( setup )
    return
end
if nargs < 3 || isempty( Overwrite )
    Overwrite   = -2;
end
% get the number of the merged file (if any)
[ ~, fname ]        = fileparts( filebase );        % get the file name
idx                 = strfind( fname, '_' );        % search for underscore
fnum                = [];
if ~isempty( idx ) && length( fname ) > idx
    fnum            = fname( ( idx + 1 ) : end );
    fnum            = str2num( fnum );
end
%
[ data dat ]        = get_channel_setup( setup, fnum );
%[ data dat ]        = get_channel_setup( setup );
if isempty( data )
    return
end
par                 = create_prm_file( filebase, data, Overwrite, dat );
if isempty( par )
    par             = LoadXml( [ filebase '.prm.xml' ] );
end
[ fnums fnames ]    = get_source_filenums( filebase );
datdir              = fileparts( fileparts( filebase ) );
n = length( fnums );
if n < 2
    return
end
pars = repmat( par, [ 1 n ] );
for i = 1 : n
    fbase           = sprintf( '%s/%s/%s', datdir, fnames{ i }, fnames{ i } );
    dataI           = get_channel_setup( setup, fnums( i ) );
    if exist( fileparts( fbase ), 'dir' )
        pars( i )   = create_prm_file( fbase, dataI, Overwrite, dat );
    else
        fprintf( '%s: missing directory %s\n', upper( mfilename ), fileparts( fbase ) )
    end
end

return

% EOF
