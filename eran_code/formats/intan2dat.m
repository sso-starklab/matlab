% intan2dat         convert Intan data to organized *dat files
% 
% call              [ newfilenames dirnames rc msg ] = intan2dat( pathToRawData, rawPrefix, varargin )
%
% gets              
% required:
%                   pathToRawData               valid path, e.g. '/Volumes/slab1/mC41/raw'. In it, directories starting with rawPrefix will be used
%                   rawPrefix                   common name, e.g. 'mC41_180819'. All directories with this prefix will be used
%
% optional (argument/value pairs):
%
%   for intan2dat:
%                   rawFilename                 {'all_in_one.dat'}  by default
%                   fromBFilename               {''}                if given (e.g. 'all_in_one.orig.dat'), the dat file (in pathToRawData) will be reconstructed from this file
%                   prefix                      {'m001'}            common prefix (e.g. mouse name)
%                   filenum1                    {1}                 how to call the first file in the set
%                   pathToDat                   {''}                if missing, will be created within pathToRawData
%
%   for massagedatfile:
%                   newOrder                    {[]}                argument passed to massagedatfile/reorderchannels. by default, channel order remains 
%                   chansToRemove               {[]}                argument passed to massagedatfile/removechannels. by default, all channels remain 
%                   nChannelsOrig               {[]}                either this or newOrder is required
%                   lineChannel                 {[]}                argument to removeline
%                   chansToClean                {[]}                argument to removeline
%
%   for mergefiles:
%                   rawPrefixSecond             {''}                similar to rawPrefix, e.g. if rawPrefix is 'mC400_1_200825' and rawPrefixSecond is 'mC400_2_200825',
%                                                                       the two files (from distinct directories) will be merged
%                   rawFilenameSecond           {'all_in_one.dat'}  similar to rawFilename (applies to file from a distinct directory)
%                   filenameSecond              {'null'}            applies to a file from the same directory as the first
%                                                                       e.g. 'digitalin.dat'
%                   nchansSecond                {0}                 applies to both cases
%                   precisionSecond             {'int16'}           applies to both cases
%
%                   sync1                       {[]}                argument to warpfile
%                   sync2                       {[]}                argument to warpfile
%                   shiftbits                   {[]}                argument to warpfile
%
% does:
%  This is a wrapper for massagedatfile
%  If mfilenum is provided, will also generate a merged (concatenated!) file in the pathToDat directory
%
% call example:
% 
%       pathToRawData	= '/Volumes/slab1/mC41/raw';
%       pathToDat       = '/Volumes/slab1/mC41/dat';
%       rawPrefix       = 'mC41_180819';
%       prefix          = 'mC41';
%       filenum1        = 149;
%       filenums        = [];
%       mfilenum        = 24;
%       newOrder        = [ 9    15    17     7    10    16    11     8    18     6    20    12     2     4    14     5    21    19    52     3     1    35    23    51    22    26    53    49    13    34    37    24    30    41    48    32    36   43    50    28    25    39    46    29    40    47    27    45    38    33    31    44    42    54    55    56    57    58    59    60    61    62    63    64    65 ];
%       chansToRemove 	= [ 1 24 ];
%       nChannelsOrig 	= 65;
% 
%       [ newfilenames dirnames rc msg ] = intan2dat( pathToRawData, rawPrefix, 'pathToDat', pathToDat, 'prefix', prefix, 'filenum1', filenum1, 'newOrder', newOrder, 'chansToRemove', chansToRemove, 'nChannelsOrig', nChannelsOrig, 'mfilenum', mfilenum );
%
%   In this recording from a Stark64 6-shank probe in which shank2 is missing,
%   two channels were deemed as not suitable for analysis (top channels in shanks 2 and 3).
%   Furthermore, the channels must be re-arranged
%   There were 53/6/2/3/1 neuronal/CS/MC/AM/digital channels
%   The first file is 149
%
% dependencies:
%   direct:         ParseArgPairs, num3str, massagedatfile, concatfiles
%   indirect:       removesamples, removechannels, reorderchannels
%                   removeline, removedcoffset, mergefiles

% 08-feb-18 ES

% revisions
% 25-jul-18 added chansToRemove
% 20-aug-18 added documentation
%           added fromBFilename
% 07-sep-18 modified makeBackup logic 
% 11-jul-19 added option to concatenate
% 26-aug-20 added support for merging second file from a distinct directory
% 30-aug-20 added support for warping second file from a distinct directory

% to do:
% (1) document code
% (2) add copying the *.001.xml to the first file (in the dat directory)        DONE 08-feb-18
% (3) add merging of two simultaneously-recorded files (intan1 and intan2)      DONE 26-aug-20

function [ newfilenames, dirnames, rc, msg ] = intan2dat( pathToRawData, rawPrefix, varargin )

% constants
suffix                              = 'dat';
newline                             = sprintf( '\n' );  % char( 9 ) is \t, not sure what is \n
precision                           = 'int16';

% initialize
newfilenames                        = cell( 1 );
dirname                             = cell( 1 );
rc                                  = cell( 1 );
msg                                 = cell( 1 );
t0                                  = clock;

% arguments
nargs = nargin;
if nargs < 2, error( 'pathToRawData and rawPrefix required' ), end
if ~exist( pathToRawData, 'dir' )
    msg{ 1 }                        = sprintf( 'missing directory %s', pathToRawData );
    rc{ 1 }                         = 1;
    if verbose
        fprintf( '%s\n', msg{ 1 } )
    end
    return
end
[ rawFilename, fromBFilename ...
    , prefix, filenum1, pathToDat, mfilenum ...
    , newOrder, chansToRemove, nChannelsOrig ...
    , sync1, sync2, shiftbits ...
    , lineChannel, chansToClean ...
    , rawPrefixSecond, rawFilenameSecond, filenameSecond, nchansSecond, precisionSecond ] = ParseArgPairs(...
    { 'rawFilename', 'fromBFilename' ...
    , 'prefix', 'filenum1', 'pathToDat', 'mfilenum' ...
    , 'newOrder', 'chansToRemove', 'nChannelsOrig' ...
    , 'sync1', 'sync2', 'shiftbits' ...
    , 'lineChannel', 'chansToClean' ...
    , 'rawPrefixSecond', 'rawFilenameSecond', 'filenameSecond', 'nchansSecond', 'precisionSecond' } ...
    , {  'all_in_one.dat', '' ...
    , 'm001', 1, '', [] ...
    , [], [], [] ...
    , [], [], [ 2^15 2^15 ] ...
    , [], [] ...
    , '', 'all_in_one.dat', 'null', 0, 'int16' } ...
    , varargin{ : } );
if ispc
    sep                             = '\';
    cpCmd                           = 'copy';
    rmCmd                           = 'del';
    mvCmd                           = 'rename';
    mkCmd                           = 'md';
else
    sep                             = '/';
    cpCmd                           = 'cp';
    rmCmd                           = 'rm';
    mvCmd                           = 'mv';
    mkCmd                           = 'mkdir';
end
if isempty( pathToDat )
    pathToDat = sprintf( '%s%sdat', pathToRawData, sep );
end
if ~exist( pathToDat, 'dir' )
    cmd                             = sprintf( '%s %s/%s', mkCmd, pathToDat );
    [ sts, txt ]                    = system( cmd );
    if sts
        fprintf( 1, 'NOTE: issue creating %s\n\t%s\n', pathToDat, txt )
    end
end
if isempty( newOrder ) && isempty( nChannelsOrig )
    error( 'must specify either newOrder or nChannelsOrig' )
end
    
% (1) get an aggolomerated char array with full paths to files (OS dependent)
if isunix
    root                            = [ sep strtok( pathToRawData, sep ) ];
    commonprefix                    = sprintf( '%s/%s', pathToRawData, rawPrefix );
    cmd                             = sprintf( 'ls -dl %s* | grep -o "%s.*"', commonprefix, root );
    [ sts, txt ]                    = system( cmd );
    if ~isempty( rawPrefixSecond )
        commonprefix2               = sprintf( '%s/%s', pathToRawData, rawPrefixSecond );
        cmd2                        = sprintf( 'ls -dl %s* | grep -o "%s.*"', commonprefix2, root );
        [ sts2, txt2 ]              = system( cmd2 );
    end
else
    error( 'must be used with a unix family OS (linux, mac, unix)' )
end
if sts || isempty( txt )
    error( 'check path: %s', commonprefix )
end
if ~isempty( rawPrefixSecond ) && ( sts2 || isempty( txt2 ) )
    error( 'check path: %s', commonprefix2 )
end

% (2) extract the individual directory names (OS independent)
nchars                              = length( txt );
i0                                  = zeros( nchars, 1 );
j                                   = 1;
i0( j )                             = 0;
for i                               = 1 : nchars
    if isequal( txt( i ), newline )
        j                           = j + 1;
        i0( j )                     = i;
    end
end
i1                                  = i0( 1 : ( j - 1 ) ) + 1;
i2                                  = i0( 2 : j ) - 1;
nfiles                              = j - 1;
dirnames                            = cell( nfiles, 1 );
for i                               = 1 : nfiles
    idx                             = i1( i ) : i2( i );
    dirnames{ i }                   = txt( idx );
end
filenums                            = filenum1 + ( 0 : nfiles - 1 );
% do the same for the second file (to be merged with the first, from a different directory)
if ~isempty( rawPrefixSecond )
    nchars                          = length( txt2 );
    i0                              = zeros( nchars, 1 );
    j                               = 1;
    i0( j )                         = 0;
    for i                           = 1 : nchars
        if isequal( txt( i ), newline )
            j                       = j + 1;
            i0( j )                 = i;
        end
    end
    i1                              = i0( 1 : ( j - 1 ) ) + 1;
    i2                              = i0( 2 : j ) - 1;
    nfiles2                         = j - 1;
    dirnames2                       = cell( nfiles2, 1 );
    for i = 1 : nfiles2
        idx                         = i1( i ) : i2( i );
        dirnames2{ i }              = txt2( idx );
    end
    filenums2                       = filenum1 + ( 0 : nfiles2 - 1 );
    if ~isequal( filenums, filenums2 )
        rawPrefixSecond             = ''; % cannot merge if directories are missing
    end
end

% (3) massage the dat files (OS independent):
nfiles                              = length( dirnames );
newfilenames                        = cell( nfiles, 1 );
for i                               = 1 : nfiles
    filename                        = sprintf( '%s%s%s', dirnames{ i }, sep, rawFilename );
    newfilename                     = sprintf( '%s%s%s.%s.%s', dirnames{ i }, sep, prefix, num3str( filenums( i ) ), suffix );
    makeBackup                      = 1;
    if ~isempty( fromBFilename )
        bfilename                   = sprintf( '%s%s%s', dirnames{ i }, sep, fromBFilename );
        if exist( bfilename, 'file' ) && ~exist( filename, 'file' )
            cmd                     = sprintf( '! %s %s %s', cpCmd, bfilename, filename );
            t0i                     = clock;
            fprintf( 1, 'Recovering %s from %s...', filename, bfilename )
            eval( cmd )
            et                      = etime( clock, t0i );
            fprintf( 1, 'done (%0.3g sec).\n', et )
            makeBackup              = 0;
        end
    end
    if isempty( rawPrefixSecond )
        filename2                   = filenameSecond;
    else
        filename2                   = sprintf( '%s%s%s', dirnames2{ i }, sep, rawFilenameSecond );
    end
    [ rcI, msgI ]                   = massagedatfile( filename, 'makeBackup', makeBackup, 'chansToRemove', chansToRemove ...
        , 'nChannelsOrig', nChannelsOrig, 'newOrder', newOrder, 'newfilename', newfilename ...
        , 'sync1', sync1, 'sync2', sync2, 'shiftbits', shiftbits ...
        , 'lineChannel', lineChannel, 'chansToClean', chansToClean ...
        , 'filenameSecond', filename2, 'nchansSecond', nchansSecond, 'precisionSecond', precisionSecond );
    rc{ 1 }{ i }                    = rcI;
    msg{ 1 }{ i }                   = msgI;
    newfilenames{ i }               = newfilename;
end

% (4) move all files to a common directory
tfilenames                          = cell( nfiles, 1 );
for i = 1 : nfiles
    sfilename                       = newfilenames{ i };
    [ ~, f1, f2 ]                   = fileparts( sfilename );
    tfilenames{ i }                 = sprintf( '%s%s%s%s', pathToDat, sep, f1, f2 );
    cmd                             = sprintf( '!%s %s %s', mvCmd, sfilename, tfilenames{ i } );
    eval( cmd )
end

% (5) concatenate files if requested 
if ~isempty( mfilenum )
    mfilename                       = sprintf( '%s%s%s_%d.%s', pathToDat, sep, prefix, mfilenum, suffix );
    if isempty( newOrder )
        nchans                      = nChannelsOrig;
    else
        nchans                      = length( newOrder ) - length( chansToRemove );
    end
    concatfiles( tfilenames, mfilename, nchans, precision );
end

% (6) copy *xml file if such exists 
sxmlfile                            = sprintf( '%s%s%s.001.xml', pathToDat, sep, prefix );
txmlfile                            = sprintf( '%s%s%s.%s.xml', pathToDat, sep, prefix, num3str( filenums( 1 ) ) );
if exist( sxmlfile, 'file' )
    cmd                             = sprintf( '!%s %s %s', cpCmd, sxmlfile, txmlfile );
    eval( cmd )
else
    %fprintf( 1, 'NOTE: missing source *xml file, cannot evaluate\n\t%s \n', cmd )
    fprintf( 1, 'NOTE: missing source *xml file, cannot copy\n' )
end

% (7) report
et                                  = etime( clock, t0 );
fprintf( 1, 'done massaging %d files (%0.3g sec).\n', nfiles, et( 1 ) )

return

% EOF


