% CONCATFILES         concatenate multiple binary files
%
% RC = CONCATFILES( SOURCEFILES, NEWFILE, NCHANS, PRECISION )
%
% concatenate multiple binary files into a single file head to tail
% all files should have the same channel count and precision
%
% use example:
% cd /Volumes/Data/phaser4/mouse371/
% concatfiles( { 'm371r2.035p1.dat', 'm371r2.035p2.dat', 'm371r2.035p3.dat', 'm371r2.035p4.dat' }, 'm371r2.035.dat', 56 )
%
% see also: EXTRACTFILE, PARTITION

% 30-may-12 ES

% revisions
% 11-jul-19 added error if any sourcefile is missing

function concatfiles( sourcefiles, newfile, nchans, precision )

% constants
BLOCKSIZE       = 2^20; % number of elements/block (not bytes)

% input arguments
nargs       = nargin;
if nargs < 2 || isempty( sourcefiles ) || isempty( newfile ) || isempty( nchans )
    error( 'missing input parameters' )
end
if nchans <= 0 || nchans ~= round( nchans )
    error( 'nchans should be a non-negative integer' )
end
if nargs < 4 || isempty( precision )
    precision = 'int16';
end

% build the type casting string
precisionstr    = sprintf( '*%s', precision );

% determine number of bytes/sample/channel
a               = ones( 1, 1, precision );
sourceinfo      = whos( 'a' );
nbytes          = sourceinfo.bytes;

% open file for writing
fp1             = fopen( newfile, 'w' );
if fp1 == -1
    error( 'fopen error for output file (check permissions!)' )
end

% actually concatenate
for fnum        = 1 : length( sourcefiles )

    % check input file
    sourcefile = sourcefiles{ fnum };
    if ~exist( sourcefile, 'file' )
        error( 'missing file %s\n', sourcefile )
    end
    fileinfo    = dir( sourcefile );
    nelements   = floor( fileinfo( 1 ).bytes / nbytes );
    
    % open file for reading
    fp0 = fopen( sourcefile, 'r' );
    if fp0 == -1
        error( 'fopen error for input file %s (check permissions!)', sourcefile )
    end

    % go over the sourcefile in blocks and write to the newfile
    nblocks         = ceil( nelements / BLOCKSIZE );
    for bnum        = 1 : nblocks
        if bnum == nblocks
            toload  = nelements - ( nblocks - 1 ) * BLOCKSIZE;
        else
            toload  = BLOCKSIZE;
        end
        data1       = fread( fp0, toload, precisionstr );
        fwrite( fp1, data1, precision );
    end
    
    % close source file
    fclose( fp0 );
end

% close output file
fclose( fp1 );

return

% EOF
