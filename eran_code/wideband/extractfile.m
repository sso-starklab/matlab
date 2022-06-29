% EXTRACTFILE       extract a segment from a binary file and write into another file
%
% call              RC = EXTRACTFILE( SOURCEFILE, NEWFILE, PERIODS, NCHANS, PRECISION )
%
% gets              sourcefile              full path
%                   newfile                 full path
%                   periods                 two-element vector [samples]
%                   nchans                  integer
%                   precision               {'int16'}
%
% does              extract the segment from periods(1,1) to period(1,2) (in samples) of file sourcefile 
%                   and save it in newfile. both would have the same channel count and precision
%
% calls             nothing
%
% see also          PARTITION (to extract multiple segments), CONCATFILES (to merge back a partitioned file)
%
% 
% usage example:
% extract the 10 seconds starting at 0:02:00.000 to a new file; assume 20 kHz sampling rate and 56 channels
% >> cd /Volumes/Data/phaser4/mouse371/
% >> rc = extractfile( 'm371r2.034.dat', 'm371r2.034p2.dat', [ 2 * 60 * 20000 +  1 ( 2 * 60 + 10 ) * 20000 ], 56 )
% 

% 30-may-12 ES

% revisions
% 08-jul-20 cleaned up

function rc             = extractfile( sourcefile, newfile, periods, nchans, precision )

% input arguments
nargs                   = nargin;
if nargs < 3 || isempty( sourcefile ) || isempty( newfile ) || isempty( periods )
    error( 'missing input parameters' )
end
if sum( sum( periods ~= round( periods ) ) ) || sum( sum( periods <= 0 ) ) || ~ismember( size( periods, 2 ), [ 0 2 ] )
    error( 'periods should be a 2-element vector of non-negative integers' )
end
if nargs < 4 || isempty( nchans )
    nchans              = 32;
end
if nchans <= 0 || nchans ~= round( nchans )
    error( 'nchans should be a non-negative integer' )
end
if nargs < 5 || isempty( precision )
    precision           = 'int16';
end

% constants
BLOCKSIZE               = 2^20; % number of elements/block (not bytes)

% parse periods
durs                    = diff( periods, 1, 2 ) + 1; % only the first row of periods is used here
row                     = 1;
datasize                = [ nchans durs( row ) ];

% build the type casting string
precisionstr            = sprintf( '*%s', precision );

% determine number of bytes/sample/channel
a                       = ones( 1, 1, precision );
sourceinfo              = whos( 'a' );
nbytes                  = sourceinfo.bytes;

% divide into blocks
nelements               = prod( datasize );
nblocks                 = ceil( nelements / BLOCKSIZE );

% open files for reading and writing
fp0                     = fopen( sourcefile, 'r' );
if fp0 == -1
    error( 'fopen error' )
end
fp1                     = fopen( newfile, 'w' );
if fp1 == -1
    error( 'fopen error' )
end

% skip to starting position of sourcefile
startposition           = nbytes * nchans * ( periods( row, 1 ) - 1 ); % start at requested sample of requested channel
rc                      = fseek( fp0, startposition, 'bof' );
if rc
    error( 'fseek error' )
end

% go over the sourcefile in blocks and write the newfile
for bnum                = 1 : nblocks
    if bnum == nblocks
        toload          = nelements - ( nblocks - 1 ) * BLOCKSIZE;
    else
        toload          = BLOCKSIZE;
    end
    data1               = fread( fp0, toload, precisionstr );
    fwrite( fp1, data1, precision );
end

% close files
fclose( fp0 );
fclose( fp1 );

return

% EOF
