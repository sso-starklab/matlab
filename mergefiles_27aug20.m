% mergefiles        two files that have exactly the same duration
%
% call              [ RC, MSG ] = MERGEFILES( SOURCEFILE1, SOURCEFILE2, NEWFILE, NCHANS, PRECISION, CASTMODE, CLIP2 )
%
% does              merge the data in SOURCEFILE1 and SOURCEFILE2
%                   and save it in a NEWFILE
% 
% note              the two channels must have the same duration, but may have a different
%                   number of channels: NCHANS( 1 ) and NCHANS( 2 )
%
% other arguments   PRECISION   is an input to fread/fwrite
%                               it is a cell array of 2 strings, default is { 'int16', 'uint16' }
%                               the NEWFILE will be of PRECISION{ 1 }
%
%                   CASTMODE is how to do type casting if the two files consist of different
%                               data types; relevant only for uint16->int16 case. the options are -1/0/{1}:
%                               -1:   the C/MATLAB default for uint16->int16: clip numbers above 2^15-1 (loses data!!)
%                               0     the "information preserving" option: translate (shift) to the -2^15 to 2^15-1 range (changes values!!)
%                               {1}:    the "intuitive" option: compress to the 0-2^15-1 range (lose resolution!!)
%
%                   CLIP2       {0}     if number of samples differs, do not do anything
%                               1       clips extra data from sourcefile2
%                                       (if sourcefile1 is longer, do not do anything)
%
% e.g. if file1, sample1 is:
%       c(1)t(1) c(2)t(1) .. c(n)t(1)
% and file2, sample1 is:
%       c(n+1)t(1) c(n+2)t(1) .. c(n+m)t(1)
% then fileNew, sample1 will simply be:
%       c(1)t(1) c(2)t(1) .. c(n)t(1) c(n+1)t(1) c(n+2)t(1) .. c(n+m)t(1)
% 
%
% calls             nothing
%
% see also          CONCATFILES, EXTRACTFILE, PARTITION (all work on the time dimension)
%                   SPLITFILE (works on the channel dimension)

% 17-mar-16 ES

% revisions
% 05-may-16 actually written the routine (previously just conceptualized)
% 03-dec-17 (1) initialized output; added msg
%           (2) modified handling of blocksize2
% 27-aug-20 (1) cleaned up
%           (2) implemented nbytes
%           (3) allowed merging if nsamples1 < nsamples2

function [ rc, msg ] = mergefiles( sourcefile1, sourcefile2, newfile, nchans, precision, castmode, clip2 )

% initialize output
rc                              = 0;
msg                             = '';

% constants
nfiles                          = 2;
BLOCKSIZE                       = 2^20; % number of elements/block (not bytes)

% input arguments
nargs                           = nargin;
if nargs < 3 || isempty( sourcefile1 ) || isempty( sourcefile2 ) || isempty( newfile )
    error( 'missing input parameters' )
end
if nargs < 4 || isempty( nchans ) || length( nchans ) ~= 2
    nchans                      = [ 32 4 ];
end
if any( nchans <= 0 ) || any( nchans ~= round( nchans ) )
    error( 'nchans should be non-negative integers' )
end
if nargs < 5 || isempty( precision ) || length( nchans ) ~= 2
    precision                   = { 'int16', 'uint16' };
end
if ~strcmp( precision{ 1 }, precision{ 2 } ) 
    if ~strcmp( precision{ 1 }, 'int16' ) && ~strcmp( precision{ 2 }, 'uint16' )
        error( 'unsupported' )
    end
end
if nargs < 6 || isempty( castmode )
    castmode                    = 1;
end
if ~ismember( castmode, [ -1 0 1 ] )
    error( 'unsupported' )
end
if nargs < 7 || isempty( clip2 )
    clip2                       = 0;
end

% build the type casting strings
precisionstr                    = cell( 1, nfiles );
for i = 1 : nfiles
    precisionstr{ i }           = sprintf( '*%s', precision{ i } );
end

% determine number of bytes/sample/channel
nbytes                          = zeros( nfiles, 1 );
for i                           = 1 : nfiles
    a                           = ones( 1, 1, precision{ i } );
    sourceinfo                  = whos( 'a' );
    nbytes( i )                 = sourceinfo.bytes;
end

% get the number of time samples/channel:
a1                              = dir( sourcefile1 );
nsamples1                       = a1.bytes / nchans( 1 ) / nbytes( 1 );
a2                              = dir( sourcefile2 );
nsamples2                       = a2.bytes / nchans( 2 ) / nbytes( 2 );
if ~isequal( nsamples1, round( nsamples1 ) )
    fprintf( 1, 'sourcefile1 incomplete\n' )
end
if ~isequal( nsamples2, round( nsamples2 ) )
    fprintf( 1, 'sourcefile2 incomplete\n' )
end
if nsamples1 ~= nsamples2
    if clip2 && ( nsamples1 < nsamples2 )
        msg                     = sprintf( 'file1: %d; file2: %d; file2 longer by %d samples' ...
            , nsamples1, nsamples2, nsamples2 - nsamples1 );
        fprintf( 1, '%s\n', msg )
    else
        msg                     = 'file duration/nchans mismatch';
        rc                      = -1;
        return
    end
end

% determine block sizes
blocksize                       = floor( BLOCKSIZE / max( nchans ) );
blocksize1                      = blocksize * nchans( 1 );
blocksize2                      = blocksize * nchans( 2 );

% divide into blocks
datasize1                       = [ nchans( 1 ) nsamples1 ];
datasize2                       = [ nchans( 2 ) nsamples1 ];
nelements1                      = prod( datasize1 );
nelements2                      = prod( datasize2 );
nblocks1                        = ceil( nelements1 / blocksize1 );
nblocks2                        = ceil( nelements2 / blocksize2 );
if nblocks1 ~= nblocks2
    error( 'debugging error' )
end
blocks1                         = [ 1 : blocksize1 : blocksize1 * nblocks1; blocksize1 : blocksize1 : blocksize1 * nblocks1 ]';
blocks1( nblocks1, 2 )          = nelements1;
blocks2                         = [ 1 : blocksize2 : blocksize2 * nblocks2; blocksize2 : blocksize2 : blocksize2 * nblocks2 ]';
blocks2( nblocks2, 2 )          = nelements2;
nblocks                         = nblocks1;

% open files for reading and writing
fp0                             = fopen( newfile, 'w' );
fp1                             = fopen( sourcefile1, 'r' );
fp2                             = fopen( sourcefile2, 'r' );
if fp0 == -1
    error( 'fopen error (%s)', newfile )
end
if fp1 == -1
    error( 'fopen error (%s)', sourcefile1 )
end
if fp2 == -1
    error( 'fopen error (%s)', sourcefile2 )
end

% go over the sourcefile in blocks and write the newfile
for bnum                        = 1 : nblocks

    % collect the data
    if bnum == nblocks
        toload1                 = nelements1 - ( nblocks - 1 ) * blocksize1;
        toload2                 = nelements2 - ( nblocks - 1 ) * blocksize2;
    else
        toload1                 = blocksize1;
        toload2                 = blocksize2;
    end
    data1                       = fread( fp1, toload1, precisionstr{ 1 } );
    data2                       = fread( fp2, toload2, precisionstr{ 2 } );
    
    % convert uint16 to int16 (specifically for uint16 in sourcefile2)
    if strcmp( precision{ 1 }, 'int16' ) && strcmp( precision{ 2 }, 'uint16' )
        switch castmode
            case -1 % 'clip': 0:2^15-1 remains; 2^15:2^16-1 are clipped to 2^16-1
                data2           = int16( data2 );
            case 0 % 'shift': 0:2^16-1 is mapped 1:1 to -2^15:2^16-1
                data2           = int16( double( data2 ) - 2^15 );
            case 1 %'compress': 0:2^16-1 are mapped to 0:2^15-1
                data2           = int16( floor( double( data2 ) / 2 ) );
        end
    end
    
    % merge
    data1hat                    = reshape( data1, nchans( 1 ), size( data1, 1 ) / nchans( 1 ) );
    data2hat                    = reshape( data2, nchans( 2 ), size( data2, 1 ) / nchans( 2 ) );
    data0hat                    = [ data1hat; data2hat ];
    data0                       = data0hat( : );
  
    % write out
    fwrite( fp0, data0, precision{ 1 } );
    
end

% close files
rc0                             = fclose( fp0 );
rc1                             = fclose( fp1 );
rc2                             = fclose( fp2 );
rc                              = [ rc0 rc1 rc2 ];

return

% EOF

