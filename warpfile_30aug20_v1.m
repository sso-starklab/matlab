% warpfile              time warps (offsets and squeezes) file2 with respect to file1
%
% call                  rc = warpfile( sourcefile1, sourcefile2, sync1, sync2, nchans1, nchans2, newfile )
%                       [ ..., msg, wf, os ] = warpfile( ..., precision, shiftbits, graphics, verbose, dryrun )
%
% gets                  sourcefile1         full path to reference file (file1)
%                       sourcefile2         full path to file to be warped (file2)
%                       sync1               sync signal in file1
%                       sync2               sync signal in file2
%                       nchans1             number of channels in file1
%                       nchans2             number of channels in file2
%                       newfile             full path to new (warped) file2
%
%                       precision           {'int16','int16'}
%                       shiftbits           {[ 2^15 2^15 ]}
%                       graphics            {0}
%                       verbose             {1}
%                       dryrun              {0}
%
% formats               sync channel        [ XX YY ] (e.g. [ 40 6 ] means bit 6 in channel 40; 1-based)
%                       shiftbits           2^15 shifts by 2^15 (from -32767 to 1; relevant to int16)
%
% does                  1. loads entire duration of sync signals
%                       2. determines negative to positive (and positive to negative) transitions 
%                       3. checks conditions for time warping:
%                           - same number of sync pulses
%                           - identical (up to 1 sample difference) duration of each pulse
%                           - file2 is longer than file1
%                           - warping factor smaller than 1 (no interpolation required)
%                       4. if conditions are permissive, time warps file2 with respect to file1
%                       5. verifies the warping by comparing the warped sync signal from file2
%
% note                  if number of sync signals differs between the files, this routine will exist
%                       if the duration of any sync event differs between files, this routine will exit
%                       if file2 is shorter, this routine will exit
%                       if the wf is larger than 1, this routine will exit
%
% calls                 nothing
%
% see also              mergefiles, parseDigital

% 28-aug-20 ES

% revisions
% 30-aug-20 (1) added wf and os to all output msg

% notes:
% (1) should add shiftbits in removeline
% (2) should make sure blocksize in other routines that do not require nchans multiples (e.g. mergefiles) is a power of 2
% blocksize2                      = 2.^floor( log2( floor( BLOCKSIZE / nchans2 ) * nchans2 ) );
% (3) should not allow merging (mergefiles) if number of samples differs

function [ rc, msg, wf, os ] = warpfile( sourcefile1, sourcefile2, sync1, sync2, nchans1, nchans2, newfile, precision, shiftbits, graphics, verbose, dryrun )

%------------------------------------------------------------------------
% initialize output
rc                              = [ 1 1 ];
msg                             = '';
wf                              = NaN;
os                              = NaN;
mfname                          = mfilename;

%------------------------------------------------------------------------
% constants
BITSBYTE                        = 8;            % 8-bit bytes
NFILES                          = 2;            % hard coded
MAXDIFF                         = 1;            % [samples]
BLOCKSIZE                       = 2^20;         % number of elements/block (not bytes)

%------------------------------------------------------------------------
% input arguments
nargs                           = nargin;
if nargs < 7 
    error( 'missing input parameters' )
end
if isempty( sourcefile1 ) || isempty( sourcefile2 )
    error( 'missing source file path and names' )
end
if isempty( sync1 ) || isempty( sync2 ) || length( sync1 ) < 2 || length( sync2 ) < 2
    error( 'missing sync signal numbers' )
end
if isempty( nchans1 ) || isempty( nchans2 )
    error( 'missing number of channels' )
end
if isempty( newfile )
    error( 'missing target file path and name' )
end

% precision
if nargs < 8 || isempty( precision )
    precision                   = { 'int16', 'int16' };
end
if isa( precision, 'char' )
    precision                   = repmat( { precision }, [ 1 NFILES ] );
end
if isa( precision, 'cell' ) && length( precision ) == 1
    precision                   = repmat( precision, [ 1 NFILES ] );
end
if ~strcmp( precision{ 1 }, precision{ 2 } ) 
    if ~strcmp( precision{ 1 }, 'int16' ) && ~strcmp( precision{ 2 }, 'uint16' )
        error( 'unsupported' )
    end
end

% shiftbits
if nargs < 9 || isempty( shiftbits )
    shiftbits                   = [ 2^15 2^15 ];
end
if length( shiftbits ) == 1
    shiftbits                   = repmat( shiftbits, [ 1 NFILES ] );
end

% graphics
if nargs < 10 || isempty( graphics )
    graphics                    = 0;
end
if nargs < 11 || isempty( verbose )
    verbose                     = 1;
end
if nargs < 12 || isempty( dryrun )
    dryrun                      = 0;
end

%------------------------------------------------------------------------
% collect information from arguments

% build the type casting strings
precisionstr                    = cell( 1, NFILES );
for i = 1 : NFILES
    precisionstr{ i }           = sprintf( '*%s', precision{ i } );
end

% determine number of bytes/sample/channel
nbytes                          = zeros( NFILES, 1 );
for i                           = 1 : NFILES
    a                           = ones( 1, 1, precision{ i } );
    sourceinfo                  = whos( 'a' );
    nbytes( i )                 = sourceinfo.bytes;
end

% check the channel and bit numbers
if nchans1 < 1 || nchans1 ~= round( nchans1 )
    error( 'nchans1 should be a non-negative integer' )
end
if nchans2 < 1 || nchans2 ~= round( nchans2 )
    error( 'nchans2 should be a non-negative integer' )
end
chan1                           = sync1( 1 );
bit1                            = sync1( 2 );
chan2                           = sync2( 1 );
bit2                            = sync2( 2 );
if chan1 > nchans1 || chan1 < 1
    error( 'chan1 should be a non-negative integer smaller than nchans1' )
end
if chan2 > nchans2 || chan2 < 1
    error( 'chan2 should be a non-negative integer smaller than nchans2' )
end
if bit1 > ( nbytes( 1 ) * BITSBYTE ) || bit1 < 1
    error( 'bit1 should be a non-negative integer smaller than %d', ( nbytes( 1 ) * BITSBYTE ) )
end
if bit2 > ( nbytes( 2 ) * BITSBYTE ) || bit2 < 1
    error( 'bit2 should be a non-negative integer smaller than %d', ( nbytes( 2 ) * BITSBYTE ) )
end
shiftbits1                      = shiftbits( 1 );
shiftbits2                      = shiftbits( 2 );

%------------------------------------------------------------------------
% collect information about the files
if verbose
	fprintf( 1, '%s: determining number of samples per file...', mfname );
end

% get the number of time samples/channel:
a1                              = dir( sourcefile1 );
nsamples1                       = a1.bytes / nchans1 / nbytes( 1 );
a2                              = dir( sourcefile2 );
nsamples2                       = a2.bytes / nchans2 / nbytes( 2 );
if ~isequal( nsamples1, round( nsamples1 ) )
    fprintf( 1, 'sourcefile1 incomplete\n' )
end
if ~isequal( nsamples2, round( nsamples2 ) )
    fprintf( 1, 'sourcefile2 incomplete\n' )
end
if nsamples1 ~= nsamples2
    if nsamples1 < nsamples2
        msg                     = sprintf( '%sfile1: %d; file2: %d; file2 is longer by %d samples\n' ...
            , msg, nsamples1, nsamples2, nsamples2 - nsamples1 );
        if verbose
            fprintf( 1, '%s', msg )
        end
    else
        msg                     = sprintf( '%sfile1: %d; file2: %d; file2 is shorter by %d samples (cannot be warped)\n' ...
            , msg, nsamples1, nsamples2, nsamples1 - nsamples2 );
        if verbose
            fprintf( 1, '%s', msg )
        end
        rc                      = -1;
        return
    end
end

if verbose && isempty( msg )
	fprintf( 1, 'file1: %d samples; file2: %d samples\n', nsamples1, nsamples2 );
end

%------------------------------------------------------------------------
% get the sync signal from each file and determine if warping by sync is possible
if verbose
	fprintf( 1, '%s: detecting sync events...', mfname );
end

m1                              = memmapfile( sourcefile1, 'Format', precision{ 1 }, 'writable', false );
d1                              = single( m1.Data( chan1 : nchans1 : end ) );
d1                              = rem( floor( ( d1 + shiftbits1 ) * pow2( 1 - bit1 ) ), 2 );
m2                              = memmapfile( sourcefile2, 'Format', precision{ 2 }, 'writable', false );
d2                              = single( m2.Data( chan2 : nchans2 : end ) );
d2                              = rem( floor( ( d2 + shiftbits2 ) * pow2( 1 - bit2 ) ), 2 );

% get the transitions
dvec                            = diff( d1 );
hi                              = find( dvec == 1 );
lo                              = find( dvec == -1 );
nhi                             = length( hi ); 
nlo                             = length( lo );
n                               = min( nhi, nlo );
mat1                            = [ hi( 1 : n ) lo( 1 : n ) ];

dvec                            = diff( d2 );
hi                              = find( dvec == 1 );
lo                              = find( dvec == -1 );
nhi                             = length( hi ); 
nlo                             = length( lo );
n                               = min( nhi, nlo );
mat2                            = [ hi( 1 : n ) lo( 1 : n ) ];

% check validity of syncs
nsyncs1                         = size( mat1, 1 );
nsyncs2                         = size( mat2, 1 );
if nsyncs1 ~= nsyncs2 
    msg                         = sprintf( '%sfile1: %d syncs; file2: %d syncs; cannot warp by syncs\n' ...
        , msg, nsyncs1, nsyncs2 );
    if verbose
        fprintf( 1, '%s', msg )
    end
    return
end
durs1                           = diff( mat1, [], 2 ) + 1;
durs2                           = diff( mat2, [], 2 ) + 1;
ddurs                           = abs( durs1 - durs2 );
maxdd                           = max( ddurs );
if maxdd > MAXDIFF
    msg                         = sprintf( '%smax difference: %d samples; cannot warp by syncs\n' ...
        , msg, maxdd );
    if verbose
        fprintf( 1, '%s', msg )
    end
    return
end

if verbose
	fprintf( 1, 'file1: %d syncs; file2: %d syncs\n', nsyncs1, nsyncs2 );
end

%------------------------------------------------------------------------
% determine warping factor and offset for file2 
if verbose
	fprintf( 1, '%s: determining warping factor and offset...', mfname );
end

dur1                            = mat1( n, 2 ) - mat1( 1, 1 );
dur2                            = mat2( n, 2 ) - mat2( 1, 1 );
os                              = mat2( 1, 1 ) - mat1( 1, 1 );              % offset (0: no offset)
%wf                              = nsamples1 / nsamples2;                   % could be derived from the number of samples - but then not synced
wf                              = dur1 / dur2;                              % warping factor (1: no warping)

if wf > 1
    msg                         = sprintf( '%swarping factor: %0.16g; offset: %d samples; this routine does not interpolate\n' ...
        , msg, wf, os );
    if verbose
        fprintf( 1, '%s', msg )
    end
    return
else
    msg                         = sprintf( '%swarping factor: %0.16g; offset: %d samples\n' ...
        , msg, wf, os );
end

% warp
t2                              = ( 1 : nsamples2 )';
[ ~, t2w ]                      = unique( round( t2 * wf ) );               % keep only a subset of samples
%t2w                             = round( interp1( t2, t2 * wf, 1 : nsamples1 ) ); % replicate some samples
%nw                              = length( t2w );                           % (nw/nsamples2/wf) will be the empirical wf after warping

% offset
% to impose (for debugging purposes only):
% os = -os; 
nans                            = [ 0 0 ];
if os > 0
    t2w( 1 : os )               = [];                                       % remove first os samples from file2
elseif os < 0
    t2w                         = [ NaN( -os, 1 ); t2w ];                    % add some null samples (NaN indices) at the beginning of file2
    nans( 1 )                   = -os;
end

% match samples with file1 
n2                              = length( t2w );
if n2 > nsamples1
    t2w( nsamples1 + 1 : n2 )   = [];                                       % remove last samples from file2
elseif n2 < nsamples1                                                           
    t2w                         = [ t2w; NaN( nsamples1 - n2, 1 ) ];        % append some null samples (NaNs indices) at the end of file2
    nans( 2 )                   = nsamples1 - n2;
end
% to impose (for debugging only):
% t2w( 1 : 100 ) = []; t2w = [ t2w; NaN( 100, 1 ) ]; os = t2w( 1 ); nans( 2 ) = 100;

if length( t2w ) ~= nsamples1
    error( 'internal (coding) bug' )
end

if verbose
	fprintf( 1, 'warping factor: %0.16g; offset: %d samples\n', wf, os );
end

%------------------------------------------------------------------------
% actually warp the sync signal of file2; plot

if graphics || dryrun

    if sum( nans ) > 0
        d2w                     = zeros( size( t2w ) );
        nnans                   = ~isnan( t2w );
        d2w( nnans )            = d2( t2w( nnans ) );
    else
        d2w                     = d2( t2w );
    end
    
end

if graphics
    
    nplot                       = 4 * 20000;
    if nsamples1 < nplot
        nplot                   = nsamples1;
    end
    xidx1                       = ( 1 : nplot );
    xidx2                       = length( d1 ) - ( nplot + 1 ) + ( 1 : nplot );

    figure
    subplot( 2, 3, 1 )
    plot( xidx1, d1( xidx1 ), 'b', xidx1, d2( xidx1 ), '--r' ), 
    title( 'before warping' )
    
    subplot( 2, 3, 2 )
    plot( xidx2, d1( xidx2 ), 'b', xidx2, d2( xidx2 ), '--r' )
    title( 'before warping' )
    
    subplot( 2, 3, 3 )
    plot( mat2( :, 1 ) - mat1( :, 1 ), 'b' )
    title( 'before warping' )
    xlabel( 'sync event' )
    ylabel( 'onset difference [samples]' )
    axis tight
    
    subplot( 2, 3, 4 )
    plot( xidx1, d1( xidx1 ), 'b', xidx1, d2w( xidx1 ), '--r' ), 
    title( 'after warping' )
    
    subplot( 2, 3, 5 )
    plot( xidx2, d1( xidx2 ), 'b', xidx2, d2w( xidx2 ), '--r' )
    title( 'after warping' )
    
    for i = 1 : 5
        subplot( 2, 3, i )
        set( gca, 'tickdir', 'out', 'box', 'off' );
    end
    
end

%------------------------------------------------------------------------
% prepare to warp, offset, and write out file2 (blockwise)
if verbose
	fprintf( 1, '%s: actually warping and offsetting...\n', mfname );
end

% determine block size
blocksize2                      = floor( BLOCKSIZE / nchans2 ) * nchans2;   % must be an integer multiple of the number of channels
nelements2                      = nsamples2 * nchans2;
nblocks                         = ceil( nelements2 / blocksize2 );

% partition into reading blocks
blocks                          = [ 1 : blocksize2 : ( blocksize2 * nblocks ); blocksize2 : blocksize2 : ( blocksize2 * nblocks ) ]';
blocks( nblocks, 2 )            = nelements2;

if nans( 1 ) > blocksize2 || nans( 2 ) > blocksize2
    error( 'internal (coding) bug' )
end

% protocode (dry run):
if dryrun
    
    d2new                       = zeros( nsamples1, 1 );
    d2os                        = 0;
    if nans( 1 )
        si                      = nans( 1 ) + 1;
    else
        si                      = 1;
    end
    ei                          = si;
    nsamples1hat                = ( nsamples1 - nans( 2 ) );
    for bnum                    = 1 : nblocks
        
        % collect the data
        sidx                    = ( ( blocks( bnum, 1 ) - 1 ) / nchans2 + 1 ) : ( blocks( bnum, 2 ) / nchans2 );
        d2b                     = d2( sidx );
        
        % actually warp
        while ei < nsamples1hat && t2w( ei ) < ( blocks( bnum, 2 ) / nchans2 )
            ei                  = ei + 1;
        end
        if verbose && ~mod( bnum, 50 )
            fprintf( 'dry run: block %d/%d: [ %d %d ] -> [ %d %d ]\n' ...
                , bnum, nblocks, ( blocks( bnum, 1 ) - 1 ) / nchans2 + 1, blocks( bnum, 2 ) / nchans2, si, ei )
        end
        kidx                    = t2w( si : ei ) - ( ( blocks( bnum, 1 ) - 1 ) / nchans2 );
        d2s                     = d2b( kidx );
        
        % increment
        si                      = ei + 1;

        % zero pad edges
        if bnum == 1 && nans( 1 )
            zpad                = zeros( nans( 1 ), 1 );
            d2s                 = [ zpad; d2s ];
        end
        if bnum == nblocks && nans( 2 )
            zpad                = zeros( nans( 2 ), 1 );
            d2s                 = [ d2s; zpad ];
        end
        
        % populate
        nb                      = length( d2s );
        didx                    = d2os + ( 1 : nb );
        d2new( didx )           = d2s;
        d2os                    = d2os + nb;
        
    end
    if ~isequal( d2new, d2w )
        error( 'internal (coding) bug' )
    end
    
end

%------------------------------------------------------------------------
% actually warp, offset, and write out file2 (blockwise)

% open files for reading and writing
fp0                             = fopen( newfile, 'w' );
fp2                             = fopen( sourcefile2, 'r' );
if fp0 == -1
    error( 'fopen error (%s)', newfile )
end
if fp2 == -1
    error( 'fopen error (%s)', sourcefile2 )
end

% go over the sourcefile in blocks and write the newfile
if nans( 1 )
    si                          = nans( 1 ) + 1;
else
    si                          = 1;
end
ei                              = si;
nsamples1hat                    = ( nsamples1 - nans( 2 ) );
for bnum                        = 1 : nblocks

    % collect the data
    if bnum == nblocks
        toload2                 = nelements2 - ( nblocks - 1 ) * blocksize2;
    else
        toload2                 = blocksize2;
    end
    data2                       = fread( fp2, toload2, precisionstr{ 2 } );
    
    % actually warp
    while ei < nsamples1hat && t2w( ei ) < ( blocks( bnum, 2 ) / nchans2 )
        ei = ei + 1;
    end
    if verbose && ~mod( bnum, 50 )
        fprintf( '\twet run: block %d/%d: [ %d %d ] -> [ %d %d ]\n' ... 
            , bnum, nblocks, ( blocks( bnum, 1 ) - 1 ) / nchans2 + 1, blocks( bnum, 2 ) / nchans2, si, ei )
    end
    kidx                        = t2w( si : ei ) - ( ( blocks( bnum, 1 ) - 1 ) / nchans2 );
    data2                       = reshape( data2, [ nchans2 toload2 / nchans2 ] );
    data2                       = data2( :, kidx );
    
    % increment
    si                          = ei + 1;
    
    % zero pad edges
    if bnum == 1 && nans( 1 )
        zpad                    = zeros( nchans2, nans( 1 ), precision{ 2 } );
        data2                   = [ zpad data2 ];
    end
    if bnum == nblocks && nans( 2 )
        zpad                    = zeros( nchans2, nans( 2 ), precision{ 2 } );
        data2                   = [ data2 zpad ];
    end
  
    % write out
    data2w                      = data2( : );
    fwrite( fp0, data2w, precision{ 2 } );
    
end

% close files
rc0                             = fclose( fp0 );
rc2                             = fclose( fp2 );
rc                              = [ rc0 rc2 ];

if verbose
	fprintf( 1, 'done warping and offsetting.\n' );
end

%------------------------------------------------------------------------
% check the warped sync signals from the warped file
if verbose
	fprintf( 1, '%s: verifying sync events post-warping...', mfname );
end

m2n                             = memmapfile( newfile, 'Format', precision{ 2 }, 'writable', false );
d2n                             = single( m2n.Data( chan2 : nchans2 : end ) );
d2n                             = rem( floor( ( d2n + shiftbits2 ) * pow2( 1 - bit2 ) ), 2 );

% get the transitions
dvec                            = diff( d2n );
hi                              = find( dvec == 1 );
lo                              = find( dvec == -1 );
nhi                             = length( hi ); 
nlo                             = length( lo );
n                               = min( nhi, nlo );
mat2n                           = [ hi( 1 : n ) lo( 1 : n ) ];

% check the files
devents                         = abs( mat1 - mat2n );
maxde                           = max( devents( : ) );
msg                             = sprintf( '%safter warping - max difference is %d samples\n' ...
    , msg, maxde );

if verbose
    fprintf( 1, 'maximal difference is %d samples\n', maxde )
end

return

% EOF

sourcefile1     = '/Volumes/slab1/mC400/raw/mC400_1_200825_173519/all_in_one.dat';
sourcefile2     = '/Volumes/slab1/mC400/raw/mC400_2_200825_173519/all_in_one.dat';
sync1           = [ 40 5 ];
sync2           = [ 9 5 ];
nchans1         = 40;
nchans2         = 9;
newfile         = '/Volumes/slab1/mC400/raw/mC400_2_200825_173519/all_in_one_warped.dat';
precision       = { 'int16', 'int16' };
shiftbits       = [ 0 2^15 ];
graphics        = 1;
verbose         = 1;
dryrun          = 1;

% warp sourcefile2 (longer, later syncs):
[ rc, msg, wf, os ] = warpfile( sourcefile1, sourcefile2, sync1, sync2, nchans1, nchans2, newfile, precision, shiftbits, graphics, verbose, dryrun );
