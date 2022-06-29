% removesamples        from dat file
%
% [ rc, msg ] = removesamples( fname, nchans, sampleIdx, forceRemove )
% 
% fname         full file name 
% nchans        number of channels
% sampleIdx     indices of samples to remove, in binary (not temporal) indices 
%               e.g. for 32 channels, 3 time points, the last temporal
%               sample is 3, whereas the last binary index is 96
% forceRemove   {0}. When 1, remove sampleIdx even if number of binary
%               samples is an integer multiple of nchans
%
% does
% rewrites the file without the samples to remove
%
% see also      removedcoffset, reorderchannels, removechannels

% 14-nov-17 ES

% revisions
% 03-dec-17 actually written
% 26-aug-20 cleaned up
% 29-sep-20 fixed s.t. does what it should

function [ rc, msg ]            = removesamples( fname, nchans, sampleIdx, forceRemove )

% initialize output
rc                              = 0;
msg                             = '';

% constants
nbytes                          = 2;                                        % [bytes/sample]
blocksize                       = 2^20;                                     % [samples]

% arguments
nargs                           = nargin;
if nargs < 3 || isempty( fname ) || isempty( nchans ) || isempty( sampleIdx )
    return
end
if ~issorted( sampleIdx )
    msg                         = sprintf( 'sampleIdx to remove must be sorted' );
    rc                          = -1;
    return
end
if nargs < 4 || isempty( forceRemove )
    forceRemove                 = 0;
end

if ~isa( fname, 'char' ) || ~exist( fname, 'file' ) 
    msg                         = sprintf( 'missing source %s', fname );
    rc                          = -1;
    return
end
tmpfname                        = [ fname '.tmp' ];
if ~isa( nchans, 'numeric' )
    msg                         = 'improper format for nchans';
    rc                          = -1;
    return
end

% partition into blocks
info                            = dir( fname );
nsamples                        = info.bytes / nbytes;
if isequal( nsamples / nchans, round( nsamples / nchans ) )
    msg                         = sprintf( 'note that number of samples is an integer multiple of nchans (%d) for file %s', nchans, fname );
    if ~forceRemove
        rc                      = -1;
        return
    end
end
nblocks                         = ceil( nsamples / blocksize );
blocks                          = [ 1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks ]';
blocks( nblocks, 2 )            = nsamples;

% open file for writing
fid                             = fopen( tmpfname, 'w' );
if fid == -1
    msg                         = 'cannot open file';
    rc                          = fid;
    return
end

% go over blocks and write out
j                               = 1;
for i                           = 1 : nblocks
    boff                        = ( blocks( i, 1 ) - 1 ) * nbytes;
    bsize                       = ( diff( blocks( i, : ) ) + 1 );
    m                           = memmapfile( fname, 'Format', 'int16', 'Offset', boff, 'Repeat', bsize, 'writable', false );
    if sampleIdx( j ) < blocks( i, 2 )
        disp( 'here!!' )
        idx1                    = blocks( i, 1 ) : blocks( i, 2 );
        idx2                    = sampleIdx( j ) : sampleIdx( end );
        %kidx                    = setdiff( idx1, idx2 );
        [ ~, i1, i2 ]           = intersect( idx1, idx2 );
        kidx                    = setdiff( 1 : length( idx1 ), i1 );
        d                       = m.data( kidx );
        j                       = j + length( i2 );
    else
        d                       = m.data;
    end
    fwrite( fid, d( : ), 'int16' );
    clear d m
end

% close file
rc                              = fclose( fid );
if rc == -1
    msg                         = 'cannot save new file';
    return
end

% remove original file and rename new one
if isunix
    cmd                         = sprintf( '!rm %s; mv %s %s', fname, tmpfname, fname );
elseif ispc
    cmd                         = sprintf( '!del %s; rename %s %s', fname, tmpfname, fname );
else
    msg                         = 'unrecognized system';
    rc                          = -1;
    return
end
eval( cmd )

return

% EOF

