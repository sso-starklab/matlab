% makeblocks            for batch processing
%
% blocks = makeblocks( n, blocksize, noverlap )
%
% n             data size (or a list of periods)
% blocksize     blocksize
% noverlap      {0}
%
% blocks        nblocks x 2 matrix

% 26-dec-12 ES

% revisions
% 29-feb-16 revised noverlap case
% 16-may-16 bug correction
% 20-may-20 cleaned up

function blocks         = makeblocks( n, blocksize, noverlap )

nargs                   = nargin;
if nargs < 2 || isempty( blocksize )
    blocksize           = n;
end
if nargs < 3 || isempty( noverlap )
    noverlap            = 0;
end
if noverlap < 0
    noverlap            = 0;
end
if noverlap >= blocksize
    noverlap            = blocksize - 1;
end
noverlap                = round( noverlap );

if noverlap
    start               = ( 1 : ( blocksize - noverlap ) : n )';
    blocks              = [ start start + blocksize - 1 ];
    nblocks             = size( blocks, 1 );
else
    nblocks             = ceil( n / ( blocksize - noverlap ) );
    blocks              = [ 1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks ]';
end

blocks( nblocks, 2 )    = n;
blocks( blocks > n )    = n;

return

% EOF
