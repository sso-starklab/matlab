% spk2fet               loads *.spk.# and *.res.# data, computes PC basis, projects on the basis, and saves a *.fet.# file
%
% call                  [ rc, fet, basis, varex ] = spk2fet( spkfname, resfname, fetfname, nchans )
%
% gets                  spkfname            full path and name of spk file (input)
%                       resfname            full path and name of spk file (input)
%                       fetfname            full path and name of fet file (output)
%                       nchans           number of channels in spk file
%
% optional arguments (name/value pairs)
% 
%                       nsamples            {32} samples per spike per channel
%                       nupsample           {4} upsampling factor
%                       npcs                {3} number of eigenvectors to keep
%
% returns               rc                  flag (0 is OK)
%                       fet                 nfeatures x nspikes matrix of PCA coefficients (int32)
%                                               npcs per channel (for all channels), energy for each channel, then res
%                       basis               eigenvectors per channel; npcs x ( nsamples x nupsample ) x nchans
%                       varex               variance explained by each e.v. per channel; npcs x nchans
%
% does                  loads res
%                       first pass - goes over blocks:
%                           loads spikes
%                           removes the mean of every spike
%                           fft-upsamples each spike
%                           computes covariance matrix of each channel
%                           accumulates covariance matrix over blocks
%                       second pass - goes over channels:
%                           computes the eigenvectors and keeps the first PCs
%                       third pass - goes over blocks:
%                           loads spikes
%                           removes the mean of every spike
%                           fft-upsamples each spike
%                           projects every spike on these PCs
%                           computes the energy (norm) of each (upsampled) spike
%                           builds a fet file block with PCA coefficients and energy of every spike, and spike times
%                           writes the block to the fet file
%
% note                  if there are multiple simultaneous spikes, their
%                           PCA representation will be identical
%
% calls                 ParseArgPairs, fft_upsample

% 19-may-20 LS + ES 

% revisions
% 20-may-20 (1) moved call to makeblocks out of loop
%           (2) added fetBits to list of optional arguments
% 07-mar-21 (1) modified to block structure (not memory intensive any more)

function [ rc, fet, basis, varex ] = spk2fet( spkfname, resfname, fetfname, nchans, varargin )

%------------------------------------------------------------------------
% preps
%------------------------------------------------------------------------

% constants
target                          = 'int32';                                  % temporal (res) resolution in fet file
precision0                      = 'int16';                                  % parameters for block-wise loading of data
precision1                      = 'single';
BLOCKSIZE                       = 10 * 2^20;                                % bytes/block; use 10 MB blocks

% default values
NSAMPLES                        = 32;                                       % [samples]
NUPSAMPLE                       = 4;                                        % natural number
NPCS                            = 3;                                        % number of eigenvectors to keep
FETBITS                         = 16;                                       % number of bits in each non-temporal feature

% arguments
nargs                           = nargin;
if nargs < 4 || isempty( spkfname ) || isempty( resfname ) || isempty( fetfname ) || isempty( nchans )
    return
end
[ nsamples, nupsample ...
    , npcs, fetBits ...
    , verbose ]                 = ParseArgPairs(...
    { 'nsamples', 'nupsample' ...
    , 'npcs', 'fetBits' ...
    , 'verbose' } ...
    , { NSAMPLES, NUPSAMPLE ...
    , NPCS, FETBITS ...
    , 1 } ...
    , varargin{ : } );
M                               = nsamples * nupsample;

% initialize output
rc                              = 0;
fet                             = [];
basis                           = NaN( npcs, M, nchans );
varex                           = NaN( npcs, nchans );

% check if files exist
if ~exist( spkfname, 'file' )
    fprintf( 1, '%s: missing file %s, exiting\n', upper( mfilename ), spkfname )
    return
end
if ~exist( resfname, 'file' )
    fprintf( 1, '%s: missing file %s, exiting\n', upper( mfilename ), resfname )
    return
end

% set up the block-wise loading of spiking
precision                       = sprintf( '%s=>%s', precision0, precision1 );                  % build the type casting string
a0                              = ones( 1, 1, precision0 );                                     % determine number of bytes/sample
a1                              = ones( 1, 1, precision1 );
sourceinfo                      = whos( 'a0' );
targetinfo                      = whos( 'a1' );
nbytes0                         = sourceinfo.bytes;
nbytes1                         = targetinfo.bytes;

% load res file
if verbose
    fprintf( 1, '%s: loading res file %s\n', upper( mfilename ), resfname )
end
fid                             = fopen( resfname, 'r' );
res                             = fscanf( fid, '%d' );
fclose( fid );

%------------------------------------------------------------
% pass 1/3: load spk waveforms (blockwise), compute covariance matrix for each channel
%------------------------------------------------------------

% determine number of spikes in the spk file
spksize                         = nchans * nsamples * nbytes0;
s                               = dir( spkfname );
spksource                       = s.bytes / spksize;
if spksource ~= length( res )
    error( 'input size mismatch' )
end
covxs                           = zeros( M, M, nchans );

% divide into blocks
spkblk                          = floor( BLOCKSIZE / ( nchans * nsamples * nbytes1 ) );                     % number of spikes/block in RAM
nblocks                         = ceil( spksource / spkblk );
blocks                          = [ 1 : spkblk : spkblk * nblocks; spkblk : spkblk : spkblk * nblocks ]';
blocks( nblocks, 2 )            = spksource;

% open file for reading
fp                              = fopen( spkfname, 'r' );
if fp == -1
    error( 'fopen error' )
end

% go over blocks
if verbose
    fprintf( 1, '\tpass 1/3: covariance matrices, %d blocks ', nblocks )
end

for bidx                        = 1 : nblocks
    if verbose
        fprintf( 1, '.' )
        if ~mod( bidx, 10 )
            fprintf( 1, '%d ', bidx )
        end
    end
    nspkBlock                   = blocks( bidx, 2 ) - blocks( bidx, 1 ) + 1;
    startpos                    = ( blocks( bidx, 1 ) - 1 ) * spksize;
    datasize                    = [ nchans nsamples * nspkBlock ];
    rc                          = fseek( fp, startpos, 'bof' );
    if rc
        error( 'fseek error' )
    end
    spks                        = fread( fp, datasize, precision );
    spks                        = reshape( spks, [ nchans nsamples nspkBlock ] );                           % organize output
    
    % go over channels and accumulate statistics
    for cnum                    = 1 : nchans
        
        x                       = permute( spks( cnum, :, : ), [ 2 3 1 ] ); % x: extract the data from one channel, spikes in columns     
        xmean                   = mean( x );
        xhat                    = x - ones( nsamples, 1 ) * xmean;          % xhat: remove mean
        y                       = fft_upsample( xhat, nupsample, 1 );       % y: upsample using fft
        covx                    = y * y';                                   % covx: compute covariance matrix
        covxs( :, :, cnum )     = covxs( :, :, cnum ) + covx;               % accumulate covx
        
    end     % channels
end         % blocks
fclose( fp );
covxs                           = covxs / ( spksource - 1 );
fprintf( 1, 'done covariances.\n' )

%------------------------------------------------------------
% pass 2/3: compute statistics
%------------------------------------------------------------

% loop over channels
if verbose
    fprintf( 1, '\tpass 2/3: eigenvectors, %d channels ', nchans )
end

for cnum                    = 1 : nchans
    
    if verbose
        fprintf( 1, '%d ', cnum )
    end

    % A: compute eigenvectors from the covariance matrix
    covx                    = covxs( :, :, cnum );
    [ U, S ]                = eig( covx );
    [ s, ind ]              = sort( diag( S ) );                % sort
    s                       = flipud( s );
    ind                     = flipud( ind );
    A                       = U( :, ind( 1 : npcs ) );          % take first npcs
    shat                    = s / sum( s );
    evals                   = shat( 1 : npcs );                 % keep for information

    % accumulate statistics
    basis( :, :, cnum ) 	= A';
    varex( :, cnum )        = evals;
    
end

if verbose
    fprintf( 1, 'done eigenvectors.\n' )
end

%------------------------------------------------------------
% pass 3/3: load spk waveforms (blockwise), project, and write to fet file
%------------------------------------------------------------

% determine parameters for the fet file
ncols                           = ( npcs + 1 ) * nchans + 1;

% open file for reading
fp                              = fopen( spkfname, 'r' );
if fp == -1
    error( 'fopen error' )
end

% open file for writing
formatstr                       = sprintf( '%%d%s\\n', repmat( '\t%d', [ 1 ncols - 1 ] ) );
fid                             = fopen( fetfname, 'w' );
if fid == -1
    error( 'fopen error' )
end

% go over blocks
if verbose
    fprintf( 1, '\tpass 3/3: projecting features and writing fet file, %d blocks ', nblocks )
end

for bidx                        = 1 : nblocks
    if verbose
        fprintf( 1, '.' )
        if ~mod( bidx, 10 )
            fprintf( 1, '%d ', bidx )
        end
    end
    nspkBlock                   = blocks( bidx, 2 ) - blocks( bidx, 1 ) + 1;
    startpos                    = ( blocks( bidx, 1 ) - 1 ) * spksize;
    datasize                    = [ nchans nsamples * nspkBlock ];
    rc                          = fseek( fp, startpos, 'bof' );
    if rc
        error( 'fseek error' )
    end
    spks                        = fread( fp, datasize, precision );
    spks                        = reshape( spks, [ nchans nsamples nspkBlock ] );                           % organize output
    
    % go over channels and accumulate statistics
    fet                         = zeros( ncols, nspkBlock, 'single' );
    for cnum                    = 1 : nchans

        % compute features
        x                       = permute( spks( cnum, :, : ), [ 2 3 1 ] ); % x: extract the data from one channel, spikes in columns     
        xmean                   = mean( x );
        xhat                    = x - ones( nsamples, 1 ) * xmean;          % xhat: remove mean
        y                       = fft_upsample( xhat, nupsample, 1 );       % y: upsample using fft
        A                       = basis( :, :, cnum )';
        coef                    = A' * y;                                   % coef: project X on A
        e                       = sqrt( sum( y .^ 2 ) / nupsample );        % energy: compute norm of each spike
        
        % accumulate over channels
        ridx1                   = cnum * npcs + ( -npcs + 1 : 0 );
        ridx2                   = npcs * nchans + cnum;
        fet( ridx1, : )         = coef;
        fet( ridx2, : )         = e;
        
    end     % channels
    
    % scale to fetBits and type cast to int32
    fet                         = int32( fix( fet / max( abs( fet( : ) ) ) * 2 ^ ( fetBits - 1 ) ) );
    
    % add row of res
    kidx                        = blocks( bidx, 1 ) : blocks( bidx, 2 );
    eval( sprintf( 'resb = %s( res( kidx ) );', target ) );
    ridx                        = ncols;
    fet( ridx, : )              = resb;
    
    % write to fet file
    if bidx == 1
       	fprintf( fid, '%d\n', ncols );
    end
    fprintf( fid, formatstr, fet );
   
end         % blocks
fclose( fp );
rc                          = fclose( fid );

if verbose
    fprintf( 1, 'done fet file.\n' )
end

if rc == 0
    if verbose
        fprintf( 1, '%s: written file %s\n', upper( mfilename ), fetfname )
    end
else
    fprintf( 1, '%s: failed writing file %s\n', upper( mfilename ), fetfname )
end

return

% EOF


% rc( sidx, 2 )           = spk2fet_blocked( spkfname, resfname, fetfname, ngchans ...
%             , 'fetBits', fetBits, 'verbose', 1, 'npcs', nrank );