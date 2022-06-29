% Inputs:
% filebase        full path and session name e.g '/odin/Recordings/mF105/dat/mF105_10/mF105_10';
% nchans          number of recordings saits in the shank e.g 11
% clu
% shankNum        shanknumber as a string e.g '4'
% ignoredClu



function [mspk,sspk,nspk_vec] = mWaveForm(filebase, nchans, clu, shanknum, ignoredClu)
nsamples                    = 32;

if ~exist ('ignoredClu')
    ignoredClu     = [];
end

% default values
BLOCKSIZE                   = 10 * 2^20;            % bytes/block; use 10 MB blocks
precision0                  = 'int16';              % parameters for block-wise loading of data
precision1                  = 'single';
%ignoredClu                  = [ 0 1 ]; % do not re-align noise and artifacts

% set up the block-wise loading of spiking - check if needed
precision                   = sprintf( '%s=>%s', precision0, precision1 );                  % build the type casting string
a0                          = ones( 1, 1, precision0 );                                     % determine number of bytes/sample
a1                          = ones( 1, 1, precision1 );
sourceinfo                  = whos( 'a0' );
targetinfo                  = whos( 'a1' );
nbytes0                     = sourceinfo.bytes;
nbytes1                     = targetinfo.bytes;
spkfname                    = sprintf( '%s.spk.%s', filebase, shanknum );

% determine number of spikes in the file
spksize                     = nchans * nsamples * nbytes0;
s                           = dir( spkfname );
spksource                   = s.bytes / spksize;
kidx                        = true( spksource, 1 );
fkidx                       = find( kidx );    

uclu                        = unique(clu);     
%shankclus                   =  uclu';
  
shankclus                   = setdiff( uclu( : )', ignoredClu( : )' );
nclusShank                  = length( shankclus );       

%------------------------------------------------------------
% pass 1/3: load spk waveforms (blockwise), compute clu-specific statistics
% divide into blocks
spkblk                      = floor( BLOCKSIZE / ( nchans * nsamples * nbytes1 ) );                     % number of spikes/block in RAM
nblocks                     = ceil( spksource / spkblk );
blocks                      = [ 1 : spkblk : spkblk * nblocks; spkblk : spkblk : spkblk * nblocks ]';
blocks( nblocks, 2 )        = spksource;
%fprintf( 1, '; %d blocks; Block number ', nblocks )

 % initialize variables for accumulation
 nspks_CLU                   = zeros( 1, nclusShank );
 mspks_SUM                   = single( zeros( nchans, nsamples, nclusShank ) );
 mspks_SSQ                   = single( zeros( nchans, nsamples, nclusShank ) );
 
 % open file for reading
 fp                          = fopen( spkfname, 'r' );
 if fp == -1
    error( 'fopen error' )
 end
 % go over blocks
 for bidx                    = 1 : nblocks
    if ~mod( bidx, 10 )
%        fprintf( 1, '%d ', bidx );
    end
    nspkBlock               = blocks( bidx, 2 ) - blocks( bidx, 1 ) + 1;
    startpos                = ( blocks( bidx, 1 ) - 1 ) * spksize;
    datasize                = [ nchans nsamples * nspkBlock ];
    rc                      = fseek( fp, startpos, 'bof' );
    if rc
        error( 'fseek error' )
    end
    spks                    = fread( fp, datasize, precision );
    spks                    = reshape( spks, [ nchans nsamples nspkBlock ] );                           % organize output
    spks                    = spks( :, :, kidx( blocks( bidx, 1 ) : blocks( bidx, 2 ) ) );              % select subset of spikes and their corresponding labels
    cluBlock                = clu( find( fkidx >= blocks( bidx, 1 ), 1 ) : find( fkidx <= blocks( bidx, 2 ), 1, 'last' ) );
    % go over clusters and accumulate statistics
    ci                      = 0;
    for clunum              = shankclus
        ci                  = ci + 1;
        idx                 = cluBlock == clunum;
        nspks_CLU( ci )     = nspks_CLU( ci ) + sum( idx );
        if ~sum( idx )
            continue
        end
        spk                 = spks( :, :, idx );
        mspks_SUM( :, :, ci )   = mspks_SUM( :, :, ci ) + sum( spk, 3 );              % spike amplitude (average)
        mspks_SSQ( :, :, ci )   = mspks_SSQ( :, :, ci ) + ( sum( spk .^ 2, 3 ) );
    end
 end
 
fclose( fp );
%fprintf( 1, '\n' )

%------------------------------------------------------------
% pass 2/3: compute statistics
verbose                     = 1;
mspk                        = zeros( nchans, nsamples, nclusShank );
sspk                        = zeros( nchans, nsamples, nclusShank );
nspk_vec                    = zeros(1, nclusShank);
for clunum                  = shankclus
    ci                      = find( shankclus == clunum );
    idx                     = clu == clunum;
    nspks                   = sum( idx );
    nspk_vec(ci)            = nspks; 
    if ~nspks
        if verbose
            fprintf( 1, 'no spikes for %s.%d (%d)\n', shanknum, clunum, nspks );
        end
        continue
    end
    if verbose
%        fprintf( 1, '\t%s.%d (%d)', shanknum, clunum, nspks );
    end
       
    % summarize/compute waveform statistics for that cluster
    % divide the accumulated statistics by the number of spikes in that cluster
    mspk( :, :, ci )        = mspks_SUM( :, :, ci ) / nspks_CLU( ci );                                  % spike amplitude (average, then p2p)
    sspk( :, :, ci )        = sqrt( mspks_SSQ( :, :, ci ) / nspks_CLU( ci ) - mspk( :, :, ci ) .^ 2 );  % SD = sqrt( <x^2> - <x>^2 );
       
end
%fprintf( 1, '\n' );
end
% spk = readspk([filebase1,'.spk.3'],11,32);
% clu1 =clu;
% clu1(1)=[];
% nclu =unique(clu);
% idx =clu1==nclu(10);
%mspk10 =mean(spk(:,:,idx),3);