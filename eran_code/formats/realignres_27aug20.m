% realignres                realign res file according to cluster-specific extemum
%
% call                      rc = realignres( filebase )
% 
% gets                      filebase 
% 
% optional arguments        halfLength          {8}; half length of detection
%                           method              {'extremum'} or 'template
%                           byPar               {0}; see spikes_stats
%
% requirements: 
%                           *xml for the filebase
%                           *clu.#, *res.#, and *spk.# for the nElecGps in the *xml
%
% does: 
%                           (1) loads each clu/res
%                           (2) goes over clusters and computes statistics (channel, extremum, polarity, time, mean, SD)
%                           (3) goes over all invidiual spikes in each cluster and realigns the res
%                           (4) write res to disk
%
% calls                     ParseArgPairs, LoadXml, resort
%
% see also                  uniqueres, dat2spkfet
%
% motivating problems
% 
% (1) many clusters consist of variable waveforms, some of the variability is
%       due to suboptimal alignment since the initial alignment is agnostic of
%       the cluster identity (and is done according to the extremum only)
% 
% (2) some clusters are merged from two clusters, and then not all spikes
%       have an extermum at the same polarity/time (can happen following manual
%       curation; can be of "mirrored" PCs or of "bipolar" units)
% 
% (3) in some clusters, the same waveform is detected twice, at different
%       time lags (can happen with template-matching engines)
%
% (4) in some clusters, the best match (and corresponding time) may be at an 
%       offset from the prescribed peak (e.g. 16-th) sample (can happen with template-matching engines)

% 16-aug-20 LS+ES

% revisions 
% 18-aug-20 wrote blockwise loading of spk (first pass)
% 23-aug-20 wrote statistcs (second pass)
%           started to write offset computations (third pass)
% 25-aug-20 finalized coding of algorithm 2
% 27-aug-20 cleaned up, documented

% to do:
% (1) write algorithm 1 (template matching)

function rc = realignres( filebase, varargin )

%------------------------------------------------------------
% constants
%------------------------------------------------------------
% constants
verbose                         = 1;

% default values
HALFLENGTH                      = 8;            % [samples]

BLOCKSIZE                       = 10 * 2^20;            % bytes/block; use 10 MB blocks
precision0                      = 'int16';              % parameters for block-wise loading of data
precision1                      = 'single';

ignoredClu                      = [ 0 1 ]; % do not re-align noise and artifacts

%------------------------------------------------------------
% arguments
%------------------------------------------------------------
nargs = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ halfLength ...
    , method ...
    , byPar ] = ParseArgPairs(...
    { 'halfLength' ...
    , 'method' ...
    , 'byPar' } ...
    , {  HALFLENGTH ...
    , 'extremum' ...
    ,  0 } ...
    , varargin{ : } );
if ~ismember( method, { 'template', 'extremum' } )
    error( 'unsupported method' )
end

%------------------------------------------------------------
% collect information
%------------------------------------------------------------
% *.m file name
mfname                          = upper( mfilename );

% xml file
fprintf( 1, '%s: Loading *xml file...', mfname, filebase )
par                             = LoadXml( sprintf( '%s.xml', filebase ) );
fprintf( 1, 'done.\n' )

nshanks                         = par.nElecGps;
nchans_all                      = zeros( 1, nshanks );
for si                          = 1 : nshanks
    nchans_all( si )            = length( par.SpkGrps( si ).Channels );
end
nsamples                        = par.SpkGrps( 1 ).nSamples;                                    % assume same for all spike groups
peaksample                      = par.SpkGrps( 1 ).PeakSample;                                  % assume same for all spike groups
p2pvoltage                      = par.VoltageRange;
scalefactor                     = 1 / 2^par.nBits * p2pvoltage * 1e6 / par.Amplification;       % a2d units -> microVolts

% set up the block-wise loading of spiking
precision                       = sprintf( '%s=>%s', precision0, precision1 );                  % build the type casting string
a0                              = ones( 1, 1, precision0 );                                     % determine number of bytes/sample
a1                              = ones( 1, 1, precision1 );
sourceinfo                      = whos( 'a0' );
targetinfo                      = whos( 'a1' );
nbytes0                         = sourceinfo.bytes;
nbytes1                         = targetinfo.bytes;

% extremum detection/template matching parameters
tmpidx                          = peaksample + ( -halfLength : halfLength );                    % to prevent spurious edge effects

%------------------------------------------------------------
% go over shanks
%------------------------------------------------------------
rc                              = NaN * ones( nshanks, 2 );
ushanks                         = 1 : nshanks;
for sidx                        = 1 : nshanks
    
    % determine the shank-specific spikes
    shanknum                    = ushanks( sidx );
    indfname                    = sprintf( '%s.ind.%d', filebase, shanknum );
    clufname                    = sprintf( '%s.clu.%d', filebase, shanknum );
    resfname                    = sprintf( '%s.res.%d', filebase, shanknum );
    spkfname                    = sprintf( '%s.spk.%d', filebase, shanknum );
    chans                       = sort( par.SpkGrps( sidx ).Channels ) + 1;
    nchans                      = length( chans );
    if byPar
        ridx                    = resort( par.SpkGrps( shanknum ).Channels  );                  % reorder according to par file order
    else
        ridx                    = 1 : length( par.SpkGrps( shanknum ).Channels );
    end
    
    if exist( indfname, 'file' )
        error( 'not supported yet' )
    end
    
    %------------------------------------------------------------
    % load res file
    fprintf( 1, '%s: Loading shank#%d res... ', mfname, shanknum )
    res                         = load( resfname );
    nspks                       = length( res );
    fprintf( 1, '%d spikes. ', nspks )
    res0                        = res;

    %------------------------------------------------------------
    % load clu file
    fprintf( 1, 'clu... ' )
    clu                         = load( clufname );
    nclu0                       = clu( 1 );
    clu( 1 )                    = [];
    uclu                        = unique( clu );
    nclu                        = length( uclu );
    if nclu0 ~= nclu
        warning( 'weird clu file' )
    end
    fprintf( 1, '%d clu. ', nclu )
    shankclus                   = setdiff( uclu( : )', ignoredClu( : )' );
    nclusShank                  = length( shankclus );                      % number of relevant units
    
    %------------------------------------------------------------
    % determine number of spikes in the file
    spksize                     = nchans * nsamples * nbytes0;
    s                           = dir( spkfname );
    spksource                   = s.bytes / spksize;
    kidx                        = true( spksource, 1 );
    fkidx                       = find( kidx );
    
    %------------------------------------------------------------
    % pass 1/3: load spk waveforms (blockwise), compute clu-specific statistics
    % divide into blocks
    spkblk                      = floor( BLOCKSIZE / ( nchans * nsamples * nbytes1 ) );                     % number of spikes/block in RAM
    nblocks                     = ceil( spksource / spkblk );
    blocks                      = [ 1 : spkblk : spkblk * nblocks; spkblk : spkblk : spkblk * nblocks ]';
    blocks( nblocks, 2 )        = spksource;
    fprintf( 1, '; %d blocks; Block number ', nblocks )
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
            fprintf( 1, '%d ', bidx )
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
        spks                    = spks( ridx, :, : );                                                       % reorder spikes according to par file
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
            mspks_SUM( :, :, ci )   = mspks_SUM( :, :, ci ) + sum( spk, 3 ) * scalefactor;              % spike amplitude (average)
            mspks_SSQ( :, :, ci )   = mspks_SSQ( :, :, ci ) + ( sum( spk .^ 2, 3 ) ) * scalefactor^2;
        end
    end
    fclose( fp );
    fprintf( 1, '\n' )
        
    %------------------------------------------------------------
    % pass 2/3: compute statistics
    mspk                        = zeros( nchans, nsamples, nclusShank );
    sspk                        = zeros( nchans, nsamples, nclusShank );
    extchan                     = zeros( nclusShank, 1 );
    extval                      = zeros( nclusShank, 1 );
    extsamp                     = zeros( nclusShank, 1 );
    for clunum                  = shankclus
        
        ci                      = find( shankclus == clunum );
        idx                     = clu == clunum;
        nspks                   = sum( idx );       
        if ~nspks
            if verbose
                fprintf( 1, 'no spikes for %d.%d (%d)\n', shanknum, clunum, nspks )
            end
            continue
        end
        if verbose
            fprintf( 1, '\t%d.%d (%d)', shanknum, clunum, nspks )
        end
        
        % summarize/compute waveform statistics for that cluster
        % divide the accumulated statistics by the number of spikes in that cluster
        mspk( :, :, ci )        = mspks_SUM( :, :, ci ) / nspks_CLU( ci );                                  % spike amplitude (average, then p2p)
        sspk( :, :, ci )        = sqrt( mspks_SSQ( :, :, ci ) / nspks_CLU( ci ) - mspk( :, :, ci ) .^ 2 );  % SD = sqrt( <x^2> - <x>^2 );
        % determine the extermum value, channel, and polarity
        maxspk                  = max( mspk( :, :, ci ), [], 2 );
        minspk                  = min( mspk( :, :, ci ), [], 2 );
        maxval                  = max( [ max( abs( minspk ) ) max( maxspk ) ] );
        minchan                 = find( abs( minspk ) == maxval );
        maxchan                 = find( maxspk == maxval );
        if isempty( minchan )                                                                               % Punit
            extchan( ci )       = maxchan;
            extval( ci )        = maxval;
        else                                                                                                % Nunit
            extchan( ci )       = minchan;
            extval( ci )        = -maxval;
        end
        extsamp( ci )           = find( mspk( extchan( ci ), :, ci ) == extval( ci ) );
        
        % figure, plot_spk_waveforms( cat( 3, mspk, sspk ), [], [], 1, scalefactor, [ 1 1 0 ] );
        
    end
    fprintf( 1, '\n' )
    
    %------------------------------------------------------------
    % pass 3/3: go over clusters and realign each
    
    % initialize variables for accumulation
    dres                        = NaN( size( res ) );
    % open file for reading
    fp                          = fopen( spkfname, 'r' );
    if fp == -1
        error( 'fopen error' )
    end
    % go over blocks
    for bidx                    = 1 : nblocks
        if ~mod( bidx, 10 )
            fprintf( 1, '%d ', bidx )
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
        spks                    = spks( ridx, :, : );                                                       % reorder spikes according to par file
        bcidx                   = find( fkidx >= blocks( bidx, 1 ), 1 ) : find( fkidx <= blocks( bidx, 2 ), 1, 'last' );
        cluBlock                = clu( bcidx );
        
        % go over clusters and accumulate statistics
        for clunum              = shankclus
            
            ci                  = find( shankclus == clunum );
            idx                 = cluBlock == clunum;
            spk                 = permute( spks( extchan( ci ), :, idx ), [ 2 3 1 ] ) * scalefactor;        % take only the relevant waveforms - clu, chan
            
            %figure, plot( spk, 'b' ), hold on, ph = plot( tmpl, 'k' ); set( ph, 'LineWidth', 2 ), hold on, plot( extsamp( ci ), extval( ci ), 'or' )            
            
            switch method
                case 'template'
                    
                    % method1 (ideal for noisy waveforms, double detection):
                    %   take the 16 central samples of the mean as a template
                    %   find the optimal offset for each specific spike (-7 to 8) using
                    %       template matching
                    tmpl        = permute( mspk( extchan( ci ), :, ci ), [ 2 3 1 ] );                  % compute template
                    
                case 'extremum'
                    
                    % method2 (ideal for bipolar units):
                    %   use the extremum of the mean (polarity, channel)
                    %   find the signed extremum of each specific spike
                    %   compute the time difference between the time of the extremum and the "peaksample"
                    if extval( ci ) < 0
                        [ ~, extidx ] = min( spk( tmpidx, : ), [], 1 );      % look for minimum
                    else
                        [ ~, extidx ] = max( spk, [], 1 );      % look for maximum
                    end
                    extidx              = extidx + tmpidx( 1 ) - 1;
                    dres( bcidx( idx ) ) = extidx - peaksample;
            end     % switch
            
        end         % clunum
        
    end             % blocks
    fclose( fp );
    
    fprintf( 1, '\n' )
    
    %------------------------------------------------------------
    % apply offset to res
    if any( ~ismember( unique( clu( isnan( dres ) ) ), ignoredClu ) )
        error( 'problem' )
    end
    dres( isnan( dres ) )       = 0;
    res                         = res0 + dres;

    %------------------------------------------------------------
    % save res file
    % structure: ASCII file, one row per spike time, in samples
    res                         = int32( res );           % this supports 29 hours of recordings
    fid                         = fopen( resfname, 'w' );
    [ ~ ]                       = fprintf( fid, '%d\n', res );
    rc( sidx, 1 )               = fclose( fid );
    if rc( sidx, 1 ) == 0
        fprintf( 1, '%d spikes.\n', length( res ) )
    else
        fprintf( 1, 'failed writing file: %s!!\n', resfname  )
    end
 
end

return

% EOF
