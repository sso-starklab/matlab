% spikes_stats        waveforms, ACH, and spectra
%
% CALL              [ out ] = spikes_stats( filebase )
%
% GETS              filebase          full path
%
% OPTIONAL ARGUMENTS (given as name/value pairs):
%
% shanknums         any subset; defaults to all
% clunums           any subset; LoadCluRes format
% clustr            {'clu'}
% plotflag          {0} to plot (separate figure for each clu and summary for features)
% Overwrite         {0} to load w/o computing if previously computed; 1 to recompute and overwrite; -1 to just recompute 
% periods           {-1}. periods (in samples) to keep spikes from. overloaded: 
%                   -1 indicates to remove spikes during segments specified in val files; 
%                   -XX indicates to remove only during segments in val.tXX (can be a vector for multiple triggers).
%                       right now the combination of the two is not suppported (i.e. keep spikes in certain periods but not during
%                       stimulation etc). also not relevant for merged files 
% blockwise         {1} to work in blocks
% byPar             {0} reorder according to par file order
% doFET             {1} 
%
% RETURNS           sst
% 
% DOES              loads spikes and times, computes center of mass,
%                   spectrum and waveform features of peak waveform, and
%                   various timing statistics
%
% CALLS             ParseArgPairs, LoadXml, LoadCluRes, LoadVals
%                   dilutesegments, setdiffranges
%                   resort, readspk, LoadClu, LoadRes, LoadFet, inranges
%                   my_spectrum, calc_spk_snrs, calc_com, CCG
%                   plot_spk_waveforms, wfeatures, fft_upsample
%
% see also          compute_cluster_quality, plot_ss, plot_spikes_stats

% 08-jan-12 ES

% revisions
% 16-jan-12 cluster quality measures added
% 17-jan-12 cell type measures added
% 24-jan-12 modifications in Overwrite
% 25-jan-12 waveform classification into PYR/IN added
% 30-jan-12 periods argument added (supports periods/trigger time removal)
%           freqs field filled (to support plot_ss)
% 28-feb-12 saves if overwrite or no sst file; skip if no spikes (assign
%           zero to all
% 08-jun-12 (1) changed algorithm to block-wise permitting structure for
%           handling huge *.spk.* files
%           (2) changed format to support r2009b
%           (3) handled case of no spikes at all
% 10-jun-12 (1) FET handling improved - single; blockwise mahal
% 06-sep-13 (1) p2pvoltage from xml file
% 29-sep-13 (1) reorder channels according to par channel (optional)
% 22-jul-15 (1) corrrected frate computation - now by durValid
% 20-nov-15 (1) fixup for durValid (use *eeg file if missing *dat)
% 23-feb-18 (1) clean up after EOF
%           (2) change classify_waveform to GMM
% 23-feb-19 (1) fixed temporary hack for HPF spike waveforms (line 614)
% 01-aug-19 (1) changed DefaultArgs to ParseArgPairs
% 17-aug-19 cleaned up
% 09-sep-19 (1) added support of *.ind.* files (reorganization of fet; offset of geo_com)
%           (2) added new descriptor - geo_sd - spatial dispersion (units of inter-site intervals)
% 09-jul-20 (1) added new descriptor - geo_fwhm - spatial dispersion (units of inter-site intervals)
%           (2) added new descriptor - extremum - signed amplitude at
%                   extremum of sum-of-squares of each waveform
% 25-oct-20 (1) added DTu (defaulted to DT, 8 samples) to remove all spikes
%               that occurred with an ISI <= DTu from all quality analyses
%               this is to ensure that detections based on template
%               matching are not over-penalized d.t. simultaneous detection
%               of two units on the shank. 
%
%           (2) note that for template matching algorithms with solution of
%               the overlap problem there can be spikes of the same unit
%               during the DT and thus the corrected ISIindex may be
%               over-permissive. In such cases, DT shoud be reduced to the
%               actual dead-time (e.g. 0). To account for this, both DTu
%               and DT are arguments and not constants. 
%
%           (3) added geo_comh and geo_sdh for localization and dispersion
%               in the horizontal axis (mainly relevant for grid i.e. DBC
%               probes)

% problems:
% (1) bursting cells clustered as 2 distinct clusters - can cause the first
% cluster to seem inh+exc (e.g. es25nov11_3 1.34->1.27)
% (2) FS IN: sometime form a synchronized network, then there are peaks in
% the CCHs (e.g. es26nov11_3: 2.12;3.17->1.11;3.6=42;65->10,54) that not
% always stradle zero-lag, so the cells may seem mixed or even purely exc
% (e.g. 
%
% (*) if this is a merged file, cluster statistics will be biased by
% spikes from within triggers (right now the removal is file-specific)

function [ sst, fig ] = spikes_stats( filebase, varargin )

verbose                         = 1;

%-------------------- constants  --------------------%
nfpc                            = 3;                    % number of features/channel

padbuffer                       = [ -0.01 0.01 ];       % [s]
minCluToKeep                    = 2;                    % remove artefacts, noise
wintype                         = 'hanning';            % spectrum
nfft                            = 1024;                 % spectrum

MAXLAG                          = 0.050;                % ACH
CCHBS                           = 0.001;                % ACH

ISI1                            = 0.002;                % [s] count number of ISIs smaller than 2 ms
ISI2                            = 0.020;                % [s] out of ISIs smaller than 20 ms
DT_DEFAULT                      = 8;                    % [samples] detection dead time

USF                             = 32;                   % to increase spatial resolution
nbins                           = 100;                  % for firing rate stationarity analysis

BLOCKSIZE                       = 10 * 2^20;            % bytes/block; use 10 MB blocks
precision0                      = 'int16';              % parameters for block-wise loading of data
precision1                      = 'single';

BARCOLOR                        = [ 0 0 0 ];
CALIBWIDTH                      = 2;
CALIBCOLOR                      = [ 0 0.7 0 ];
ntoplot                         = 100;                  % waveform plot
color                           = [ 0.7 0.7 1 ];        % waveform plot
colorrand                       = [ 1 0 0 ];            % white noise spectrum plot
LW                              = 2;                    % plots
mcolor                          = zeros( 1, 3 );
mcolor( color == max( color ) ) = 1;

%-------------------- par, clu, res  --------------------%
% load/save
if isa( filebase, 'cell' ) && length( filebase ) == 2
    adate                       = filebase{ 1 };
    fnum                        = filebase{ 2 };
    params                      = get_date_parameters( adate );
    par                         = params.par;
    if fnum > 0
        mfname                  = sprintf( '%s%s.%s', params.fileprefix, params.date, num3str( fnum ) );
    else
        if isa( params.merged, 'char' )
            mfname              = params.merged;
        else
            mfnums              = setdiff( unique( params.merged( :, 2 ) ), 0 );
            if sum( mfnums == -fnum )
                mfname          = sprintf( '%s%s_%d', params.fileprefix, params.date, -fnum );
            end
        end
    end
    filebase                    = sprintf( '%s/%s/dat/%s/%s', params.drive, params.date, mfname, mfname );
end
if ~exist( filebase, 'file' ) && ~exist( fileparts( filebase ), 'dir' )
    error( 'missing filebase' )
end
[ shanknums, clunums, clustr, plotflag, Overwrite ...
    , periods, blockwise, byPar, doFET, DT, DTu ] = ParseArgPairs( ...
    { 'shanknums', 'clunums', 'clustr', 'plotflag', 'Overwrite' ...
    , 'periods', 'blockwise', 'byPar', 'doFET', 'DT', 'DTu' } ...
    , { [], [], 'clu', 0, -2 ...
    , -1, 1, 0, 1, DT_DEFAULT, DT_DEFAULT }...
    , varargin{ : } );
SaveFn                          = [ filebase '.sst' ];
if Overwrite < 0 && FileExists( SaveFn )
    if verbose
        fprintf( 1, '%s: Loading %s\n', upper( mfilename ), SaveFn )
    end
    load(SaveFn,'-MAT');
    fig                         = NaN;
    return
end
par                             = LoadXml( sprintf( '%s.xml', filebase ) );
if isempty( shanknums )
    shanknums                   = 1 : par.nElecGps;
end
nchans_all                      = zeros( 1, length( shanknums ) );
for si = shanknums
    nchans_all( si )            = length( par.SpkGrps( si ).Channels );
end
nsamples                        = par.SpkGrps( 1 ).nSamples;                % assume same for all spike groups
peaksample                      = par.SpkGrps( 1 ).PeakSample;              % assume same for all spike groups
Fs                              = par.SampleRate;
p2pvoltage                      = par.VoltageRange;
CorrFactor                      = ( ISI2 - DT / Fs ) / ( ISI1 - DT / Fs );
scalefactor                     = 1 / 2^par.nBits * p2pvoltage * 1e6 / par.Amplification;       % a2d units -> microVolts
switch wintype
    case 'hanning'
        win                     = hanning( nsamples );
    case 'rect'
        win                     = ones( nsamples, 1 );
    case 'hamming'
        win                     = hamming( nsamples );
end
if ~isa( clustr, 'char' )
    error( 'input type mismatch: CLUSTR' )
end

% clu and res files (with the 0,1 clusters to ensure matching with the spk files)
[ Res, Clu, Map ]               = LoadCluRes( filebase, shanknums, clunums, clustr, 0 );
if isempty( Res ) || sum( ismember( Clu, Map( Map( :, 3 ) >= minCluToKeep, 1 ) ) ) == 0
    fprintf( 1, 'No relevant spikes in %s, shanks %s\n', filebase, num2str( shanknums ) )
    sst                         = [];
    return
end
uclu = unique( Clu );

% select spikes
if periods( 1 ) < 0
    % remove spikes during stimulation
    [ Vals, Trigs ]             = LoadVals( filebase );                                         % load all val files
    if periods == -1 % any trigger
        tidx                    = true( size( Trigs ) );
    else % specific subset
        tidx                    = ismember( Trigs, abs( periods ) );
    end
    if ~isempty( tidx )
        uvals                   = dilutesegments( Vals( tidx, 1 : 2 ), 0, max( abs( padbuffer ) ) * Fs );
        uvals                   = uvals + ones( size( uvals, 1 ), 1 ) * [ floor( padbuffer( 1 ) * Fs ) ceil( padbuffer( 2 ) * Fs ) ];
    else
        uvals = [];
    end
end

% determine effective duration
if exist( [ filebase '.dat' ], 'file' )
    a                           = memmapfile( [ filebase '.dat' ], 'Format', 'int16' );
    nAll                        = length( a.Data ) / par.nChannels; % [samples]
    clear a
else
    a                           = memmapfile( [ filebase '.eeg' ], 'Format', 'int16' );
    nAll                        = round( length( a.Data ) / par.nChannels * par.SampleRate / par.lfpSampleRate ); % [samples]
    clear a
end
validPeriods                    = setdiffranges( [ 1 nAll ], uvals );
durValid                        = sum( diff( validPeriods, 1, 2 ) + 1 );

% set up the block-wise loading of spiking
% build the type casting string
precision                       = sprintf( '%s=>%s', precision0, precision1 );
% determine number of bytes/sample
a0                              = ones( 1, 1, precision0 );
a1                              = ones( 1, 1, precision1 );
sourceinfo                      = whos( 'a0' );
targetinfo                      = whos( 'a1' );
nbytes0                         = sourceinfo.bytes;
nbytes1                         = targetinfo.bytes;

% bins for stationarity analysis
fbin                            = round( ( max( Res ) - min( Res ) ) / ( nbins + 1 ) );
edges                           = min( Res ) : fbin : max( Res );
edges                           = edges( [ 1 : nbins; 2 : nbins + 1 ] )';
edges( nbins, 2 )               = edges( nbins, 2 ) + 1;

% initiliaze loop
sst.filebase                    = filebase;
shankclu                        = Map( :, 2 : 3 );
shankclu( shankclu( :, 2 ) < minCluToKeep, : ) = [];
sst.shankclu                    = shankclu;
sst.Tsec                        = ( max( Res ) - min( Res ) ) / Fs;
sst.freqs                       = [];
nclus                           = size( shankclu, 1 );

z                               = zeros( nclus, 1 );
sst.nspks                       = z; 
sst.mean                        = cell( nclus, 1 );
sst.sd                          = cell( nclus, 1 );
sst.spec                        = single( zeros( nfft / 2, nclus ) );
sst.ach                         = zeros( MAXLAG/CCHBS*2+1, nclus );
sst.max                         = single( zeros( nsamples, nclus ) );
sst.frateb                      = zeros( nbins, nclus );
sst.maxp2p                      = z;
sst.fmax                        = z;
sst.geo_com                     = single( z );
sst.geo_sd                      = single( z );
sst.geo_comh                    = single( z );
sst.geo_sdh                     = single( z );
sst.snr                         = single( z );
sst.ID                          = z;
sst.Lratio                      = z;
sst.frate                       = z;
sst.ISIindex                    = z;
sst.ISIratio                    = z;
sst.ach_com                     = z;
sst.hw                          = z;
sst.tp                          = z;
sst.asy                         = z;
sst.pyr                         = ~true( nclus, 1 );
sst.geo_fwhm                    = single( z );
sst.extremum                    = single( z );

% go over shanks
i                               = 0;
ushanks                         = unique( Map( uclu, 2 ) )';
for shanknum                    = ushanks
    
    %-------------------- allocate shank-specific - clu, res, fet, ind, spk indices--------------------%
    % determine the clu numbers to analyze from this shank
    shankclus                   = Map( Map( :, 2 ) == shanknum, 3 );
    shankclus                   = shankclus( shankclus >= minCluToKeep );
    shankclus                   = shankclus( : )';
    % general parameters
    nchans                      = nchans_all( shanknum );
    fprintf( 1, 'loading %s, shank %d, ', filebase, shanknum )
    if byPar
        ridx                    = resort( par.SpkGrps( shanknum ).Channels  );                  % reorder according to par file order
    else
        ridx                    = 1 : length( par.SpkGrps( shanknum ).Channels );
    end
    sourcefile = sprintf( '%s.spk.%d', filebase, shanknum );    
    if blockwise
        % determine number of spikes in the file
        spksize                 = nchans * nsamples * nbytes0;
        s                       = dir( sourcefile );
        spksource               = s.bytes / spksize;
    else
        % load the entire file at once:
        spks                    = readspk( sprintf( '%s.spk.%d', filebase, shanknum ), nchans, nsamples );
        spksource               = size( spks, 3 );
        spks                    = spks( ridx, :, : );
    end    
    fprintf( 1, '%d spikes', spksource )
    % keep the relevant (shank-specific) clu/res
    idx                         = ismember( Clu, Map( Map( :, 2 ) == shanknum, 1 ) );
    clu                         = Clu( idx );
    res                         = Res( idx );
    % load the channel map 
    indfname                    = sprintf( '%s.ind.%d', filebase, shanknum );
    byInd                       = 1;
    if exist( indfname, 'file' )
        ind                     = LoadInd( indfname );
    else
        byInd                   = 0;
    end
    % if clu is a partial list, must load to match spks
    if spksource ~= size( clu, 1 ) || byInd
        idxMap                  = 0;
        clu                     = LoadClu( sprintf( '%s.%s.%d', filebase, clustr, shanknum ) );
        res                     = LoadRes( sprintf( '%s.%s.%d', filebase, 'res', shanknum ) );
    else
        idxMap                  = 1;
    end
    % check channel map for consistency
    if byInd
        neurochans              = get_stimchans( filebase, [], 'neuronal' );
        indchans                = setdiff( unique( ind( : ) ), 0 ); % exclude 0
        if ~isempty( setdiff( indchans, neurochans ) )
            fprintf( '%s exists and loaded, but includes non-neuronal channels - ignored!\n', indfname )
            byInd               = 0;
        end
        if max( clu ) > size( ind, 1 )
            fprintf( '%s exists and loaded, but does not include mapping for some of the clu - ignored!\n', indfname )
            byInd               = 0;
        end
    end
    % load the FET files
    if doFET
        fets                    = LoadFet( sprintf( '%s.fet.%d', filebase, shanknum ) );
        fets( :, size( fets, 2 ) + ( -nchans : 0 ) ) = [];                                      % remove non-PCA features
        fets                    = single( fets );                                               % actually int16
    end
    % decide which spikes to keep:
    kidx                        = ~ismember( clu, 0 );                                          % remove 0 clusters (keep noise for fet)
    if periods( 1 ) < 0 && ~isempty( uvals )
        kidx( inranges( res, uvals ) ) = 0;                                                     % remove any spike that is in any segment defined by a val file
    elseif ~isempty( periods ) && size( periods, 2 ) == 2                                       % alternatively, keep only spikes during specified periods
        gidx                    = false( size( kidx ) );
        gidx( inranges( res, periods ) ) = 1;
        kidx                    = kidx & gidx;                                                  % logical
    end
    % also remove spikes that occurred twice within DT (can happen in template matching algorithms)
    ridx1                       = find( diff( res ) <= DTu );
    ridx2                       = sort( [ ridx1; ridx1 + 1 ] );
    ridxd                       = false( size( kidx ) );
    ridxd( ridx2 )              = 1;
    fprintf( 1, '; removing %d/%d d.t. periods and', sum( ~kidx ), length( kidx ) )
    kidx                        = kidx & ~ridxd;
    fprintf( 1, ' total of %d/%d', sum( ~kidx ), length( kidx ) )
    % actually dilute the clu/res/fet files
    res                         = res( kidx );
    clu                         = clu( kidx );
    if doFET
        fets                    = fets( kidx, : );
    end
    if ~blockwise
        spks                    = spks( :, :, kidx );
    end
    % now we have the clu, res, fet, (ind), and the indices of the proper spikes
    if byInd
        % check for a row of zeros, and if there is one, assume that these
        % are non-neuronal noise, and remove them from further analyses.
        isArtifact              = sum( ind, 2 ) == 0;
        rmv                     = ismember( clu, find( isArtifact ) );
        clu( rmv )              = [];
        res( rmv )              = [];
        fets( rmv, : )          = [];
        if ~blockwise
            spks( :, :, rmv )   = [];
        end
    end
    
    %-------------------- waveform analyses (blockwise) --------------------%
    if blockwise
        fkidx                   = find( kidx );
        % divide into blocks
        spkblk                  = floor( BLOCKSIZE / ( nchans * nsamples * nbytes1 ) );                     % number of spikes/block in RAM
        nblocks                 = ceil( spksource / spkblk );
        blocks                  = [ 1 : spkblk : spkblk * nblocks; spkblk : spkblk : spkblk * nblocks ]';
        blocks( nblocks, 2 )    = spksource;
        fprintf( 1, '; %d blocks; Block number ', nblocks )
        % initialize variables for accumulation
        nclusShank              = length( shankclus );
        nspks_CLU               = zeros( 1, nclusShank );
        mspks_SUM               = single( zeros( nchans, nsamples, nclusShank ) );
        mspks_SSQ               = single( zeros( nchans, nsamples, nclusShank ) );
        spec_SUM                = single( zeros( nfft / 2, nchans, nclusShank ) );
        snr_SUM                 = zeros( nchans, nclusShank );
        % open file for reading
        fp                      = fopen( sourcefile, 'r' );
        if fp == -1
            error( 'fopen error' )
        end
        % go over blocks
        for bidx                = 1 : nblocks
            fprintf( 1, '%d ', bidx )
            nspkBlock           = blocks( bidx, 2 ) - blocks( bidx, 1 ) + 1;
            startpos            = ( blocks( bidx, 1 ) - 1 ) * spksize;
            datasize            = [ nchans nsamples * nspkBlock ];
            rc                  = fseek( fp, startpos, 'bof' );
            if rc
                error( 'fseek error' )
            end
            spks                = fread( fp, datasize, precision );
            spks                = reshape( spks, [ nchans nsamples nspkBlock ] );                           % organize output
            spks                = spks( :, :, kidx( blocks( bidx, 1 ) : blocks( bidx, 2 ) ) );              % select subset of spikes and their corresponding labels
            spks                = spks( ridx, :, : );                                                       % reorder spikes according to par file
            cluBlock            = clu( find( fkidx >= blocks( bidx, 1 ), 1 ) : find( fkidx <= blocks( bidx, 2 ), 1, 'last' ) );
            % go over clusters and accumulate statistics
            ci                  = 0;
            for clunum          = shankclus
                ci              = ci + 1;
                if idxMap
                    idx         = cluBlock == find( ismember( Map( :, 2 : 3 ), [ shanknum clunum ], 'rows' ) );
                else
                    idx         = cluBlock == clunum;
                end
                nspks_CLU( ci ) = nspks_CLU( ci ) + sum( idx );
                if ~sum( idx )
                    continue
                end
                spk             = spks( :, :, idx );
                mspks_SUM( :, :, ci )   = mspks_SUM( :, :, ci ) + sum( spk, 3 ) * scalefactor;              % spike amplitude (average, then p2p)
                mspks_SSQ( :, :, ci )   = mspks_SSQ( :, :, ci ) + ( sum( spk .^ 2, 3 ) ) * scalefactor^2;
                for chan                = 1 : nchans
                    wav                 = squeeze( spk( chan, :, : ) ) * scalefactor;                       % spectrum (per spike, per chan)
                    [ pow2, f ]         = my_spectrum( wav, nfft, Fs, win, 0, 0 );
                    spec_SUM( :, chan, ci ) = spec_SUM( :, chan, ci ) + sum( pow2( 2 : end, : ), 2 );
                    % mean spike snr (per spike, per chan)
                    snr_SUM( chan, ci ) = snr_SUM( chan, ci ) + sum( calc_spk_snrs( shiftdim( spk( chan, :, : ), 1 ) ) );
                end
            end
        end
        fclose( fp );
    end
    fprintf( 1, '\n' )
    
    %-------------------- other statistics --------------------% 
    ci                          = 0;
    for clunum                  = shankclus
        ci                      = ci + 1;                                                           % in-shank counter
        i                       = i + 1;                                                            % global counter
        ri                      = find( ismember( Map( :, 2 : 3 ), [ shanknum clunum ], 'rows' ) ); % index in map
        %-------------------- allocate unit-specific  --------------------%
        if idxMap
            idx                 = clu == ri;
        else
            idx                 = clu == clunum;
        end
        nspks                   = sum( idx );       
        st                      = res( idx );
        if verbose
            fprintf( 1, 'computing for %d.%d (%d)', shanknum, clunum, nspks )
        end
        if ~nspks
            fprintf( 1, '\n' )
            continue
        end
        
        % summarize/compute waveform statistics for that cluster
        if blockwise
            % divide the accumulated statistics by the number of spikes in that cluster
            mspk                = mspks_SUM( :, :, ci ) / nspks_CLU( ci );                          % spike amplitude (average, then p2p)
            sspk                = sqrt( mspks_SSQ( :, :, ci ) / nspks_CLU( ci ) - mspk .^ 2 );      %SD = sqrt( <x^2> - <x>^2 );
            p2p                 = max( mspk, [], 2 ) - min( mspk, [], 2 );
            [ maxp2p, chan ]    = max( p2p );
            mp                  = spec_SUM( :, chan, ci ) / nspks_CLU( ci );                        % spectrum (per spike, on the largest amplitude channel, then average)
            snr                 = snr_SUM( chan, ci ) / nspks_CLU( ci );                            % mean spike snr (per spike, then average; for the peak amp channel only)
        else
            spk                 = spks( :, :, idx );
            mspk                = mean( spk, 3 ) * scalefactor;                                     % spike amplitude (average, then p2p)
            sspk                = std( spk, [], 3 ) * scalefactor;
            p2p                 = max( mspk, [], 2 ) - min( mspk, [], 2 );
            [ maxp2p, chan ]    = max( p2p );
            wav                 = squeeze( spk( chan, :, : ) ) * scalefactor;                       % spectrum (per spike, on the largest amplitude channel, then average)
            [ pow2, f ]         = my_spectrum( wav, nfft, Fs, win, 0, 0 );
            mp                  = mean( pow2( 2 : end, : ), 2 );
            snr                 = mean( calc_spk_snrs( shiftdim( spk( chan, :, : ), 1 ) ) );        % mean spike snr (per spike, then average)
        end
        [ ~, midx ]             = max( mp );
        fmax                    = f( midx );                                                        % peak of spectrum
        
        % spatial metrics
        sw                      = sum( mspk.^2, 2 );                                                % spatial signature of each waveform
        [ ~, midx ]             = max( sw );                                                        % determine if positive/negative
        ext                     = mspk( midx, peaksample );                                         % [microVolts]
        [ geo_com, geo_sd ]     = calc_com( ( 1 : nchans_all( shanknum ) )', sw );                  % spike COM and spatial dispersion
        [ geo_comh, geo_sdh ]   = calc_com( mod( ( 1 : nchans_all( shanknum ) )' + 1, 2 ), sw );    % spike COM and dispersion in horizontal direction
        if byInd                                                                                    % shift geo_com according to ind
            geo_com             = geo_com + ind( clunum, 1 ) - 1;
        end
        zu                      = fft_upsample( sw( : ), USF, 1 );                                  % geometrical dispersion by FWHM
        hm                      = ( max( zu ) - min( zu ) ) / 2 + min( zu );
        geo_fwhm                = sum( zu >= hm ) / USF;                                            % [sites]
        
        % morphological isolation
        dofet                   = 0;
        if doFET
            [ spkfets, nfets ]  = size( fets );
            fet                 = fets( idx, : );
            if sum( idx ) > nfets
                dofet           = 1;
            end
        end
        if dofet
            fidx                = find( sum( fet ) ~= 0 );
            if blockwise
                spkblk          = floor( BLOCKSIZE / ( nfets * nbytes1 ) );                         % number of spikes/block in RAM
                nblocks         = ceil( spkfets / spkblk );
                blocks          = [ 1 : spkblk : spkblk * nblocks; spkblk : spkblk : spkblk * nblocks ]';
                blocks( nblocks, 2 ) = spkfets;
                fprintf( 1, '; %d FET blocks; Block number ', nblocks )
                d2              = zeros( spkfets - sum( idx ), 1 );
                e0              = 0;
                for bidx        = 1 : nblocks
                    fprintf( 1, '%d ', bidx )
                    lidx        = false( spkfets, 1 );
                    lidx( blocks( bidx, 1 ) : blocks( bidx, 2 ) ) = 1;
                    if byInd                                                                        % reorganize fets as fetsI according to ind (idxMap is 0)
                        cluB        = clu( lidx );
                        fetsB       = fets( lidx, : );
                        fetsI       = single( zeros( sum( lidx ), size( fets, 2 ) ) );
                        chans0      = ind( clunum, : );
                        for clunumI = shankclus
                            idxI    = cluB == clunumI;
                            chansI  = ind( clunumI, : );
                            [ ~, i1, i2 ] = intersect( chans0, chansI );
                            cidx1   = ones( nfpc, 1 ) * ( i1( : ).' - 1 ) * nfpc + ( 1 : nfpc )' * ones( 1, length( i1 ) );
                            cidx1   = cidx1( : );
                            cidx2   = ones( nfpc, 1 ) * ( i2( : ).' - 1 ) * nfpc + ( 1 : nfpc )' * ones( 1, length( i2 ) );
                            cidx2   = cidx2( : );
                            fetsI( idxI, cidx1 ) = fetsB( idxI, cidx2 );
                        end
                        idxI        = cluB == clunum;
                        d2( e0 + ( 1 : sum( ~idxI ) ) ) = mahal( fetsI( ~idxI, fidx ), fet( :, fidx ) );
                        e0          = e0 + sum( ~idxI );
                    else
                        lidx( idx ) = 0;
                        d2( e0 + ( 1 : sum( lidx ) ) ) = mahal( fets( lidx, fidx ), fet( :, fidx ) );
                        e0          = e0 + sum( lidx );
                    end
                end
            else
                if byInd                                                                            % reorganize fets as fetsI according to ind (idxMap is 0)
                    fetsI       = single( zeros( size( fets ) ) );
                    chans0      = ind( clunum, : );
                    for clunumI = shankclus
                        idxI    = clu == clunumI;
                        chansI  = ind( clunumI, : );
                        [ ~, i1, i2 ] = intersect( chans0, chansI );
                        cidx1   = ones( nfpc, 1 ) * ( i1( : ).' - 1 ) * nfpc + ( 1 : nfpc )' * ones( 1, length( i1 ) );
                        cidx1   = cidx1( : );
                        cidx2   = ones( nfpc, 1 ) * ( i2( : ).' - 1 ) * nfpc + ( 1 : nfpc )' * ones( 1, length( i2 ) );
                        cidx2   = cidx2( : );
                        fetsI( idxI, cidx1 ) = fets( idxI, cidx2 );
                    end
                    fetI        = fetsI( idx, : );
                    d2          = mahal( fetsI( ~idx, fidx ), fetI( :, fidx ) );
                else
                    d2          = mahal( fets( ~idx, fidx ), fet( :, fidx ) );
                end
            end
            if sum( ~idx ) > sum( idx )
                sd2             = sort( d2 );
                ID              = sd2( sum( idx ) );
            else
                ID              = NaN;
            end
            Lratio              = sum( 1 - chi2cdf( d2, length( fidx ) ) ) / sum( idx );
        else
            ID                  = NaN;
            Lratio              = NaN;
        end
        
        % compute timing statistics
        [ cch, cchbins ]        = CCG( st, ones( size( st ) ), CCHBS * Fs, MAXLAG / CCHBS, Fs, 1, 'hz' );   % ACH parameters
        frate                   = nspks / durValid * Fs;                                                    % firing rate
        dt                      = diff( st ) / Fs;                                                          % ISI paramters
        ISIindex                = CorrFactor * sum( dt < ISI1 ) / sum( dt < ISI2 );                         % observed vs. expected spikes below ISI1
        ISIratio                = sum( dt < ISI1 ) / ( nspks - 1 );                                         % fraction of ISIs below ISI1
        
        % time-resolved firing rate:
        frateb                  = zeros( nbins, 1 );
        for j                   = 1 : nbins
            frateb( j, : )      = sum( st >= edges( j, 1 ) & st < edges( j, 2 ) ) / fbin * Fs;
        end
        
        %-------------------- plot  --------------------%
        if plotflag
            fig( i )            = figure;
            
            % plot all waveforms (scaled)
            subplot( 1, 3, 1 )
            plot_spk_waveforms( spk, ntoplot, color, 1, scalefactor, [ 1 1 0 ] );
            title( sprintf( 's%d.c%d waveforms', shanknum, clunum ) )
            ylabel( 'Amplitude ({\mu}V)' ), xlabel( 'Time (ms)' )
            lh                  = line( [ -0.5 -0.5 -0.25 ], [ -50 -100 -100 ] ); % 0.25 ms/50 uV calibration
            set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
            axis off
            
            % plot mean waveform for peak channel
            subplot( 2, 3, 2 )
            wfeatures( mspk( chan, : )', plotflag );
            xlabel( 'Time (samples)' )
            lh                  = line( nsamples * [ 0.25 0.25 0.5 ], [ -25 -75 -75 ] );
            set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
            ylabel( 'Amplitude ({\mu}V)' ),
            title( sprintf( 'Channel %d waveform', chan ) )
            
            % plot spectrum for that channel
            subplot( 2, 3, 5 )
            line( f( 2 : end ), mp, 'color', mcolor, 'Linewidth', LW );
            set( gca, 'xscale', 'log', 'xlim', [ 10 Fs/2 ] );
            ylabel( 'Power ({\mu}V^2)' ), xlabel( 'Frequency (Hz)' )
            title( sprintf( 'Spectrum; peak at %d Hz', round( f( midx ) ) ) )
            
            % overlay spectrum of white noise
            gnoise              = normrnd( 0, 1, nsamples, nspks );
            powrand             = mean( my_spectrum( gnoise, nfft, Fs, win, 0, 0 ), 2 );
            line( f( 2 : end ), powrand( 2 : end ) / max( powrand ) * diff( ylim )...
                , 'color', colorrand, 'LineStyle', '--', 'LineWidth', LW );
            
            % plot ACH
            subplot( 2, 3, 3 );
            cchbins             = [ fliplr( -CCHBS : -CCHBS : -MAXLAG )  0 : CCHBS : MAXLAG ];
            cchbins             = cchbins * 1000; % ms
            bh                  = bar( cchbins, cch, 1 );
            xlim( cchbins( [ 1 end ] ) - diff( cchbins( 1 : 2 ) ) * [ -1 1 ] )
            ylims               = [ 0 max( max( cch ), 1 ) * 1.2 ];
            ylim( ylims )
            set( bh, 'EdgeColor', BARCOLOR, 'FaceColor', BARCOLOR )
            ylabel( 'Spikes/s' )
            xlabel( 'Time (ms)' )
            set( gca, 'box', 'off', 'tickdir', 'out' )
            lh                  = line( [ -25 -15 ], mean( ylims ) * [ 1 1 ] ); % 10 ms calibration
            set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
            axis off
            title( sprintf( '%0.3g spikes/s', frate ) )
            
        end
        
        %-------------------- summarize --------------------%
        sst.nspks( i, : )       = nspks;
        sst.mean{ i, : }        = mspk;                     % all mean waveforms
        sst.sd{ i, : }          = sspk;                     % all SD waveforms
        sst.spec( :, i )        = mp;
        sst.ach( :, i )         = cch;
        sst.max( :, i )         = mspk( chan, : )';         % channel with peak
        sst.frateb( :, i )      = frateb;
        sst.maxp2p( i, : )      = maxp2p/1e3;               % max P2P (milliVolts)
        sst.fmax( i, : )        = fmax;                     % peak spectrum (max P2P chan)
        sst.geo_com( i, : )     = geo_com;                  % geometric COM
        sst.geo_sd( i, : )      = geo_sd;                   % spatial dispersion
        sst.geo_comh( i, : )    = geo_comh;                  % geometric COM in horizontal direction 
        sst.geo_sdh( i, : )     = geo_sdh;                   % spatial dispersion in horizontal direction 
        sst.snr( i, : )         = snr;                      % mean snr (max P2P chan)
        sst.ID( i, : )          = ID;                       % Isolation distance
        sst.Lratio( i, : )      = Lratio;                   % L-ratio
        sst.frate( i, : )       = frate;                    % mean firing rate
        sst.ISIindex( i, : )    = ISIindex;                 % ISI index (observed/expected)
        sst.ISIratio( i, : )    = ISIratio;                 % ISI ratio (observed/total)
        sst.geo_fwhm( i, : )    = geo_fwhm;                 % spatial dispersion by FWHM
        sst.extremum( i, : )    = ext;                      % signed value [microV] of peak-to-peak

        fprintf( 1, '\n' )
        
    end % clunum
    
end % shanknum

% compute other statistics
ach_com                         = calc_com( [ 0 : 1 : round( MAXLAG / CCHBS ) ], sst.ach( round( MAXLAG / CCHBS ) + 1 : end, : ) );
if plotflag
    fig( i + 1 )                = figure;
end
[ hw, asy, tp ]                 = wfeatures( sst.max, plotflag );
hw                              = hw / Fs * 1000;               % [ms]
tp                              = tp / Fs * 1000;               % [ms]

sst.ach_com                     = ach_com( : );                 % COM of ACH
sst.hw                          = hw( : );                      % waveform width at half-height (PYR wide)
sst.tp                          = tp( : );                      % trough-to-peak (PYR long)
sst.asy                         = asy( : );                     % waveform asymmetry (PYR asymmetric)
[ sst.pyr, sst.pval ]           = classify_waveform( [ sst.tp 1./sst.fmax * 1000 ] );   % logical vector for PYR/INT classification
if ~exist( 'f', 'var' )
    sst.freqs                   = [];
else
    sst.freqs                   = f( 2 : end );
end

if Overwrite || ~FileExists( SaveFn )
    if verbose
        fprintf( 1, 'Saving %s\n', SaveFn )
    end
    save( SaveFn, 'sst', '-v6' );
end

return

% EOF
