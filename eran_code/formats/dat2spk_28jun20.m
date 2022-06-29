% DAT2SPK               detrend and make *.spk.# file from *dat; optionally subsample
%
% CALL                  [ RC, SPK ] = DAT2SPK( SOURCEFNAME, SPKFNAME, RES, NCHANNELS, CHANNELS )
%
% GETS                  SOURCEFNAME     full path to binary file, e.g. *dat
%                       SPKFNAME        full path to target file
%                       RES             at samples of SOURCEFNAME (res)
%                       NCHANNELS       number of channels in SOURCEFNAME
%                       CHANNELS        for which to extract the waveforms (1-based)
%
% OPTIONAL ARGUMENTS (given as name/value pairs):
%
%                       NSAMPLES        {32}; for each waveform
%                       PEAKSAMPLE      {16}; sample number to coincide with SPIKETIMES
%                       DETRENDFLAG     {1}; whether or not to detrend
%                       VERBOSE         {0}; whether to report to STDOUT
%                       CLU             {[]}; identity of each spike
%                       IND             {[]}; source channels for each spike (1-based)
%
% DOES                  -loads the spikes at the prescribed times (res)
%                       -optionally: subsamples channels according to clu/ind
%                       -optionally: removes DC and trends
%                       -saves the spikes to the target file
%
% CALLS                 mydetrend, readbin, ParseArgPairs
%
% NOTES                 1. uses a procesing buffer of 1000 spikes
%                       2. channel order is kept as in the vector CHANNELS (if IND not specified)
%                               or as in the array IND (if specified)

% 21-apr-19 ES

% revisions
% 27-jul-19 (1) modified argument handling
%           (2) extended to extract different channels for different spikes
% 20-may-20 renamed from make_spk_file to dat2spk
% 28-jun-20 (1) bug fix for ind

function [ rc, spk ] = dat2spk( sourcefname, spkfname, res, nchannels, channels, varargin )

% constants
source                      = 'int16';      % how data are stored
target                      = 'single';     % how data are used in MATLAB workspace
blocksize                   = 1000;         % spikes; e.g. 10 x 32 x 1000 x 4 bytes = 1.25 MB

% default values
PEAKSAMPLE                  = 16;           % [samples]
NSAMPLES                    = 32;           % [samples]

% arguments
nout                        = nargout;
nargs                       = nargin;
if nargs < 5 || isempty( sourcefname ) || isempty( spkfname ) || isempty( res ) || isempty( nchannels ) || isempty( channels ) 
    return
end
[ nsamples, peaksample ...
    , clu, ind ...
    , detrendflag, verbose ] = ParseArgPairs(...
    { 'nsamples', 'peaksample' ...
    , 'clu', 'ind' ...
    , 'detrendflag', 'verbose' } ...
    , { NSAMPLES, PEAKSAMPLE ...
    , [], [] ...
    ,  1, 0 } ...
    , varargin{ : } );

% output
rc                          = 0;
spk                         = [];

% file names
if ~exist( sourcefname, 'file' )
    fprintf( 1, '%s: missing file %s, existing\n', upper( mfilenme ), sourcefname )
    return
end

% load parameters
nchans                      = length( channels );

% get spike times
nspikes                     = size( res, 1 ); 
if nspikes == 0
    fprintf( 1, '%s: no spikes\n', upper( mfilenme ), sourcefname )
    return
end
periods                     = res * [ 1 1 ] + ones( nspikes, 1 ) * [ -peaksample + 1 nsamples - peaksample ];

% determine channel allocations
if ~isempty( clu ) && ~isempty( ind )
    nsClu                   = size( clu, 1 );
    nsRes                   = size( res, 1 );
    if nsClu ~= nsRes
        clu                 = [];
    end
    uclu                    = unique( clu );
    nclu                    = length( uclu );
    nzcols                  = sum( ind ) ~= 0;
    ind                     = ind( :, nzcols );
    nind                    = sum( nzcols );
    if nclu ~= nind 
        ind                 = [];
    end
    uchans                  = unique( ind( : ) );
    echans                  = setdiff( uchans, channels );
    if ~isempty( echans )
        ind                 = [];
    end
end
if ~isempty( clu ) && ~isempty( ind )
    subsample               = 1;
    nichans                 = size( ind, 1 );
else
    subsample               = 0;
    nichans                 = nchans;
end

% work in blocks
t0                          = clock;
if verbose
    fprintf( 1, '%s: %d spikes', spkfname, nspikes )
end
nblocks                     = ceil( nspikes / blocksize );
fp                          = fopen( spkfname, 'wb' );
if fp == -1
    fprintf( 'error opening file %s for writing', spkfname )
    rc                      = rc + 1;
    return
end
for i                       = 1 : nblocks
    if verbose
        fprintf( 1, '.' )
    end
    bidx                    = ( ( i - 1 ) * blocksize + 1 )  : min( i * blocksize, nspikes );
    nspikesb                = length( bidx );
    % load spikes
    if i == 1 && periods( bidx( 1 ), 1 ) < 0 % first spike
        spk0                = readbin( datfname, channels, nchannels, [ 1 periods( bidx( 1 ), 2 ) ], source, target );
        spk1                = [ zeros( nchans, 1 - periods( bidx( 1 ), 1 ), target ) spk0 ];
        spk                 = readbin( datfname, channels, nchannels, periods( bidx( 2 : end ), : ), source, target );
        spk                 = reshape( spk, [ nchans nsamples nspikesb - 1 ] );
        spk                 = cat( 3, spk1, spk );
    else % any other spike
        spk                 = readbin( sourcefname, channels, nchannels, periods( bidx, : ), source, target );
        spk                 = reshape( spk, [ nchans nsamples nspikesb ] );
    end
    % subsample according to clu and ind
    if subsample
        sidx                = clu( bidx );
        cidx                = ind( :, sidx );
        tmp                 = zeros( nichans, nsamples, nspikesb );
        for j               = 1 : nspikesb
            tmp( :, :, j )  = spk( cidx( :, j ), :, j );
        end
        spk                 = tmp;
    end
    % detrend
    if detrendflag
        for ci              = 1 : nichans
            spk( ci, :, : ) = mydetrend( squeeze( spk( ci, :, : ) ) );
        end
    end
    % write to binary file
    spk                     = reshape( spk, [ nichans nsamples * nspikesb ] );
    count                   = fwrite( fp, spk, source );
    if count                ~= ( nspikesb * nsamples * nichans )
        rc                  = rc + 2;
        return
    end
end
fclose( fp );
if verbose
    fprintf( 1, 'done (%0.3g sec)\n', etime( clock, t0 ) )
end

% load the detrended spikes
if nout > 1
    precision               = sprintf( '%s=>%s', source, target );
    fp                      = fopen( spkfname, 'r' );
    spk                     = fread( fp, [ nichans inf ], precision );
    fclose( fp );
    nspks                   = size( spk, 2 ) / nsamples;
    spk                     = reshape( spk, [ nichans nsamples nspks ] );
end

return

% EOF
