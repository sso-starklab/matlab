% hfoAnalysisSpikingCount       ripple-triggerd PSTH and CC between spike count and ripple power.
% 
% call                          [ stats, fig ] = hfoAnalysisSpikingCount( filebase )
% 
% does                          computes ripple-triggerd PSTH and CC between spike count and ripple power
%
% return                        stats           structure
%                               fig             handle to single figure
%
% calls                         (blab)        CCG, LoadXml
%                               (formats)     LoadStims, LoadVals
%                               (general)     ParseArgPairs, replacetok, sortcols
%                               (graph)       fig_out, imagescbar, myjet, patch_band, replacetok, textf 
%                               (lfp)         rips_select
%                               (sets)        inranges, resampleranges, setdiffranges
%                               (spikes)      burstf, load_spikes, multipeth, warpedpeth
%                               (stats)       calc_spearman, inrange
%                               (ssp)         firfilt, makegaussfir
% 
% see also                      hfoAnalysisSpiking

% 24-sep-13 ES

% revisions
% 29-oct-13 batch processing supported - stats output, Overwrite option
% 10-nov-13 ACH computation added
% 23-dec-13 warped PETH added
% 17-aug-19 cleaned a bit and renamed hfoAnalysisSpikingCountNew
% 11-mar-21 cleaned properly and renamed back as hfoAnalysisSpikingCount

function [ stats, fig ] = hfoAnalysisSpikingCount( filebase, varargin )

stats                               = [];
fig                                 = [];

nargs                               = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ channel, graphics, ilevel...
    , minISI, binsizeSEC, halfwinSEC...
    , nbinsCenter, nbinsPre...
    , offsets, windur...
    , achbinsizeSEC, achhalfwinSEC, savetype, Overwrite ...
    ]                               = ParseArgPairs(...
    { 'channel', 'graphics', 'ilevel'...
    , 'minISI', 'binsizeSEC', 'halfwinSEC'...
    , 'nbinsCenter', 'nbinsPre' ...
    , 'offsets', 'windur'...
    , 'achbinsizeSEC', 'achhalfwinSEC', 'savetype', 'Overwrite',...
    }...
    , { [], 0, 'B'...
    , 0.3, 0.001, 0.2...
    , 4, 12 ...
    , -0.2 : 0.025 : 0.2, 0.05...
    , 0.001, 0.03, 'png', -2 ...
    }...
    , varargin{ : } );

delim                               = strfind( filebase, '/dat/' );
if isempty( delim )
    fprintf( '%s: Cannot save fig and/or data\n', upper( mfilename ) )
end
if isa( graphics, 'char' ) && exist( graphics, 'dir' )
    figdir                          = graphics;
    graphics                        = 1;
else
    figdir                          = [ filebase( 1 : delim ) 'figs/hfo' ];
    if ~exist( fileparts( figdir ), 'dir' )
        mkdir( fileparts( fileparts( figdir ) ), 'figs' )
    end
    if ~exist( figdir, 'dir' )
        mkdir( fileparts( figdir ), 'hfo' )
    end
end
matdir                              = [ filebase( 1 : delim ) 'mat/hfo' ];
if ~exist( fileparts( matdir ), 'dir' )
    mkdir( fileparts( fileparts( matdir ) ), 'mat' )
end
if ~exist( matdir, 'dir' )
    mkdir( fileparts( matdir ), 'hfo' )
end
[ ~, filename, extname ]            = fileparts( filebase );
filename                            = [ filename extname ];

%------------------------------------------------------------------------
% load data
%------------------------------------------------------------------------
par                                 = LoadXml( filebase );

% select a single ripple channel
sps                                 = load( [ filebase '.sps' ], '-mat' );
if isempty( channel )
    vidx                            = find( inrange( sps.stats( :, 5 ), [ 110 210 ] ) );
    [ ~, maxidx ]                   = max( sps.stats( vidx, 7 ) );
    channel                         = sps.stats( vidx( maxidx ), 2 );
end

if isempty( channel )
    return
end

figname                             = [ figdir '/' filename '.hfo_PSTH_ampCC.' num3str( channel ) ];
savename                            = [ matdir '/' filename '.hfo_PSTH_ampCC.' num3str( channel ) ];

% check saved
if Overwrite < 0 && exist( savename, 'file' )
    fprintf( '%s: loading %s...\n', upper( mfilename ), savename )
    load( savename, '-mat', 'stats' )
    return
end
load( [ filebase '.spw.' num3str( channel ) ], '-mat', 'rips' );

% load spikes
shanknums                           = [];
s                                   = load_spikes( filebase, shanknums, ilevel );
if isempty( s.map )
    return
end
clucat                              = [ s.map( :, 1 ) s.shankclu( :, 3 ) ];

vals                                = LoadStims( filebase );
if isempty( vals )
    vals                            = LoadVals( filebase );
end
vals                                = resampleranges( vals( :, 1 : 2 ), par.lfpSampleRate, par.SampleRate );

%------------------------------------------------------------------------
% select the ripples
%------------------------------------------------------------------------

% remove ripples during stim
[ ~, idx ]                          = setdiffranges( rips.edges, vals );
sidx                                = false( size( rips.trigs ) );
sidx( idx )                         = 1;

% also burst filter ripples
[ ~, ridx ]                         = burstf( rips.trigs, minISI * par.lfpSampleRate );
sidx( ridx )                        = 0;

rips                                = rips_select( rips, sidx );
trigs                               = resampleranges( rips.trigs, par.SampleRate, par.lfpSampleRate );
mat                                 = resampleranges( rips.edges, par.SampleRate, par.lfpSampleRate );
amp                                 = rips.sd;
nrips                               = size( mat, 1 );
onoff                               = mean( mat - trigs * [ 1 1 ] ) / par.SampleRate;

%------------------------------------------------------------------------
% compute all PSTH:
%------------------------------------------------------------------------
% no plotting, no smoothing, 200 ms before/after, with 1 ms bins:
fprintf( '%s: Computing peri-ripple PETHS...\n', upper( mfilename ) )
binsize                             = binsizeSEC * par.SampleRate;
halfwin                             = halfwinSEC / binsizeSEC;
[ peth, bins ]                      = multipeth( s.clu, s.res, ones( size( trigs ) ), trigs...
    , 'scale', 'hz', 'binsize', binsize, 'halfwin', halfwin...
    , 'Fs', par.SampleRate, 'clucat', clucat );

%------------------------------------------------------------------------
% compute all warped PSTH:
%------------------------------------------------------------------------
% no plotting, no smoothing, 9 bins during ripple, 12 before/after:
fprintf( '%s: Computing warped peri-ripple PETHS...\n', upper( mfilename ) )
sts                                 = load( [ filebase '.sts.sws' ] ); % already at eegFs
sts                                 = resampleranges( sts, par.SampleRate, par.lfpSampleRate ); % @ spkFs
wstats                              = warpedpeth( s.clu, s.res, trigs, mat, 'nbinsCenter', nbinsCenter, 'nbinsPre', nbinsPre...
    , 'Fs', par.SampleRate, 'map', s.map( :, 1 ), 'sts', sts, 'vals', vals );

%------------------------------------------------------------------------
% compute all ACH during ripples:
%------------------------------------------------------------------------
fprintf( '%s: Computing ripple ACHs...\n', upper( mfilename ) )
mat                                 = resampleranges( rips.edges, par.SampleRate, par.lfpSampleRate );
achbinsize                          = achbinsizeSEC * par.SampleRate;
achhalfwin                          = achhalfwinSEC / achbinsizeSEC;
idx                                 = inranges( s.res, mat );
clu                                 = s.clu( idx );
uclu                                = unique( clu );
nclu                                = length( uclu );
ccg                                 = CCG( s.res, s.clu, achbinsize, achhalfwin, par.SampleRate, uclu, 'count', mat );
ach                                 = zeros( 2 * achhalfwin + 1, nclu );
for i                               = 1 : nclu
    ach( :, i )                     = ccg( :, i, i );
end

%------------------------------------------------------------------------
% compute trial-by-trial correlation with ripple magnitude (SD)
%------------------------------------------------------------------------
fprintf( '%s: Computing spike-ripple correlations...\n', upper( mfilename ) )

% spike count (during ripples, can be shifted) of each unit during each event:
uclu                                = unique( s.clu );
edges                               = [ uclu - 0.5 ; uclu( end ) + 0.5 ];
nclu                                = length( uclu );
noffsets                            = length( offsets );
acc                                 = zeros( nclu, noffsets );
for j                               = 1 : noffsets
    cmat                            = [ mat( :, 1 ) + ceil( offsets( j ) * par.SampleRate ) mat( :, 1 ) + ceil( ( offsets( j ) + windur ) * par.SampleRate ) ];
    [ idx, seg ]                    = inranges( s.res, cmat );
    useg                            = unique( seg );
    count                           = zeros( nrips, nclu );
    for i                           = 1 : length( useg )
        aseg                        = useg( i );
        aidx                        = idx( seg == aseg );
        h                           = histc( s.clu( aidx ), edges );
        count( aseg, : )            = h( 1 : nclu );
    end
    cc                              = calc_spearman( count, amp * ones( 1, nclu ) );
    acc( :, j )                     = cc;
end

%------------------------------------------------------------------------
% summarize output, save
%------------------------------------------------------------------------
stats                               = struct( 'filename', filename );
stats.filename                      = repmat( { filename }, [ nclu, 1 ] );
stats.nevents                       = nrips * ones( nclu, 1 );
stats.shankclu                      = s.shankclu;
stats.bins                          = bins;                                 % [s]
stats.onoff                         = onoff;                                % [s]
stats.acc                           = acc;
stats.peth                          = peth';
stats.wgain                         = wstats.wgain;
stats.wbins                         = wstats.wbins;
stats.phi                           = wstats.phi;
stats.plo                           = wstats.plo;
stats.mi                            = wstats.mi;
stats.gain                          = wstats.gain;

if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( savename, 'file' ) )
    fprintf( '%s: saving %s...\n', upper( mfilename ), savename )
    save( savename, 'stats', '-v6' )
end

%------------------------------------------------------------------------
% plot
%------------------------------------------------------------------------
if graphics
    [ peth, bins ]                  = multipeth( s.clu, s.res, ones( size( trigs ) ), trigs...
        , 'scale', 'hz', 'binsize', binsize, 'halfwin', halfwin...
        , 'Fs', par.SampleRate, 'graphics', 0, 'clucat', clucat, 'sdGauss', 1 );
    gwin                            = makegaussfir( 1, 1 );
    peth                            = firfilt( peth, gwin );
    
    fig                             = figure;
    subplot( 2, 2, 1 ),
    mathat                          = [];
    for ct                          = 0 : 1
        uidx                        = s.shankclu( :, 3 ) == ct;
        mat                         = sortcols( peth( :, uidx ) );
        mathat                      = [ mathat mat ];
    end
    %imagesc( bins, 1 : size( mathat, 2 ), mathat' ), axis xy
    bords                           = [ sum( s.shankclu( :, 3 ) == 0 ) sum( ismember( s.shankclu( :, 3 ), [ 0 1 ] ) ) ] + 0.5;
    [ ~, ah ]                       = imagescbar( bins, 1 : size( mathat, 2 ), mathat );
    tstr                            = sprintf( '%s, ch%d, %d events', filename, channel, nrips );
    subplot( ah( 1 ) ), 
    ylabel( 'Cell #' )
    alines( onoff, 'x', 'color', [ 1 0 0 ] );
    alines( bords, 'y', 'color', [ 1 0 0 ] );
    subplot( ah( 2 ) )
    xlabel( 'Spikes/s' )
    alines( bords, 'y', 'color', [ 1 0 0 ] );

    subplot( 2, 2, 2 ),
    wbins                           = stats.wbins;
    mat                             = stats.wgain';
    mat                             = firfilt( mat, [ 0.25 0.5 0.25 ] );    % 3-bin smoother
    mathat                          = [];
    for ct                          = 0 : 1
        uidx                        = s.shankclu( :, 3 ) == ct;
        mathat                      = [ mathat sortcols( mat( :, uidx ) ) ];
    end
    bords                           = [ sum( s.shankclu( :, 3 ) == 0 ) sum( ismember( s.shankclu( :, 3 ), [ 0 1 ] ) ) ] + 0.5;
    [ ~, ah ]                       = imagescbar( wbins, 1 : size( mathat, 2 ), mathat );
    subplot( ah( 1 ) ), 
    ylabel( 'Cell #' )
    alines( [ onoff + diff( wbins( 1 : 2 ) ) * [ -1 1 ] / 2 0 ], 'x', 'color', [ 1 1 1 ], 'linestyle', '--' );
    alines( bords, 'y', 'color', [ 1 0 0 ] );
    subplot( ah( 2 ) ), xlabel( 'Gain' )
    alines( bords, 'y', 'color', [ 1 0 0 ] );
    alines( 1, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );

    subplot( 2, 1, 2 )
    patch_band( offsets, mean( acc( s.shankclu( :, 3 ) == 0, : ), 1 ) ...
        , calc_sem( acc( s.shankclu( :, 3 ) == 0, : ), 1 ), [ 0 0 0.7 ] );
    patch_band( offsets, mean( acc( s.shankclu( :, 3 ) == 1, : ), 1 ) ...
        , calc_sem( acc( s.shankclu( :, 3 ) == 1, : ), 1 ), [ 1 0 0 ] );
    xlabel( 'Offset (50 ms window relative to ripple peak power) [s]' )
    ylabel( 'CC (spike count vs. ripple amp, SD)' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    alines( onoff, 'x', 'color', [ 1 0 0 ] );
    pos                             = get( ah( 1 ), 'position' );
    pos1                            = get( gca, 'position' );
    set( gca, 'position', [ pos1( 1 : 2 ) pos( 3 ) pos1( 4 ) ] )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    textf( 0.5, 0.975, replacetok( tstr, '\_', '_' ) )
    colormap( myjet )
    
    if ~isempty( savetype ) && ~all( isnan( savetype ) )
        if ~isa( savetype, 'cell' )
            savetype = { savetype };
        end
        for j = 1 : length( savetype )
            for i = 1 : length( fig )
                fig_out( fig( i ), 1, [ figname '.part.' num2str( i ) '.' savetype{ j } ], savetype{ j } );
            end
        end
    end
    
end

return

% EOF