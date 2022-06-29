% opticalTagging            determines activation/silencing during stim 
%
% call                      [ pstats pinfo peth bins ] = opticalTagging( Clu, Res, Map, Stim )
% 
% receives                  Clu       identities (output of load_spikes)
%                           Res       [samples] (output of load_spikes)
%                           Map       cross-ref matrix (output of load_spikes)
%                           Stim      [samples], periods for each stimulus
%
% optional arguments (name/value pairs)
% 
%                           blackout        {[]}; periods to exclude from baseline computation; defaults to same as a padBuffered Stim
%                           spkFs           {20000}; sampling rate of Res, Stim, and blackout
%                           padBufferSec    { [ -0.01 0.01 ] }
%
%                           nbins           {100}, for PSTH and latency computations
%                           nT              {5} number of stim periods to expand before and after stimulus onset (used for latency estimation)
%                           ispyr           { all 1 }; for graphical purposes only
%
%                           savetype        { 'png' }
%                           filebase        for saving etc
%                           corename        {''}; to be added to saved files
%                           Overwrite       { -2 }
%                           graphics        { 1 }
%
% note          this routine does not load data from disk - all data appear as
%               arguments, and it is the user's responsibility to supply them properly
%
% returns 
%           pstats:               matrix with the following row for each unit:
%                                     [ nout nin rout rin pAct pSup latAct ]
%                                     nout:   number of spikes during baseline
%                                     nin:    number of spikes during stim
%                                     rout:   spike rate during baseline
%                                     rin:    spike rate during stim
%                                     pAct:   Poisson probability to see nin or more
%                                     pSup:   Poisson probability to see nin or less
%                                     latAct: latency [s] of first rate in/decrease
%
%           pinfo:                matrix with the following row for each unit:
%                                     [ shank clu resptype nt mdur ]
%                                     resptype:   1 if activated (pAct < 0.001)
%                                                -1 if silenced (pSup < 0.001)
%                                     nt:         number of stimuli
%                                     mdur:       mean stim duration
%
%           peth, bins:           for each unit
%
%
% For baseline, this routine uses the entire range during which the Stim periods were applied.
% Specifically, without the stims themselves (expanded to the sides by padbuffer) and
% without blacked-out periods (e.g. other stims, also expanded to the sides). 
%
% calls:                calc_gain, calc_sem, ParseArgPairs (general)
%                       barwerror, alines, textf, replacetok, fig_out (graph)
%                       inranges, setdiffranges, uniteranges (sets) 
%                       multipeth, multipeth_plot poissonTest (spikes)
%                       firfilt, makegaussfir (ssp)
%                       bounds (stats)
%
% see also              celltypeClassification, DCanalysis, selectIntensity

% 27-sep-13 ES

% revisions
% 15-nov-13 many modifications (rewritten)
% 11-dec-19 problem with hcolorbar handle, removed
% 15-dec-19 partial cleanup
% 16-dec-19 finalized cleanup
% 17-dec-19 nbins moved to argument list
% 19-dec-19 ~empty( fig ) check added

function [ pstats, pinfo, peth, bins ] = opticalTagging( Clu, Res, Map, Stim, varargin )

% constants
pTH                         = 0.01;             % for latency estimation (will be bonf. corr)
pTHclass                    = 0.001;            % for post-hoc classification
minSpikes                   = 100;              % cannot determine optical response with less than that

% initialize output
pstats                      = [];
pinfo                       = [];
peth                        = [];
bins                        = [];

% arguments
nargs = nargin;
if nargs < 4 || isempty( Clu ) || isempty( Res ) || isempty( Map ) || isempty( Stim )
    error( 'critical arguments missing!!\n' )
end
[ blackout, ispyr, nT, nbins ...
    , spkFs, padBufferSec ...
    , savetype, filebase, corename, Overwrite, graphics ] = ParseArgPairs(...
    { 'blackout', 'ispyr', 'nT', 'nbins' ...
    , 'spkFs', 'padBufferSec'...
    , 'savetype', 'filebase', 'corename', 'Overwrite', 'graphics' }...
    , { [], true( size( Map, 1 ), 1 ), 5, 100 ...
    , 20000, [ -0.01 0.01 ]...
    , 'png', '', '', -2, 1 }...
    , varargin{ : } );
ispyr( isnan( ispyr ) )     = 1; % NaN->PYR
ispyr                       = logical( ispyr );

% directories for figures, data
filename                    = [];
if ~isempty( filebase ) && exist( fileparts( filebase ), 'dir' )
    delim                   = strfind( filebase, '/dat/' );
    if isempty( delim )
        fprintf( '%s: Cannot save fig and/or data\n', upper( mfilename ) )
    end
    if isa( graphics, 'char' ) && exist( graphics, 'dir' )
        figdir              = graphics;
        graphics            = 1;
    else
        figdir              = [ filebase( 1 : delim ) 'figs/dc' ];
        if ~exist( fileparts( figdir ), 'dir' )
            mkdir( fileparts( fileparts( figdir ) ), 'figs' )
        end
        if ~exist( figdir, 'dir' )
            mkdir( fileparts( figdir ), 'dc' )
        end
    end
    matdir                  = [ filebase( 1 : delim ) 'mat/dc' ];
    if ~exist( fileparts( matdir ), 'dir' )
        mkdir( fileparts( fileparts( matdir ) ), 'mat' )
    end
    if ~exist( matdir, 'dir' )
        mkdir( fileparts( matdir ), 'dc' )
    end
    [ ~, filename, extname ] = fileparts( filebase );
    filename                = [ filename extname ];
else
    figdir                  = '';
    matdir                  = '';
    Overwrite               = 0;
end

% names of saved data and figure files
if isempty( corename )
    corename                = sprintf( '%s.opticalTagging', filename );
else
    corename                = sprintf( '%s.opticalTagging', corename );
end
if ~isempty( figdir ) && exist( figdir, 'dir' )
    figname                 = [ figdir '/' corename ];
else
    figname                 = '';
end
if ~isempty( matdir ) && exist( matdir, 'dir' )
    savename                = [ matdir '/' corename ];
else
    savename                = '';
end

if Overwrite < 0 && exist( savename, 'file' )
    fprintf( '%s: loading %s...', upper( mfilename ), savename )
    load( savename, 'filebase', 'pstats', 'pinfo', 'peth', 'bins', '-mat' );
    return
end

%----------------------------------------------------------------%
% preps
%----------------------------------------------------------------%

padBuffer                   = [ floor( padBufferSec( 1 ) * spkFs ) ceil( padBufferSec( 1 ) * spkFs ) ];
durs                        = diff( Stim, 1, 2 ) + 1;
stimBuf                     = bsxfun( @plus, Stim, padBuffer );
if isempty( blackout )
    excludedRanges          = uniteranges( stimBuf );
else
    blackBuf                = bsxfun( @plus, blackout, padBuffer );
    excludedRanges          = uniteranges( stimBuf, blackBuf );
end
mdur                        = mean( durs );
bperiods                    = setdiffranges( minmax( stimBuf ), excludedRanges );
tout                        = sum( diff( bperiods, 1, 2 ) + 1 );
tin                         = sum( diff( Stim, 1, 2 ) + 1 );
nclu                        = size( Map, 1 );

fprintf( '%d events; tin: %0.3g sec; tout: %0.3g sec\n', size( Stim, 1 ), tin / spkFs, tout / spkFs )

%----------------------------------------------------------------%
% analysis
%----------------------------------------------------------------%

% compute PSTH for each unit
clucat                      = [ Map( :, 1 ) ispyr ];
binsize                     = ceil( mdur * nT / nbins );
ntrig                       = size( Stim, 1 );
trg                         = ones( ntrig, 1 );

t0                          = clock;
[ peth, bins ]              = multipeth( Clu, Res, trg, Stim( :, 1 )...
    , 'binsize', binsize, 'halfwin', nbins, 'Fs', spkFs...
    , 'graphics', 0, 'clucat', clucat, 'scale', 'count' );
et1                         = etime( clock, t0 );

% preparations for Poisson tests and bonferroni corrections
idxIn                       = inranges( Res, Stim );
idxOut                      = inranges( Res, bperiods );
nsbins                      = ceil( mdur / binsize ); 
pTHcorr                     = pTH / nsbins;

% actually go over units and compute:
sh                          = size( peth( :, 1 ) );
pstats                      = zeros( nclu, 7 );
surprise                    = zeros( nbins * 2 + 1, nclu );
t0 = clock;
for uidx                    = 1 : nclu
    fprintf( '.' )
    
    % compute the Poisson prob. to see nin or more spikes based on no-light periods
    clunum                  = Map( uidx, 1 );
    nin                     = sum( Clu( idxIn ) == clunum );
    nout                    = sum( Clu( idxOut ) == clunum );
    [ pAct, pSup ]          = poissonTest( nout / tout, nin, tin );
    rin                     = nin / tin * spkFs;
    rout                    = nout / tout * spkFs;
    
    % compute latency based on PSTH and global considerations (this assumes similar durations of all stims)
    h                       = peth( :, uidx );
    [ pvalsInc, pvalsDec ]  = poissonTest( repmat( nout / tout, sh ), h, repmat( binsize * ntrig, sh ) );
    latAct                  = find( pvalsInc( bins >= 0 ) < pTHcorr | pvalsDec( bins >= 0 ) < pTHcorr, 1, 'first' ) * binsize / spkFs; %[s]
    if isempty( latAct )
        latAct              = NaN;
    end
    
    % keep the results
    surprise( :, uidx )     = log10( ( pvalsDec + eps )./ ( pvalsInc + eps ) );
    pstats( uidx, : )       = [ nout nin rout rin pAct pSup latAct ];
    
end
fprintf( '\n' )

et2                         = etime( clock, t0 );
fprintf( 'Total time : %0.3g+%0.3g sec\n' , et1, et2  )

% summarize pinfo
nt                          = ones( nclu, 1 ) * length( durs );
resptype                    = NaN * ones( nclu, 1 );
resptype( pstats( :, 5 ) <= pTHclass ) = 1;
resptype( pstats( :, 6 ) <= pTHclass ) = -1;
resptype( sum( pstats( :, 1 : 2 ), 2 ) < minSpikes ) = NaN;
pinfo                       = [ Map( :, 2 : 3 ) resptype nt mdur * ones( nclu, 1 ) ];

if Overwrite == 1 || Overwrite < 0 && ~exist( savename, 'file' )
    fprintf( '%s: saving %s...\n', upper( mfilename ), savename )
    save( savename, 'filebase', 'pstats', 'pinfo', 'peth', 'bins' );
end

%----------------------------------------------------------------%
% graphics
%----------------------------------------------------------------%

if ~graphics
    return
end

% plot the summary PSTH:
gwin                = makegaussfir( 1, 1 );
pethHz              = firfilt( peth / ntrig / binsize * spkFs, gwin );
[ fig, ah ]         = multipeth_plot( pethHz, bins, [ pinfo( :, 1 : 2 ) clucat( :, 2 ) ], 'trigs', trg, 'durs', durs/spkFs, 'vals', trg, 'savetype', NaN );

% shift to upper left:
pos                 = get( ah( 1 ), 'position' );
set( ah( 1 ), 'position', [ pos( 1 ) pos( 2 ) + pos( 4 ) / 2 pos( [ 3 4 ] ) / 2 ] );
pos                 = get( ah( 2 ), 'position' );
set( ah( 2 ), 'position', [ pos( 1 ) pos( 2 ) + pos( 4 ) pos( [ 3 4 ] ) / 2 ] );
pos                 = get( ah( 3 ), 'position' );
set( ah( 3 ), 'position', [ pos( 1 ) - pos( 3 ) / 2 pos( 2 ) + pos( 4 ) / 2 pos( [ 3 4 ] ) / 2 ] );
pos                 = get( ah( 4 ), 'position' );
set( ah( 4 ), 'position', [ pos( 1 ) - pos( 3 ) / 2 pos( 2 ) + pos( 4 ) pos( [ 3 4 ] ) / 2 ] );

% summarize the effects anatomically
edges               = [ -inf -1.05 : 0.1 : 1.05 inf ];
gbins               = -1.1 : 0.1 : 1.1;
ushanks             = unique( pinfo( :, 1 ) );
h                   = zeros( length( edges ), max( ushanks ) );
gain                = log10( calc_gain( pstats( :, [ 4 3 ] ) ) );

for i               = ushanks.'
    sidx            = pinfo( :, 1 ) == i;
    h( :, i )       = histc( gain( sidx ), edges  );
    mm( :, i )      = [ nanmean( gain( sidx & ispyr == 0 ) ) nanmean( gain( sidx & ispyr == 1 ) ) ]';
    ss( :, i )      = [ calc_sem( gain( sidx & ispyr == 0 ) ) calc_sem( gain( sidx & ispyr == 1 ) ) ]';
    nn( :, i )      = [ nansum( sidx & ispyr == 0 ) nansum( sidx & ispyr == 1 ) ]';
end
h( end, : )         = [];

ah( 5 )             = axes( 'position', [ 0.575 0.7125 0.425 * 3 / 4 0.425 / 2 ] );
subplot( ah( 5 ) )
imagesc( 1 : size( h, 2 ), gbins, scale( h ) ), axis xy
alines( -1 : 1, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
set( gca, 'ytick', -1 : 1, 'yticklabel', 10.^( -1 : 1 ) )
set( gca, 'tickdir', 'out', 'box', 'off' )

ah( 6 )             = axes( 'position', [ 0.575 0.5 0.425 * 3 / 4 0.425 / 2 ] );
subplot( ah( 6 ) )
barwerror( 1 : size( h, 2 ), mm( 1, : ), ss( 1, : ), [ 0 0 0.7 ], 0.8, [ 0 0 0 ], 0.5 );
hold on
barwerror( 1 : size( h, 2 ), mm( 2, : ), ss( 2, : ), [ 1 0 0   ], 0.8, [ 0 0 0 ], 0.5 );

for i               = 1 : size( h, 2 )
    text( i, min( ylim ) * 0.8, sprintf( '%d', nn( 1, i ) ), 'HorizontalAlignment', 'center', 'color', [ 0 0 0.7 ] );
    text( i, max( ylim ) * 0.8, sprintf( '%d', nn( 2, i ) ), 'HorizontalAlignment', 'center', 'color', [ 1 0 0 ] );
end
ylims               = ylim;
alines( -1 : 1, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
set( gca, 'ytick', -1 : 1, 'yticklabel', 10.^( -1 : 1 ) )
ylim( ylims );
set( gca, 'tickdir', 'out', 'box', 'off' ),
xlabel( 'Shank#' ), ylabel( 'Gain' )

sf                  = firfilt( surprise, gwin );
for celltypeidx = 0 : 1
    uidx = ispyr == celltypeidx;
    if ~sum( uidx )
        continue
    end
    if celltypeidx == 0
        ah( 7 )     = axes( 'position', [ 0.075       0.075       0.2125       0.2125 ] );
    else
        ah( 8 )     = axes( 'position', [ 0.2875       0.075       0.2125       0.2125 ] );
    end
    str             = sprintf( '%d,%d/%d (sup,act/tot)'...
        , sum( pinfo( uidx, 3 ) == -1 ), sum( pinfo( uidx, 3 ) == 1 ), sum( uidx ) );
    mat             = sf( :, uidx );
    clim            = bounds( mat( : ), 0.01 );
    imagesc( bins, 1 : sum( uidx ), sf( :, uidx )', clim( : ).' )%, axis xy,
    sh              = pinfo( uidx, 1 );
    borders         = find( diff( sh ) > 0 ) + 0.5;
    shnums          = sh( [ 1; ceil( borders ) ] );
    alines( borders, 'y', 'color', [ 1 0 0 ], 'linewidth', 1 );
    set( gca, 'tickdir', 'out', 'box', 'off' ),
    alines( [ 0 mdur / spkFs ], 'x', 'color', [ 1 0 0 ], 'linestyle', '--' );
    set( gca, 'ytick', [], 'yticklabel', [] )
    if celltypeidx == 0
        xlabel( 'Time [sec]' )
    end
    title( str )
end

textf( 0.5, 0.975, replacetok( corename, '\_', '_' ) );

%----------------------------------------------------------------%
% save the graphics
%----------------------------------------------------------------%
if ~isempty( savetype ) && ( isa( savetype, 'cell' ) || ~all( isnan( savetype ) ) ) && ~isempty( figname ) && ~isempty( fig )
    if ~isa( savetype, 'cell' )
        savetype = { savetype };
    end
    for j = 1 : length( savetype )
        figsave = [ figname '.' savetype{ j } ];
        fig_out( fig, 1, figsave, savetype{ j } );
        fprintf( '%s: saved %s\n', upper( mfilename ), figsave )
    end
end

return

% EOF
