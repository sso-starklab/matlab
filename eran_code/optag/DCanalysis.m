% DCanalysis            determine optical response at every possible intensity and soma-diode distance
%
% call                  [ pstats, pinfo ] = DCanalysis( filebase )
% 
% optional arguments (given as name/value pairs)
%
%                       s                 {[]]}; spike structure, output of load_spikes
%                       ilevel            {'B'}; isolation level
%                       nmax              {inf}; maximal number of trials to use
%                       dxs               {0}; [inter-shank interval] range of distances (soma to diode)
%
%                       argument for get_triggers:
%                       stimType          { { 'PULSE', 'PSINE' } }
%                       freqRange         {[]};             
%                       valRange          {[ 0.06 0.07 ]}
%                       durRange          {[ 0.04 0.08 ]}
%                       uStim             {1}
%                       cmp               {'eq'}
%                       simOnly           {0}
%
%                       stim channel selection:
%                       wavRange          {[ -inf inf  ]}
%                       sourcetypes       {''}
% 
%                       summary parameters:
%                       toplot            {1}
%                       savetype          {'png'}
%                       figdir            {[]}
%                       matdir            {[]}
%                       Overwrite         {2}
%
% returns               pstats          [ nout nin rout rin pr NaN NaN lat pSup ]
%                       pinfo           [ shank clu celltype nt mdur dx ]
%                       
%
% calls                 LoadXml                                     (blab)
%                       get_stimchans, get_triggers, LoadStims      (formats) 
%                       calc_gain, makeedges, minmax, ParseArgPairs (general) 
%                       alines, textf, replacetok, fig_out          (graph)
%                       sortranges, setdiffranges                   (sets)
%                       determine_units, load_spikes, make_psth     (spikes)
%                       calc_index, inrange                         (stats)
%
% see also              celltypeClassification, opticalTagging, selectIntensity

% 13-may-13 ES

% revisiosn
% 28-may-13 wavelentgth, sourcetype, and p-suppression added
% 17-dec-19 cleaned up, removed SALT and SR1
% 22-dec-19 (1) defaultDurRange constant added
%           (2) if durRange differs from default, add duration to title.
%           this allows using multiple/non-default durations for the same channel

function [ pstats, pinfo ] = DCanalysis( filebase, varargin  )

% initialize output
pstats                  = [];
pinfo                   = [];

%---------------------------------------------------------------%
% constants
%---------------------------------------------------------------%
colors                  = [ 0 0 0.7; 1 0 0 ];
blackColor              = [ 0 0 0 ];
pTH                     = 0.01; % just for plotting
bias                    = 1e-2;
defaultDurRange         = [ 0.04 0.08 ];

%---------------------------------------------------------------%
% arguments
%---------------------------------------------------------------%
nargs                   = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ s, ilevel...
    , stimType, freqRange, valRange, durRange, wavRange, sourcetypes ...
    , uStim, cmp, simOnly...
    , nmax, dxs, nT, pbinsize ...
    , toplot, savetype, figdir, matdir, Overwrite ] = ParseArgPairs(...
    { 's', 'ilevel'...
    , 'stimType', 'freqRange', 'valRange', 'durRange', 'wavRange', 'sourcetypes' ...
    , 'uStim', 'cmp', 'simOnly'...
    , 'nmax', 'dxs', 'nT', 'pbinsize' ...
    , 'toplot', 'savetype', 'figdir', 'matdir', 'Overwrite' }...
    , { [], 'B' ...
    , { 'PULSE', 'PSINE' }, [], [ 0.06 0.07 ], defaultDurRange, [ -inf inf ], ''...
    , 1, 'eq', 0 ...
    , inf, 0, 5, 20 ...
    , 1, 'png', [], [], -2 }...
    , varargin{ : } );


par                     = LoadXml( filebase );
spkFs                   = par.SampleRate;
[ pathname, filename, extname ] = fileparts( filebase );
filename = [ filename extname ];

if isa( stimType, 'char' )
    stimstr             = stimType;
elseif isa( stimType, 'cell' ) && isa( stimType{ 1 }, 'char' )
    stimstr             = stimType{ 1 };
else
    error( 'stimType mismatch' )
end

%---------------------------------------------------------------%
% determine the channels, stims, spikes
%---------------------------------------------------------------%
% stimuli
[ stimchans, stimshanknums, ~, stimsources, stimwavelengths ] = get_stimchans( par );
kidx                    = inrange( stimwavelengths, wavRange );
if ~isempty( sourcetypes )
    kidx                = kidx & ismember( stimsources, sourcetypes );
end
stimchans( ~kidx )      = [];
stimshanknums( ~kidx )  = [];
if isempty( stimchans )
    fprintf( '%s: no relevant stim channels!\n', upper( mfilename ) );
    return
end

% optionally load the spikes too 
if isempty( s )
    shankclu            = determine_units( filebase, stimshanknums, ilevel );
else
    shankclu            = s.shankclu;
end

% determine the stim parameters
allvals                 = LoadStims( filebase );
vals                    = sortranges( allvals( :, 1 : 2 ) );
bperiods                = setdiffranges( [ 1 inf ], vals ); % at spkFs
params                  = { 'types', stimType, 'durs', durRange, 'vals', valRange, 'freqRange', freqRange };

% title name
epochname               = '';
if isempty( epochname )
    estr                = '';
else
    estr                = [ '_' epochname ];
end
shanknums               = unique( shankclu( :, 1 ) );
disstr                  = ( [ repmat( 'D', size( dxs( : ) ) ) num2str( dxs( : ) ) ]' );
disstr                  = disstr( : ).';
if isequal( disstr, 'D0D1D2D3' )
    disstr              = '';
else
    disstr              = [ '_' disstr ];
end
shstr                   = ( [ repmat( 'S', size( shanknums( : ) ) ) num2str( shanknums( : ) ) ]' );
shstr                   = sprintf( '%s_c%s', shstr( : ).', ilevel );
tstr                    = ( [ repmat( 'T', size( stimchans( : ) ) ) num2str( stimchans( : ) ) ]' );
tstr                    = tstr( : ).';
if isempty( freqRange )
    fstr                = '';
else
    mf                  = mean( freqRange, 2 );
    if mf( 1 ) > mf( 2 )
        f12             = round( [ max( freqRange( 1, : ) ) min( freqRange( 2, : ) ) ] );
    else
        f12             = round( [ min( freqRange( 1, : ) ) max( freqRange( 2, : ) ) ] );
    end
    fstr                = sprintf( '_f%df%d', f12( 1 ), f12( 2 ) );   
end
if isequal( durRange, defaultDurRange )
    durstr              = '';
else
    durstr              = sprintf( 't%dt%d_', round( durRange( 1 ) * 1000 ), round( durRange( 2 ) * 1000 ) );
end
stmstr                  = sprintf( 'v%0.3gv%0.3g_%s%su%d%s', valRange( 1 ), valRange( 2 )...
    , durstr, fstr, uStim, estr );
stmstr                  = replacetok( stmstr, '#', '.' );
dstr                    = sprintf( '%s_%s_%s%s', shstr, tstr, stmstr, disstr );
savebase                = sprintf( '%s.DC.%s', filename, stimstr );
corename                = sprintf( '%s.%s', savebase, dstr );

% saving target
basedir                 = fileparts( fileparts( pathname ) );
if isempty( figdir ) || isequal( figdir, 0 )
    figdir              = 0;
else
    if isa( figdir, 'char' ) && exist( figdir, 'dir' )
    else
        figdir          = sprintf( '%s/figs/wn', basedir );
    end
    figname             = [ figdir '/' corename ];
end
if isempty( matdir ) || isequal( matdir, 0 )
    matdir              = 0;
else
    if isa( matdir, 'char' ) && exist( matdir, 'dir' )
    else
        matdir          = sprintf( '%s/mat/wn', basedir );
    end
    savename            = [ matdir '/' corename ];
end

%---------------------------------------------------------------%
% actually compute
%---------------------------------------------------------------%
if Overwrite < 0 && exist( savename, 'file' )
    load( savename, 'filebase', 'pinfo', 'pstats', 'params', '-mat' )
    fprintf( '%s: loaded %s\n', upper( mfilename ), savename );
else
    if isempty( s )
        s               = load_spikes( filebase, stimshanknums, ilevel );
    end
    if isempty( s.res )
        return
    end
    bperiods( end )     = s.res( end ) + 1;
    
    pinfo               = [];
    pstats              = [];
    T                   = zeros( 1, 3 );
    
    for i               = 1 : length( stimchans )
        
        % get the stim times for a given channel
        trigchan        = stimchans( i ); % stimchans
        [ ~, ~, ~, ~, stims ] = get_triggers( filebase, trigchan, uStim, cmp, simOnly, params );
        if isempty( stims )
            continue
        end
        smat            = stims.times;
        if nmax < inf
            tidx        = 1 : min( [ size( smat, 1 ) nmax ] );
            smat        = smat( tidx, : );
        end
        nt              = size( smat, 1 );
        if nt == 0
            continue
        end
        mdur            = mean( ( diff( smat, 1, 2 ) + 1 ) / spkFs, 1 );
        
        % determine the relevant units
        shanknum        = stimshanknums( i );
        dxi             = abs( s.map( :, 2 ) - shanknum );
        
        % go over different distances
        for di          = 1 : length( dxs )
            
            dx          = dxs( di );
            umap        = s.map( dxi == dx, 1 );
            nclu        = length( umap );
            
            % base on an extended baseline period
            tout        = sum( diff( bperiods, 1, 2 ) + 1 );
            tin         = sum( diff( smat, 1, 2 ) + 1 );
            stts        = zeros( length( umap ), 9 );
            t1          = zeros( length( umap ), 3 );

            % go over different units
            for uidx = 1 : nclu
                
                t0      = clock;
                
                % prepare the spikes for the given unit
                res     = s.res( s.clu == umap( uidx ) );
                
                % compute the Poisson prob. to see nin or more spikes based on no-light periods
                nin     = length( inranges( res, smat ) );
                nout    = length( inranges( res, bperiods ) );
                lambda  = nout / tout * tin;
                pr      = 1 - poisscdf( nin - 1, lambda ); % prob to get %% NIN OR MORE %%
                pSup    = poisscdf( nin, lambda );
                rin     = nin / tin * spkFs;
                rout    = nout / tout * spkFs;
                
                % compute latency based on PSTH and global considerations
                win     = [ -1 1 ] * nT * mdur * spkFs;
                [ ~, ~, rr ] = make_psth( res, smat( :, 1 ), win, pbinsize, 0, 0 );
                nbins   = size( rr, 1 );
                r1      = full( rr( ceil( nbins / 2 ) : ( nbins - 1 ), : )' );
                h       = sum( r1, 1 );
                pvals   = 1 - poisscdf( h - 1, nout / tout * pbinsize );
                nsbins  = floor( mdur * spkFs / pbinsize );
                latAct  = find( pvals < pTH / nsbins, 1, 'first' ) * pbinsize / spkFs;
                if isempty( latAct )
                    latAct = NaN; 
                end

                t1( uidx, 1 ) = etime( clock, t0 );
                stts( uidx, : ) = [ nout nin rout rin pr NaN NaN latAct pSup ];
                
            end
            
            T           = T + sum( t1 );
            info        = [ s.shankclu( dxi == dx, : ) ones( nclu, 1 ) * [ nt mdur dx ] ];
            pinfo       = [ pinfo; info ];
            pstats      = [ pstats; stts ];
            
        end
        
    end
    
    if ~isequal( matdir, 0 ) && ( Overwrite == 1 || ( Overwrite ~= -1 && ~exist( savename, 'file' ) ) )
        save( savename, 'filebase', 'pinfo', 'pstats', 'params', '-v6' )
        fprintf( '%s: saved %s\n', upper( mfilename ), savename );
    end
end

if isempty( pinfo ) || ~toplot
    return
end

%---------------------------------------------------------------%
% plot
%---------------------------------------------------------------%

% plot and save
figtitle                = corename;

% generatl parameters
pidx                    = pinfo( :, 3 ) == 1;
rates0                  = pstats( :, 3 : 4 );
rates0( rates0 == 0 )   = bias;
rates                   = log10( rates0 ); % log scale
pvals                   = pstats( :, 5 );
sig                     = pvals < pTH | pvals > ( 1 - pTH );

% indices
midx                    = calc_index( pstats( :, [ 4 3 ] ) );
bins                    = -1 : 0.2 : 1;
edges                   = makeedges( bins );
nbins                   = length( bins );

fig                     = figure;

% rate ratio
for di                  = 1 : length( dxs )
    
    dx                  = dxs( di );
    didx                = pinfo( :, 6 ) == dx;
    
    subplot( 2, 2, di )
    
    % plot the dots
    hold on,
    cidx                = ~pidx & didx;
    if sum( cidx )
        ph( 1 )         = plot( rates( cidx, 1 ), rates( cidx, 2 ), '.' );
        set( ph( 1 ), 'color', colors( 1, : ) )
    end
    cidx                = pidx & didx;
    if sum( cidx )
        ph( 2 )         = plot( rates( cidx, 1 ), rates( cidx, 2 ), '.' );
        set( ph( 2 ), 'color', colors( 2, : ) )
    end
    cidx                = ~pidx  & didx & sig;
    if sum( cidx )
        ph( 1 ) = plot( rates(  cidx, 1 ), rates( cidx, 2 ), 'o' );
        set( ph( 1 ), 'color', colors( 1, : ) )
    end
    cidx                = pidx  & didx & sig;
    if sum( cidx )
        ph( 2 ) = plot( rates(  cidx, 1 ), rates( cidx, 2 ), 'o' );
        set( ph( 2 ), 'color', colors( 2, : ) )
    end
    
    % take care of axis scaling to log
    xlims               = xlim;
    ylims               = ylim;
    lims                = minmax( [ xlims( : ); ylims( : ) ] );
    if lims( 2 ) > 0
        lims( 2 )       = lims( 2 ) * 1.1;
    else
        lims( 2 )       = lims( 2 ) * 0.9;
    end
    if lims( 1 ) > 0
        lims( 1 )       = lims( 1 ) * 0.9;
    else
        lims( 1 )       = lims( 1 ) * 1.1;
    end
    set( gca, 'xlim', lims, 'ylim', lims )
    axis square
    ytick               = 10.^( floor( log10( 10^ min( get( gca, 'ytick' ) ) ) ) : ceil( log10( 10^ max( get( gca, 'ytick' ) ) ) ) );
    set( gca, 'ytick', log10( ytick ), 'yticklabel', ytick, 'xtick', log10( ytick ), 'xticklabel', ytick )
    line( lims, lims, 'linestyle', '--', 'color', [ 0 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' )
    if di == 1
        ylabel( 'Spikes/s (induced)' )
    end
    
    title( sprintf( '%d/%d INT/PYR, %d/%d shanks; %d trials, %0.3g ms (dx=%d)'...
        , sum( ~pidx & didx ), sum( pidx & didx )...
        , length( unique( pinfo( ~pidx & didx, 1 ) ) )...
        , length( unique( pinfo(  pidx & didx, 1 ) ) ) ...
        , round( mean( pinfo( didx, 4 ) ) )...
        , round( 1000 * mean( pinfo( didx, 5 ) ) ) ...
        , dx ...
        ) )

    % histogram of indices (all units)
    pos                 = get( gca, 'position' );
    newpos              = [ pos( 1 : 2 ) + pos( 3 : 4 ) * 0.1 pos( 3 : 4 ) * 0.9 ];
    set( gca, 'position', newpos );
    newsub              = [ pos( 1 ) + pos( 3 ) * 0.1 pos( 2 ) + pos( 4 ) * 0.1 - pos( 4 ) * 0.31 pos( 3 ) * 0.9 pos( 4 ) * 0.3 ];
    ah                  = axes( 'position', newsub );
    h0                  = histc( midx( ~pidx & didx ), edges ); 
    h1                  = histc( midx( pidx & didx ), edges ); 
    if isempty( h0 )
        h0              = zeros( nbins + 1, 1 );
    end
    if isempty( h1 )
        h1              = zeros( nbins + 1, 1 );
    end
    h0                  = h0( 1 : nbins ); 
    h1                  = h1( 1 : nbins );
    h                   = [ h0( : ) h1( : ) ];
    bh                  = bar( bins, bsxfun( @rdivide, h, sum( h ) ) );
    set( bh( 1 ), 'EdgeColor', colors( 1, : ), 'FaceColor', colors( 1, : ) )
    set( bh( 2 ), 'EdgeColor', colors( 2, : ), 'FaceColor', colors( 2, : ) )
    xlim( bins( [ 1 end ] ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    set( gca,'xticklabel', '' )
    alines( 0, 'x', 'color', blackColor, 'linestyle', '--' );
    if di == 1
        ylabel( 'Count' )
        xlabel( 'Modulation' )
    end

end

if di <= 2
    
    logscale            = 1;
    gain                = calc_gain( pstats( :, [ 4 3 ] ), [], 2 );
    lat                 = pstats( :, 8 ) * 1000;
    xlims               = [ 0.01 100 ];
    ylims               = [ 1 256 ];
    xsep                = 1; % gain of 1, i.e. no chance
    ysep                = 8; % 8 ms
    
    if logscale
        gain            = log10( gain );
        xlims           = log10( xlims ); 
        xsep            = log10( xsep );
        lat             = log2( lat );
        ylims           = log2( ylims );
        ysep            = log2( ysep );
    end
    
    for di              = 1 : length( dxs )
        dx              = dxs( di );
        didx            = pinfo( :, 6 ) == dx;
        
        subplot( 2, 2, di + 2 )
        
        % plot the dots
        hold on,
        cidx            = ~pidx & didx;
        if sum( cidx )
            ph( 1 )     = plot( gain( cidx, 1 ), lat( cidx, 1 ), '.' );
            set( ph( 1 ), 'color', colors( 1, : ) )
        end
        cidx            = pidx & didx;
        if sum( cidx )
            ph( 2 )     = plot( gain( cidx, 1 ), lat( cidx, 1 ), '.' );
            set( ph( 2 ), 'color', colors( 2, : ) )
        end
        cidx            = ~pidx  & didx & sig;
        if sum( cidx )
            ph( 1 )     = plot( gain(  cidx, 1 ), lat( cidx, 1 ), 'o' );
            set( ph( 1 ), 'color', colors( 1, : ) )
        end
        g( 1 )          = nanmedian( gain(  cidx, 1 ) );
        l( 1 )          = nanmedian( lat( cidx, 1 ) );
        cidx            = pidx  & didx & sig;
        if sum( cidx )
            ph( 2 )     = plot( gain(  cidx, 1 ), lat( cidx, 1 ), 'o' );
            set( ph( 2 ), 'color', colors( 2, : ) )
        end
        g( 2 )          = nanmedian( gain(  cidx, 1 ) );
        l( 2 )          = nanmedian( lat( cidx, 1 ) );
        
        if logscale
            g           = 10.^g;
            l           = 2.^l;
        end
        istr            = sprintf( 'g=%0.2g; %d ms', g( 1 ), round( l( 1 ) ) ); 
        pstr            = sprintf( 'g=%0.2g; %d ms', g( 2 ), round( l( 2 ) ) );
        
        xlims0          = xlim; 
        ylims0          = ylim;
        xlim( minmax( [ xlims0; xlims ] ) ); 
        ylim( minmax( [ ylims0; ylims ] ) );
        
        xtick           = 10.^( floor( log10( 10^ min( get( gca, 'xtick' ) ) ) ) : ceil( log10( 10^ max( get( gca, 'xtick' ) ) ) ) );
        set( gca, 'xtick', log10( xtick ), 'xticklabel', xtick )
        ytick           = 2.^( floor( log2( 2^ min( get( gca, 'ytick' ) ) ) ) : ceil( log2( 2^ max( get( gca, 'ytick' ) ) ) ) );
        set( gca, 'ytick', log2( ytick ), 'yticklabel', ytick )
        
        alines( xsep, 'x', 'color', blackColor, 'linestyle', '--' );
        alines( ysep, 'y', 'color', blackColor, 'linestyle', '--' );
        
        xlabel( 'Gain' ), 
        ylabel( 'Latency (ms)' ), 
        axis square, 
        set( gca, 'tickdir', 'out', 'box', 'off' )
        
        th              = text( min( xlim ) + 0.8 * diff( xlim ), min( ylim ) + 0.9 * diff( ylim ), istr );
        set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', colors( 1, : ) );
        th              = text( min( xlim ) + 0.8 * diff( xlim ), min( ylim ) + 0.8 * diff( ylim ), pstr );
        set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'color', colors( 2, : ) );
        
    end
    
end

textf( 0.5, 0.975, replacetok( figtitle, '\_', '_' ) );

if figdir ~= 0
    fig_out( fig, 1, [ figname '.' savetype ], savetype );
    fprintf( '%s: saved %s\n', upper( mfilename ), figname );
end

return

% EOF
