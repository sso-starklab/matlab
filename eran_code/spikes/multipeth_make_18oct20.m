% multipeth_make    mega-wrapper for multiple-unit/multiple event PETH
%
% [ peth, bins, trigs, tims, durs, shankclu, fig, ahx ] = multipeth_make( filebase )
%
% ARGUMENTS:
% filebase          full path + name
% 
% all other arguments are optional, specified by parameter name/value pairs:
% ilevel            {'B'}; see determine_units
% shanknums         {[]}; 
% data              {[]}; load_spikes.m structures to use data from memory
%
% can trigger on externally-specified events:
% periods           {[]}; [sec]; periods (start/end pairs) for which to compute PETH
% trigs             {[]}; identity (classification) of each period
% pooltrigs         {0}; if 1, pools all triggers together
%
% alternatively, can use stimuli:
% stimTypes         {'PULSE'}; if empty, will use only the external
%                       periods; if both are non-empty, will intersect
% valRange          {[ 0.001 6 ]}, [V]
% durRange          {[0.04 0.08}, [sec]
% freqRange         {[]}; [Hz]
% dfreqRange        {[]}; [Hz]
% channels          {[]}; see get_triggers
% uflag             {1}; flag passed to get_triggers; if on, use only
%                       triggers that appeared on a single channel
% multi             {'eq'}; comparison passed to get_triggers
% simonly           {0}; 1 considers only simultaneous events
% minNtrigs         {2}; minimal number of events/PETH
%
% parameters to control the histogram binning and duration
% nT                {2}; number of stimulus durations (symmetrically around onset)
%                           to generate non-symmetrically, nT should be a
%                           2-element vector
% nbins             {100}; number of bins in each histogram
% scale             {'hz'}
% sdGauss           {1}; bins
% 
% alternative control over histogram generation:
% Fs                {20000}, [Hz], sampling rate
% binsize           {[]}
% halfwin           {[]} number of bins, assummed symmetrical
%
% general parameters
% toplot            {1}
% savetype          {'png'}
% savef             {1}; save the plot in ../../figs/ or in the provided directory
%
% OUTPUT
% peth              3D array, nbins x nclu x ntrig; see multipeth.m
% bins              vector of times [sec] relative to trigger
% trigs             trigger labels; see get_triggers.m
% tims              trigger (onset) times [samples]
% durs              event durations [sec]
% shankclu          nclu x 3 corresponding to the columns of peth
% fig               handle to figure
%
% USAGE:
% to look at what's there to decide:
% >> [ vals chans stimsRaw ] = LoadStims( filebase ); 
% >> stim_plot( stimsRaw );
% 
% to plot multiple different PETH sets with the same set of units:
% >> s = load_spikes( filebase, shanknums, ilevel );
% >> [ peth bins trigs tims ] = multipeth_make( filebase, 'data', s );
%
% to get the corresponding spikes aligned and partitioned:
% >> r = get_rasters( s.clu, s.res, trigs, tims );
%
% calls             ParseArgPairs, replacetok, uhist                        (general)
%                   determine_units, get_spikes, multipeth, multipeth_plot  (spikes)
%                   sortranges                                              (sets)
%                   get_triggers                                            (formats)
%
% See also          get_rasters, load_spikes

% 26-feb-13 ES

% revisions
% 04-mar-13 argument handling improved, ZAP/sines supported
% 05-mar-13 external period determination added
% 01-may-13 pooling option added
% 06-may-13 periods support improved
% 21-oct-14 modified defaultfigdir
% 02-feb-16 ahx returned as well
% 18-aug-19 cleaned up
% 13-nov-19 added simonly to figure title and name (only when simonly ~= 0)
% 28-jan-20 added explanation of pooltrigs to help
% 18-oct-20 allow call with data only (without filebase)

% solve this:
% [ peth, bins, trigs, fig ] = multipeth_make( datenum2filebase( { 'm531r1', -32 } ), 'B', [], 'PULSE', [ 0.001 6 ], [ 0.04 0.08 ], 33 : 36, 1, 'eq' );

% to do: 
% (1) determine how to call s.t. only simultaneous events will be included
% (check in get_triggers)
% (2) modify nT argument s.t. it may be a non-symmetric window; in that
% case call multipeth/CCG with the max( abs( nT ) ), and modify before
% the multipeth_plot call

function [ peth, bins, trigs, tims, durs, shankclu, fig, ahx ] = multipeth_make( filebase, varargin )

%--------------------------- initialize ---------------------------%
peth                        = [];
bins                        = [];
trigs                       = [];
tims                        = [];
durs                        = [];
shankclu                    = [];
fig                         = [];
ahx                         = [];

%--------------------------- arguments ---------------------------%
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    filebase                = '';
end

[ ilevel, shanknums, data...
    , Fs, binsize, halfwin...
    , stimTypes, valRange, durRange, timeRange, freqRange, dfreqRange, channels, uflag, multi, simonly, minNtrigs...
    , periods, trigs, pooltrigs...
    , nbins, scale, sdGauss, nT...
    , toplot, savetype, savef ] = ParseArgPairs(...
    { 'ilevel', 'shanknums', 'data'...
    , 'Fs', 'binsize', 'halfwin'...
    , 'stimTypes', 'valRange', 'durRange', 'timeRange', 'freqRange', 'dfreqRange', 'channels', 'uflag', 'multi', 'simonly', 'minNtrigs'...
    , 'periods', 'trigs', 'pooltrigs' ...
    , 'nbins', 'scale', 'sdGauss', 'nT'...
    , 'toplot', 'savetype', 'savef' }...
    , { 'B', [], []...
    , 20000, [], []...
    , 'PULSE', [ 0.001 6 ], [ 0.04 0.08 ], [ 0 inf ], [], [], [], 1, 'eq', 0, 2 ...
    , [], [], 0 ...
    , 100, 'hz', 1, 2 ...
    , 1, 'png', 1 } ...
    , varargin{ : } );
clear varargin
if ( isempty( periods ) || isempty( data ) ) && isempty( filebase )
    return
end

%--------------------------- events ---------------------------%
if ~isempty( periods ) && ( isempty( stimTypes ) || isequal( stimTypes, 0 ) )
    
    periods                 = sortranges( periods );                    % assume already in s
    nperiods                = size( periods, 1 );
    if isempty( trigs )
        trigs               = ones( nperiods, 1 );
    end
    periods                 = round( periods * Fs );                    % [samples]
    durs                    = ( diff( periods, [], 2 ) + 1 ) / Fs;      % [s]
    tims                    = round( periods( :, 1 ) );                 % [samples]
    vals                    = zeros( nperiods, 1 );
    pstr                    = sprintf( 'Tex_K%d_n%d_t%d', length( unique( trigs ) )...
        , nperiods, round( median( durs ) * 1000 ) );

elseif isequal( stimTypes, 0 )
    return
    
else
        
    if ~isa( stimTypes, 'cell' ) && isa( stimTypes, 'char' )
        stimTypes           = { stimTypes };
    end
    params                  = { 'types', stimTypes, 'durs', durRange, 'vals', valRange...
        , 'franges', freqRange, 'dfranges', dfreqRange, 'times', timeRange };
    
    % get those triggers
    [ trigs, tims, durs, vals ] = get_triggers( filebase, channels, uflag, multi, simonly, params );
    if pooltrigs
        trigs = ones( size( trigs ) );                                  % pool triggers together
    end
    
    % prune triggers (those with too few repetitions)
    [ ntrig, utrig ]            = uhist( trigs );
    ridx                        = ismember( trigs, utrig( ntrig < minNtrigs ) );
    trigs( ridx, : )            = [];
    tims( ridx, : )             = [];
    durs( ridx, : )             = [];
    vals( ridx, : )             = [];
    
    % keep only events with onset during specified periods (e.g. ripples, animal location, brain state, ..)
    if ~isempty( periods )
        periods                 = sortranges( periods );            % assume already in sec
        periods                 = round( periods * Fs );            % [samples]
        kidx                    = inranges( tims, periods );
        trigs                   = trigs( kidx, : );
        tims                    = tims( kidx, : );
        durs                    = durs( kidx, : );
        vals                    = vals( kidx, : );
    end
    
    if isempty( trigs )
        fprintf( '%s: No such triggers for %s\n', upper( mfilename ), filebase )
        return
    end
    
    % build a corresponding parameter string
    pstr                        = '';
    for i                       = 1 : length( channels )
        pstr                    = sprintf( '%sT%d', pstr, channels( i ) );
    end
    pstr                        = sprintf( '%sU%dC%s_', pstr, uflag, multi );
    for k                       = 1 : length( stimTypes )
        pstr                    = sprintf( '%s%s_', pstr, stimTypes{ k } );
    end
    pstr( end )                 = [];
    pstr                        = sprintf( '%s_v%0.3gv%0.3g_t%dt%d', pstr, valRange( 1 ), valRange( 2 )...
        , round( durRange( 1 ) * 1000 ), round( durRange( 2 ) * 1000 ) );
    pstr                        = replacetok( pstr, '#', '.' );
    if simonly ~= 0
        pstr                    = sprintf( '%s_s%d', pstr, simonly );
    end

end

%--------------------------- spikes ---------------------------%
if isa( data, 'struct' ) && all( isfield( data, { 'shankclu', 'clu', 'res', 'map' } ) )
    shankclu                    = data.shankclu;
    ilevel                      = data.ilevel;
    clu                         = data.clu;
    res                         = data.res;
    map                         = data.map;
    if isempty( shankclu )
        return
    end
    if ~isempty( shanknums )
        if isvector( shanknums ) && ( numel( shanknums ) > 2 || size( shanknums, 2 ) == 1 )
            shanks              = shanknums( : );
            shankclus           = [];
        elseif size( shanknums, 2 ) == 2 || size( shanknums, 2 ) == 3
            shanks              = unique( shanknums( :, 1 ) );
            shankclus           = shanknums( :, 1 : 2 );
        else
            shanks              = [];
            shankclus           = [];
        end
        if ~isempty( shankclus )
            kidx                = ismember( shankclu( :, 1 : 2 ), shankclus, 'rows' );
        elseif ~isempty( shanks )
            kidx                = ismember( shankclu( :, 1 ), shanks );
        else
            kidx                = [];
        end
        if ~isempty( kidx )
            aidx                = ismember( clu, map( kidx, 1 ) );
            clu                 = clu( aidx );
            res                 = res( aidx );
            map                 = map( kidx, : );
            shankclu            = shankclu( kidx, : );
        end
    end
else
    shankclu                    = determine_units( filebase, shanknums, ilevel );
    [ clu, res, map ]           = get_spikes( filebase, shankclu );
    [ ~, ~, i2 ]                = intersect( map( :, 2 : 3 ), shankclu( :, 1 : 2 ), 'rows' );
    shankclu                    = shankclu( i2, : );
end

% generate a string for this:
shanknums                       = unique( shankclu( :, 1 ) );
sstr                            = '';
for i                           = 1 : length( shanknums )
    sstr                        = sprintf( '%sS%d', sstr, shanknums( i ) );
end
sstr                            = sprintf( '%s_c%s', sstr, ilevel );

%--------------------------- peth ---------------------------%
% compute the peth such that they are symmetrical around onset time
% binsize = round( 2 * nT * median( durs ) / nbins * Fs );  % samples
% halfwin = nbins / 2; % bins
if ~isempty( binsize ) && ~isempty( halfwin )
    nT                          = [ -1 1 ];
    shiftPlot                   = 0;
elseif ~isempty( nT )
    if length( nT ) == 1
        nT                      = [ -1 1 ] * nT;
        shiftPlot               = 1;
    else
        shiftPlot               = 0;
    end
    binsize                     = round( diff( nT ) * median( durs ) / nbins * Fs );
    halfwin                     = ceil( max( abs( nT ) ) * 2 / diff( nT ) * nbins / 2 );
else
    error( 'missing parameters' )
end

clucat                          = [ map( :, 1 ) shankclu( :, 3 ) ];

[ peth, bins ]                  = multipeth( clu, res, trigs, tims...
    , 'binsize', binsize, 'halfwin', halfwin, 'scale', scale...
    , 'sdGauss', sdGauss, 'clucat', clucat );

%--------------------------- plot ---------------------------%
if ~toplot
    return
end
if shiftPlot
    % now shift the plots s.t. they are symmetrical around onset/offset mean
    kidx                        = ( nbins - nbins / diff( nT ) * ( diff( nT ) - 1 ) ) + 1 : length( bins );
    kidx                        = round( kidx );
    kpeth                       = peth( kidx, :, : );
    kbins                       = bins( kidx );
else
    kpeth                       = peth;
    kbins                       = bins;
end

% generate a string for the figure:
[ pathname, filename, extname ] = fileparts( filebase );
filename                        = [ filename extname ];
homedir                         = strfind( pathname, 'dat' );
if isempty( homedir )
    homedir                     = [ fileparts( pathname ) '/' ];
else
    homedir                     = pathname( 1 : homedir - 1 );
end
defaultfigdir                   = sprintf( '%sfigs', homedir );
figtitle                        = [ filename '_' sstr '_' pstr ];
if isequal( savef, 1 )
    if ~exist( defaultfigdir, 'dir' )
        mkdir( homedir, 'figs' )
    end
    figname                     = [ defaultfigdir '/' figtitle ];
elseif isa( savef, 'char' ) && exist( savef, 'dir' )
    figname                     = [ savef '/' figtitle ];
else
    figname                     = '';
end

[ fig, ahx ]                    = multipeth_plot( kpeth, kbins, shankclu...
    , 'trigs', trigs, 'durs', durs, 'vals', vals...
    , 'str', figtitle, 'figname', figname, 'savetype', savetype );

return

% EOF
