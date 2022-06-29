% selectIntensity   driving current that yields largest number of affected units
%
% call              [ slevel, stats, sstats ] = selectIntensity( filebase )
%
% optional arguments (name/value pairs)
% 
%                   celltype        { 1 }; 1 means PYR, 0 means PV-INT
%                   supFlag         { 0 }; 1 means activated, 0 silenced
%                   ilevel          { 'B' }
% 
%                   get_stimchans parameters (to select channels):
%                   stimchans       { [] }
%                   wavRange        { [ -inf inf ] }
%                   sourcetypes     { 'LED' }
%                   
%                   get_trigger parameters (to select stimuli):
%                   stimType        { 'PULSE', 'PSINE' }
%                   uStim           { 0 }
%                   durRange        { [ 0.04 0.08 ] }
%                   disRange        { 0 }; 0 means same-shank
%                   valRange        { [ 0 inf ] }
%
%                   minSpikes       { 10 }
%                   pTH             { 0.01 }
%                   toplot          { 0 }; by DCanalysis
%                   graphics        { 0 }; by this function
%                   Overwrite       { -2 }
%
% returns           slevel          range, in A
%                   stats           detailed information (see code)
%                   sstats          [ n gain latency ] (activated units)
%
% calls             LoadXml (blab)
%                   get_stimchans, get_stimlevels (formats)
%                   calc_gain, ParseArgPairs (general)
%                   myjet, replacetok (graph)
%                   DCanalysis (optag)
%                   isoverlap (sets)
%                   load_spikes (spikes)
% 
% called by         celltypeClassification, gwnAnalysisGetSegments
%
% see also          DCanalysis, opticalTagging

% 22-may-13 ES

% revisions
% 26-may-13 added valRange, defaults to [ 0 0.07 ]
% 01-oct-14 added support for 'camkii::et/tc'
% 16-apr-19 added support for 'pv::chr2 pv::jaws'
% 15-dec-19 removed argument TypeName and replaced with celltype and supFlag
% 17-dec-19 corrected typo in final selection stage (ct->celltype)
% 19-dec-19 added argument stimchans

% to do: 
% (1) udnerstand why in low-intensitied DCanalysis number of cells is lower (ntot)
% (2) integrate the ustim=-1
% (5) partition the LD/LED data properly. right now we lose about half of
% the data by choosing the optimal intensity for one source type. more
% generally, the optimization should be done by shanks (sources). but there
% is a limit to the complexity..

function [ slevel, stats, sstats ] = selectIntensity( filebase, varargin )

%------------------------------------------------------------------
% i/o
%------------------------------------------------------------------
slevel                  = [];
stats                   = [];
sstats                  = [];

nargs = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ celltype, supFlag, ilevel...
    , stimType, graphics, uStim, durRange, disRange, valRange ...
    , stimchans, wavRange, sourcetypes ...
    , minSpikes, pTH ....
    , toplot, Overwrite ] = ParseArgPairs(...
    { 'celltype', 'supFlag', 'ilevel'...
    , 'stimType', 'graphics', 'uStim', 'durRange', 'disRange', 'valRange' ...
    , 'stimchans', 'wavRange', 'sourcetypes' ...
    , 'minSpikes', 'pTH' ...
    , 'toplot', 'Overwrite' }...
    , { 1, 0, 'B' ...
    , { 'PULSE', 'PSINE' }, 0, 1, [ 0.04 0.08 ], 0, [ 0 inf ] ...
    , [], [ -inf inf ], 'LED' ...
    , 10, 0.01 ...
    , 0, -2 }...
    , varargin{ : } );
if ~ismember( celltype, [ 0 1 ] )
    error( 'celltype must be 0 or 1' )
end
switch celltype
    case 0
        cstr            = 'INT';
    case 1
        cstr            = 'PYR';
end
if ~ismember( supFlag, [ 0 1 ] )
    error( 'supFlag must be 0 or 1' )
end

%------------------------------------------------------------------
% set up - general
%------------------------------------------------------------------

% paths
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

if isa( sourcetypes, 'char' )
    sourcetypes = { sourcetypes };
end

% load the spikes
if Overwrite >= 0
    s                   = load_spikes( filebase, [], ilevel );
else
    s                   = [];
end

% check for relevant channels
par                     = LoadXml( filebase );
[ stimchans, ~, ~, stimsources, stimwavelengths ] = get_stimchans( par, stimchans );
kidx                    = inrange( stimwavelengths, wavRange );
if ~isempty( sourcetypes )
    kidx                = kidx & ismember( stimsources, sourcetypes );
end
stimchans( ~kidx )      = [];
stimsources( ~kidx )    = [];

% abort if no stim channels or if more than one source type
if isempty( stimchans )
    fprintf( '%s: no relevant stim channels!\n', upper( mfilename ) );
    return
end
usources                = unique( stimsources );
srsstr                  = usources{ 1 };
if length( usources ) > 1
    fprintf( '%s: NOTE: multiple source types for %s!\n'...
        , upper( mfilename ), filename )
    return
end

%----------------------------------------------------------------------%
% determine levels to be tested
%----------------------------------------------------------------------%
[ levels, stypes ]      = get_stimlevels( filebase, stimType );
if isempty( levels )
    return
end
if ~isempty( sourcetypes )
    lidx                = ismember( stypes, sourcetypes{ 1 } );
    if sum( lidx ) == 0
        return
    end
    levels              = levels( lidx, : );
end
pStats                  = [];
pInfo                   = [];

% dilute irrelevant range 
kidx                    = isoverlap( levels, valRange );
levels( ~kidx, : )      = [];
if ~isempty( pInfo )
    pStats( ~kidx )     = [];
    pInfo( ~kidx )      = [];
end

%----------------------------------------------------------------------%
% get the stats by calling DCanalysis
%----------------------------------------------------------------------%
nlev                    = size( levels, 1 );
ntot                    = NaN * ones( nlev, 2 );
nact                    = ntot;
gain                    = ntot;
lat                     = ntot;
nsup                    = ntot;
gainSup                 = ntot;
latSup                  = ntot;
nt                      = NaN * ones( nlev, 1 );
tdur                    = nt;
for j                   = 1 : nlev

    % load/generate the DCanalysis.m results:
    if ~isempty( pInfo )
        pinfo           = pInfo{ j };
        pstats          = pStats{ j };
    else
        [ pstats, pinfo ] = DCanalysis( filebase, 's', s, 'ilevel', ilevel...
            , 'wavRange', wavRange, 'sourcetypes', sourcetypes ...
            , 'stimType', stimType, 'valRange', levels( j, : ), 'durRange', durRange, 'uStim', uStim...
            , 'dxs', disRange, 'figdir', figdir, 'matdir', matdir, 'Overwrite', Overwrite, 'toplot', toplot );
    end
    
    % now get the stats:
    if isempty( pinfo )
        continue
    end
    
    didx                = pinfo( :, 6 ) == 0;                               % local stim
    pidx                = pinfo( :, 3 ) == 1;                               % PYR
    counts              = pstats( :, 1 : 2 );                               % count
    gg                  = calc_gain( pstats( :, [ 4 3 ] ), [], 2 );         % gain
    ll                  = pstats( :, 8 ) * 1000;                            % latency [ms]
    pvals               = pstats( :, 5 );                                   % Poiss (n or more)
    sigAct              = pvals <= pTH;                                     % activated by Poisson
    if size( pstats, 2 ) > 8
        pvalsSup        = pstats( :, 9 );
    else
        lambda          = pstats( :, 3 ) .* pinfo( :, 4 ) .* pinfo( :, 5 ); % rout * nt * mdur 
        nin             = pstats( :, 2 );
        pvalsSup        = poisscdf( nin, lambda );                          % Poiss( n or less )
    end
    sigSup              = pvalsSup <= pTH;                                  % suppressed by Poisson
    act                 = sigAct & gg > 1 & counts( :, 2 ) > minSpikes;     % activated
    sup                 = sigSup & gg < 1;                                  % suppressed
    
    for ct              = 0 : 1
        if ct == 0
            cidx = ~pidx;
        else
            cidx = pidx;
        end
        ntot( j, ct + 1 )       = sum( cidx & didx );
        
        nact( j, ct + 1 )       = sum( act & didx & cidx );
        gain( j, ct + 1 )       = nanmedian( gg( act & didx & cidx ) );
        lat( j, ct + 1 )        = nanmedian( ll( act & didx & cidx ) );
        
        nsup( j, ct + 1 )       = sum( sup & didx & cidx );
        gainSup( j, ct + 1 )    = nanmedian( gg( sup & didx & cidx ) );
        latSup( j, ct + 1 )     = nanmedian( ll( sup & didx & cidx ) );
    end
    nt( j, 1 )          = mean( pinfo( didx, 4 ) );
    tdur( j, 1 )        = mean( pinfo( didx, 5 ) );

end % levels

%----------------------------------------------------------------------%
% select the level
%----------------------------------------------------------------------%
% select the best intensity:
% -largest number of activated same-type cells
% -if tied, shortest latency responses
% -if tied, highest gain
% -if tied, highest intensity (PV) or lowest (CaMKII)

mlevels                 = mean( levels, 2 );
if supFlag 
    nact                = nsup;
    gain                = gainSup;
    lat                 = latSup;
    func                = @min;
    ostr                = 'suppressed';
else
    func                = @max;
    ostr                = 'activated';
end

nn                      = nact( :, celltype + 1 );
idx                     = nn == max( nn );
if sum( idx ) > 1 && sum( ~isnan( lat( idx, celltype + 1 ) ) )
    mlat                = min( lat( idx, celltype + 1 ) );
    idx                 = lat( :, celltype + 1 ) == mlat & idx;
end
if sum( idx ) > 1 && sum( ~isnan( gain( idx, celltype + 1 ) ) )
    mgain               = feval( func, gain( idx, celltype + 1 ) );
    idx                 = gain( :, celltype + 1 ) == mgain & idx;
end

if sum( idx ) > 1 
    switch celltype
        case 0
            idx         = find( idx, 1, 'last' ); % for INT manipulations, select the higher intensity
        case 1
            idx         = find( idx, 1, 'first' ); % for PYR manipulations, the lower
    end
end
if isempty( idx ) || isequal( idx, 0 ) || all( idx == 0 )
    idx                 = size( levels, 1 );
end
slevel                  = levels( idx, : );       
if uStim
    lstr                = 'local'; 
else
    lstr                = 'general'; 
end

%----------------------------------------------------------------------%
% summarize and compile output
%----------------------------------------------------------------------%
% UI
str                     = sprintf( '%s: [%0.2g-%0.2g] mA selected. %d/%d %s %s (%s %s stim); median gain %0.2g, median latency %0.2g ms'...
    , filename, slevel( 1 ) * 1000, slevel( 2 ) * 1000....
    , nn( idx ), ntot( idx, celltype + 1 ), ostr, cstr, lstr, srsstr...
    , gain( idx, celltype + 1 ), lat( idx, celltype + 1 ) );
fprintf( '%s: %s\n', upper( mfilename ), str );

% compile the output arguments (sstats and stats)
sstats                  = [ nn( idx ), gain( idx, celltype + 1 ), lat( idx, celltype + 1 ) ];
stats.filename          = repmat( { filename }, [ size( levels, 1 ) 1 ] );
stats.nt                = nt;
stats.tdur              = tdur;
stats.levels            = levels;
stats.ntot              = ntot;
stats.nact              = nact;
stats.gain              = gain;
stats.lat               = lat;

% get the gain at the selected intensity for all cells
if isa( idx, 'logical' )
    j                   = find( idx );
else
    j                   = idx;
end
if ~isempty( pInfo )
    pinfo               = pInfo{ j };
    pstats              = pStats{ j };
else
    [ pstats, pinfo ]   = DCanalysis( filebase, 's', s, 'ilevel', ilevel...
        , 'wavRange', wavRange, 'sourcetypes', sourcetypes ...
        , 'stimType', stimType, 'valRange', levels( j, : ), 'durRange', durRange, 'uStim', uStim...
        , 'dxs', disRange, 'figdir', figdir, 'matdir', matdir, 'Overwrite', Overwrite, 'toplot', toplot );
end
if isempty( pinfo )
    stats.filename2     = [];
    stats.shankclu      = [];
    stats.gains         = [];
else
    didx                = pinfo( :, 6 ) == 0;                              % local stim
    gg                  = calc_gain( pstats( :, [ 4 3 ] ), [], 2 );          % gain
    stats.filename2     = repmat( { filename }, [ sum( didx ) 1 ] );
    stats.shankclu      = pinfo( didx, 1 : 3 );
    stats.gains         = gg( didx, : );
end
    
if ~graphics
    return
end

%----------------------------------------------------------------------%
% plot
%----------------------------------------------------------------------%
% xy position - gain-latency
% dot size -    number of activated (suppressed) cells
% dot color -   intensity used 

map                     = myjet;
vv                      = mlevels( : ); 
maxVal                  = max( vv( ~isinf( vv ) ) );
pcolors                 = map( ceil( mlevels( isnan( gain( :, celltype + 1 ) ) ) / maxVal * size( map, 1 ) ), : );

figure
hold on
MS                      = nn / max( nn ) * 60;
line( gain( :, celltype + 1 ), lat( :, celltype + 1 ) );
for j                   = nlev: -1 : 1
    if isnan( gain( j, celltype + 1 ) )
        continue
    end
    ph                  = plot( gain( j, celltype + 1 ), lat( j, celltype + 1 ), '.' );
    set( ph, 'color', pcolors( j, : ), 'markersize', MS( j ) )
end
xlim( [ max( 0.1, min( [ 0.1 xlim ] ) ) max( [ 100 xlim ] ) ] )
ylim( [ max( 1, min( [ 1 ylim ] ) ) max( [ 256 ylim ] ) ] )
alines( 1, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
alines( 8, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
set( gca, 'xscale', 'log', 'yscale', 'log' )
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Latency (ms)' )
title( sprintf( 'Gain; max I:%0.3g mA', maxVal * 1000 ) )
colormap( map )
title( sprintf( '%s; %d/%d %s ', replacetok( filename, '\_', '_' ) ...
    , max( nact( :, celltype + 1 ) ), max( ntot( :, celltype + 1 ) ), cstr ) )
colorbar;

return

% EOF
