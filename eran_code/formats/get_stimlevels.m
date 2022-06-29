% get_stimlevels        determine the stimulation intensities actually used in a file/session
%
% CALL                  [ levels nvals ] = get_stimlevels( filebase, stimType )
% 
% GETS                  filebase
%                       stimType
% 
% CALLS
%                       ParseArgPairs, LoadStims, get_stimchans, stim_select
%                       optclust, uhist, sortranges, mergeranges, minmax
%                       cell2num, replacetok
%
% note:
% kmeans is way faster and generally OK, but GMM (in log space) is more or
% less perfect, although requires much more time.. however, GMM sometimes
% leads to wierd effects (e.g. background Gaussian etc) so in general
% K-means should be used. 

% 13-may-13 ES

% revisions
% 28-may-13 (1) distinct processing of distinct sources (still pool over
%               sources, i.e. all LEDs together but keep a pointer to LED/LD etc)
%           (2) defaults changes: CMODE1 (LED/DPSS) kmeans; filtFract2 (LD), 0.02
%           (3) grpahics modified
%           (4) merge low-count clusters (5 samples)
% 04-jun-13 (1) wavRange added
% 06-oct-14 (1) 'unique' mode added 
% 09-dec-14 (1) 'useMax' flag added (default; otherwise uses the mean, 
%               which is much more stable for WN in a noisy recording system..)
% 16-may-16 (1) minClustersForKMeans added (previously 4), maxClusters
% 17-aug-19 cleaned up

%[ levels nvals ] = get_stimlevels( { '25nov11', -3 }, { 'PULSE', 'PSINE' } );

% to do: 
% add an option to cluster individual channels separately and then combine
% this is useful when two channels have very different ranges but high
% resolition in each (e.g. 17mar12_1)

function [ levels, stypes, nvals, x ] = get_stimlevels( filebase, stimType, varargin )

%---------------------------------------------------------------%
% initialize output
%---------------------------------------------------------------%
levels                  = [];
stypes                  = [];
nvals                   = [];
x                       = [];

%---------------------------------------------------------------%
% constants
%---------------------------------------------------------------%
VFLAG                   = 0;

% clustering parameters
CMODE1                  = 'max';             % was 'SIC'
CMODE2                  = 'max';
FILTFRACT1              = 0.02;          % post-hoc: LED: each range is expanded by 2% of the geometric mean
FILTFRACT2              = 0.02;          % post-hoc: LD: each range is expanded by 2% of the range
MINCOUNT                = 5;               % post-hoc: minimal number of elements/cluster

% source types
sourcetypes             = { 'LED', 'LD', 'DPSS' };

%---------------------------------------------------------------%
% arguments
%-------------------------------- - ------------------------------%
mfname                  = upper( mfilename );

nargs                   = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( stimType )
    stimType            = 'ZAP'; 
end
[ Overwrite, graphics, durRange, valRange, valRanges, channels, wavRange ...
    , cmode1, cmode2, logt, filtFract1, filtFract2, minCount, vflag ...
    , useMax, minClustersForKMeans, maxClusters ] = ParseArgPairs(...
    { 'Overwrite', 'graphics', 'durRange', 'valRange', 'valRanges', 'channels', 'wavRange'...
    , 'cmode1', 'cmode2', 'logt', 'filtFract1', 'filtFract2', 'minCount', 'vflag' ...
    , 'useMax', 'minClustersForKMeans', 'maxClusters' }...
    , { -2, 0, [ 0.02 inf ], [ 0 inf ], [], [], [ 0 inf ]...
    , CMODE1, CMODE2, 1, FILTFRACT1, FILTFRACT2, MINCOUNT, VFLAG ...
    , 1, 4, [] }...
    , varargin{ : } );

% process the arguments
if isa( filebase, 'char' ) && exist( fileparts( filebase ), 'file' )
elseif isa( filebase, 'cell' ) && length( filebase ) >= 2
    datenum             = filebase( 1 : 2 );
    filebase            = datenum2filebase( datenum );
end
[ pathname, filename, extname ] = fileparts( filebase );
filename                = [ filename extname ];

if isa( stimType, 'char' )
    stimstr             = stimType;
elseif isa( stimType, 'cell' ) && isa( stimType{ 1 }, 'char' )
    stimstr             = stimType{ 1 };
else
    error( 'stimType mismatch' )
end
if isempty( channels )
    chstr               = '';
else
    if length( channels ) == 1
        chstr           = num3str( channels );
    else
        chstr           = num3str( channels )';
        chstr           = chstr( : )';
    end
    chstr               = [ '.' chstr ];
end
savename                = [ filebase '.lvl.' stimstr chstr ];
ntypes                  = length( sourcetypes );
if length( logt ) ~= ntypes
    logt                = repmat( logt( 1 ), [ ntypes 1 ] );
end
if isempty( valRanges ) || size( valRanges, 1 ) ~= ntypes
    valRanges           = repmat( valRange, [ ntypes 1 ] );
end

%---------------------------------------------------------------%
% get the levels
%---------------------------------------------------------------%
reDo                    = 1;
if Overwrite < 0 && exist( savename, 'file' )
    
    % load existing levels
    verb( sprintf( '%s: loading %s...', mfname, savename ), vflag )
    L                   = load( savename, '-mat' );
    if ~isfield( L, 'stypes' )
        reDo            = 1;
        Overwrite       = abs( Overwrite );
    else
        reDo            = 0;
        load( savename, '-mat' )
        x               = [];
    end
    
end

if reDo

    % cluster de-novo
    verb( sprintf( '%s: Clustering %s (%s):', mfname, filebase, stimstr ), vflag )
    
    % get the actually-used intensities in that session
    [ ~, ~, stims ]   = LoadStims( filebase );
    if isempty( channels )
        channels        = unique( [ stims.chan ] );
        % take only the relevant wavelength channels
        if ~isequal( wavRange, [ 0 inf ] )
            [ stimchans, stimtargets, stimvoltageranges, stimsources, stimwavelengths ] = get_stimchans( filebase, channels );
            widx                        = inrange( stimwavelengths, wavRange );
            stimchans( ~widx )          = [];
            stimtargets( ~widx )        = [];
            stimvoltageranges( ~widx )  = [];
            stimsources( ~widx )        = [];
            stimwavelengths( ~widx )    = [];
            channels                    = stimchans;
            kidx                        = false( length( stims ), 1 );
            for i                       = 1 : length( stims )
                if all( ismember( stims( i ).chan, channels ) )
                    kidx( i )           = 1;
                end
            end
            stims( ~kidx )              = [];
        end
    end
    
    x                                   = cell( ntypes, 1 );
    for j                               = 1 : length( stims )
        if ~all( ismember( stims( j ).chan, channels ) )
            continue
        end
        stimZ                           = stim_select( stims( j ), 'types', stimType, 'durs', durRange, 'vals', valRange );
        if useMax
            if isequal( stimZ.source, sourcetypes{ 1 } )
                x{ 1 }                  = [ x{ 1 }; stimZ.vals ];
            elseif isequal( stimZ.source, sourcetypes{ 2 } )
                x{ 2 }                  = [ x{ 2 }; stimZ.vals ];
            elseif isequal( stimZ.source, sourcetypes{ 3 } ) || isequal( stimZ.source, 'SHUT' )
                x{ 3 }                  = [ x{ 3 }; stimZ.vals ];
            end
        else
            if isequal( stimZ.source, sourcetypes{ 1 } )
                x{ 1 }                  = [ x{ 1 }; stimZ.stats( :, 1 ) ];
            elseif isequal( stimZ.source, sourcetypes{ 2 } )
                x{ 2 }                  = [ x{ 2 }; stimZ.stats( :, 1 ) ];
            elseif isequal( stimZ.source, sourcetypes{ 3 } ) || isequal( stimZ.source, 'SHUT' )
                x{ 3 }                  = [ x{ 3 }; stimZ.stats( :, 1 ) ];
            end
        end
    end
    
    % cluster them using GMM/kmeans
    nvals                               = [];
    levels                              = [];
    stypes                              = [];
    nx                                  = zeros( ntypes, 1 );
    for i                               = 1 : ntypes
        
        if isempty( x{ i } )
            continue
        end
        x{ i }( ~inrange( x{ i }, valRanges( i, : ) ) ) = [];
        nx( i )                         = length( x{ i } );
        if nx( i ) == 0
            continue
        end
        if i == 1 || i == 3
            cmode                       = cmode1; % LED, DPSS: GMM
        else
            cmode                       = cmode2; % LD: kmeans
        end
        
        % cluster
        fprintf( '%s: Clustering %d values (type%d)... ', mfname, length( x{ i } ), i )
        if logt( i )
            xi                          = x{ i };
            xi( xi <= 0 )               = eps;
            xc                          = log10( xi );
        else
            xc                          = x{ i };
        end
        if isequal( cmode, 'unique' )
            [ uvals, ~, xclu ]          = unique( xc );
        elseif isempty( maxClusters )
            xclu                        = optclust( xc, cmode, 0, [], [ 20 100 ], 0 );
        else
            xclu                        = optclust( xc, cmode, 0, [ 2 maxClusters ], [ 20 100 ], 0 );
        end
        uxclu                           = unique( xclu );
        nxclu                           = length( uxclu );
        if isequal( cmode, 'max' ) && nxclu <= minClustersForKMeans
            xclu                        = optclust( xc, 'SIC', 0, [], [], 0 );
            uxclu                       = unique( xclu );
            nxclu                       = length( uxclu );
        end
        
        % determine the indpendent levels for each cluster
        fprintf( '%d clusters; post-processing... ', nxclu )
        newlevels                       = zeros( nxclu, 2 );
        for j = 1 : nxclu
            newlevels( j, : )           = minmax( x{ i }( xclu == uxclu( j ) ) );
        end
        newvals                         = uhist( xclu )';
        [ newlevels, sidx ]            	= sortrows( newlevels, 1 );
        newvals                         = newvals( sidx );
        [ newlevels, idxy, idxx ]       = sortranges( newlevels, 1 );
        newvals = accumarray( idxx, newvals );
        
        % merge clusters with few samples to prevent overclustering
        [ newlevels, newvals ]          = mergeranges( newlevels, newvals, minCount );
        
        % expand ranges to prevent overclustering
        if i == 1 || i == 3 % LED, DPSS
            [ nlevels, idxy, idxx ]     = sortranges( newlevels + geomean( newlevels, 2 ) * filtFract1 * [ -1 1 ] );
            nnvals                      = accumarray( idxx, newvals );
        else % LD
            [ nlevels, idxy, idxx ]     = sortranges( newlevels + diff( newlevels, 1, 2 ) * filtFract2 * [ -1 1 ] );
            nnvals                      = accumarray( idxx, newvals );
        end
        
        % regularize the levels
        dx                              = diff( sort( [ 0; nlevels( : ) ] ) );
        dx                              = min( dx( dx > 0 ) );
        resolution                      = 10.^( -ceil( log10( 1 / dx ) ) - 1 );
        nlevels                         = [ remrnd( nlevels( :, 1 ), resolution, 'floor' ) ...
            remrnd( nlevels( :, 2 ), resolution, 'ceil' ) ];
        [ nlevels, idxy, idxx ]         = sortranges( nlevels, 1 );
        nnvals                          = accumarray( idxx, nnvals );
        if isequal( cmode, 'unique' )
            nlevels                     = bsxfun( @plus, nlevels, eps * [ -1 1 ] );
        end
        
        % accumulate over sources, keeping a pointer to the source type
        fprintf( 'done: %d clusters.\n', size( nlevels, 1 ) )
        levels                          = [ levels; nlevels ];
        nvals                           = [ nvals; nnvals ];
        stypes                          = [ stypes; repmat( sourcetypes( i ), [ size( nlevels, 1 ) 1 ] ) ];
        
    end

    % determine the independent levels for the filtered levels
    if isempty( levels )
        verb( sprintf( '%s: No data for %s in this file', mfname, stimstr ), vflag )
        return
    end
       
    % save 
    if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( savename, 'file' ) )
        %verb( sprintf( '%s: Done clustering. Saving %s', mfname, savename ), vflag )
        fprintf( '%s: Done clustering. Saving %s\n', mfname, savename )
        save( savename, 'levels', 'nvals', 'stypes', 'valRange', 'stimType', 'durRange', 'channels', '-v6' )
    end

end

%---------------------------------------------------------------%
% post-hoc fixes
%---------------------------------------------------------------%
tidx                        = [];
cidx                        = [];

if ~isempty( cidx ) && sum( cidx ) > 1

    if isempty( tidx )
        clevels             = minmax( levels( cidx, : ) );
        cnvals              = sum( nvals( cidx ) );
        levels              = [ levels( ~cidx, : ); clevels ];
        nvals               = [ nvals(  ~cidx ); cnvals ];
        [ levels, sidx ]    = sortrows( levels, 1 );
        nvals               = nvals( sidx );
        if exist( 'stypes', 'var' )
            cidx( find( cidx, 1, 'first' ) ) = 0; 
            stypes( cidx )  = [];
        end
    else
        % keep those aside
        klevels             = levels( ~tidx, : );
        knvals              = nvals( ~tidx, : );
        kstypes             = stypes( ~tidx, : );
        % process the relevant ones
        tlevels             = levels( tidx, : );
        tnvals              = nvals( tidx, : );
        tstypes             = stypes( tidx, : );

        clevels             = minmax( tlevels( cidx, : ) );
        cnvals              = sum( tnvals( cidx ) );
        tlevels             = [ tlevels( ~cidx, : ); clevels ];
        tnvals              =  [ tnvals(  ~cidx ); cnvals ];
        [ tlevels, sidx ]   = sortrows( tlevels, 1 );
        tnvals              = tnvals( sidx );
        tstypes( 1 : ( sum( cidx ) - 1 ) ) = [];
        % combine again
        levels              = [ klevels; tlevels ];
        nvals               = [ knvals; tnvals ];
        stypes              = [ kstypes; tstypes ];
    end
    str                     = sprintf( '%s: Post-fixed levels: ', mfname );
else
    str                     = sprintf( '%s: Raw levels: ', mfname );
end

%---------------------------------------------------------------%
% temporary fix:
%---------------------------------------------------------------%
if ( Overwrite == -3 && ~exist( 'stypes', 'var' ) ) || ( graphics && isempty( x ) )
    
    % get the actually-used intensities in that session
    [ ~, ~, stims ]         = LoadStims( filebase );
    if isempty( channels )
        channels            = unique( [ stims.chan ] );
        % take only the relevant wavelength channels
        if ~isequal( wavRange, [ 0 inf ] )
            [ stimchans, stimtargets, stimvoltageranges, stimsources, stimwavelengths ] = get_stimchans( filebase, channels );
            widx                            = inrange( stimwavelengths, wavRange );
            stimchans( ~widx )              = [];
            stimtargets( ~widx )            = [];
            stimvoltageranges( ~widx )      = [];
            stimsources( ~widx )            = [];
            stimwavelengths( ~widx )        = [];
            channels                        = stimchans;
            kidx                            = false( length( stims ), 1 );
            for i                           = 1 : length( stims )
                if all( ismember( stims( i ).chan, channels ) )
                    kidx( i )   = 1;
                end
            end
            stims( ~kidx )      = [];
        end
    end
    x                           = cell( ntypes, 1 );
    for j                       = 1 : length( stims )
        if ~all( ismember( stims( j ).chan, channels ) )
            continue
        end
        stimZ                   = stim_select( stims( j ), 'types', stimType, 'durs', durRange, 'vals', valRange );
        if isequal( stimZ.source, sourcetypes{ 1 } )
            x{ 1 }              = [ x{ 1 }; stimZ.vals ];
        elseif isequal( stimZ.source, sourcetypes{ 2 } )
            x{ 2 }              = [ x{ 2 }; stimZ.vals ];
        elseif isequal( stimZ.source, sourcetypes{ 3 } ) || isequal( stimZ.source, 'SHUT' )
            x{ 3 }              = [ x{ 3 }; stimZ.vals ];
        end
    end
    nx                          = zeros( ntypes, 1 );
    for i                       = 1 : ntypes
        if ~isempty( x{ i } )
            if exist( 'stypes', 'var' )
                arange          = minmax( levels( ismember( stypes, sourcetypes{ i } ), : ) );
            else
                arange          = minmax( levels );
            end
            if isempty( arange )
                continue
            end
            x{ i }( ~inrange( x{ i }, arange ) ) = [];
        end
        nx( i )                 = length( x{ i } );
    end
end

if Overwrite == -3 && ~exist( 'stypes', 'var' )
    
    if sum( nx > 0 ) == 1
        stypes = repmat( sourcetypes( nx > 0 ), [ size( levels, 1 ) 1 ] );
        fprintf( '%s: Re-saving %s\n', mfname, savename ) 
        save( savename, 'levels', 'nvals', 'stypes', 'valRange', 'stimType', 'durRange', 'channels', '-v6' )
    else
        fprintf( '%s: multiple sources for %s, cannot resave - should rerun!!!\n', mfname, filename )
        return
    end
    
end

%---------------------------------------------------------------%
% summary
%---------------------------------------------------------------%    
if exist( 'stypes', 'var' )
    svec                        = cell2num( stypes, sourcetypes );
else
    svec                        = [];
end
verb( str, vflag )
if vflag 
    disp( [ levels * 1000  nvals svec ] )
end

if graphics
    if isempty( x )
        for i = 1 : 3
            nx( i )             = sum( nvals( ismember( stypes, sourcetypes{ i } ) ) );
        end
    end
    figure;
    j = 0;
    for i = 1 : 3
        if nx( i ) == 0
            continue
        end
        j = j + 1;
        subplot( sum( nx > 0 ), 1, j )
        if isempty( svec )
            tidx                = true( size( nvals ) );
        else
            tidx                = ismember( stypes, sourcetypes{ i } );
        end
        tlevels                 = levels( tidx, : );
        bords                   = 1000 * mean( [ tlevels( 1 : end - 1, 2 ) tlevels( 2 : end, 1 ) ], 2 );
        if ~isempty( x )
            xc                  = 1000 * x{ i };
            if logt( i )
                bords           = log2( bords );
                xc              = log2( xc );
            end
            nbins               = min( max( 100, floor( nx( i ) / 10 ) ), 500 );
            hist( xc, nbins )
        end
        alines( bords, 'x', 'color', [ 1 0 0 ] );
        set( gca, 'tickdir', 'out', 'box', 'off' )
        if logt( i )
            lin2log( 'x', 2, 0.25 );
            xlabel( 'mV (log scale)' )
        end
        xlabel( 'mV' )
        ylabel( 'Count' )
        if j == 1
            tstr                = [ repmat( 'C', size( channels( : ) ) ) num2str( channels( : ) ) repmat( ' ', size( channels( : ) ) ) ]';
            tstr                = tstr( : ).';
            title( sprintf( '%s, %s: %s; %s', replacetok( filename, '\_', '_' )...
                , stimstr, tstr, sourcetypes{ i } ) )
        else
            title( sourcetypes{ i } )
        end
    end
    
end

return

% EOF
