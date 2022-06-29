% seqrast
%
% call                  [ mat, bins, rh, ah ] = seqrast( clu, res, periods, varargin )


function [ mat, tim, rh, ah ] = seqrast( clu, res, periods, varargin )

% initialize output
rh                                      = [];
ah                                      = [];

% default values
binsizeSEC_DEFAULT                      = 0.001;            % [s]
padSEC_DEFAULT                          = [ -0.05 0.05 ];   % [s]
spkFs_DEFAULT                           = 20000;            % [Hz]

% argument handling
nargs                           = nargin;
if nargs < 3 || isempty( clu ) || isempty( res ) || isempty( periods )
    return
end
[ shankclu, map ...
    , binsizeSEC, padSEC, spkFs ...
    , stimes, graphics ...
    , fAlpha, eAlpha, pColor ]  = ParseArgPairs(...
    { 'shankclu', 'map' ...
    , 'binsizeSEC', 'padSEC', 'spkFs' ...
    , 'stimes', 'graphics' ...
    , 'fAlpha', 'eAlpha', 'pColor' }...
    , { [], [] ...
    , binsizeSEC_DEFAULT, padSEC_DEFAULT, spkFs_DEFAULT ...
    , [], 1 ...
    , 0.1, 0.1, [ 0 0 0.7 ] } ...
    , varargin{ : } );

% (1) get the rasters as sparse martices for all units in clu
binsize                         = ceil( binsizeSEC * spkFs );                       % [samples]
h                               = [ 0 max( diff( periods, [], 2 ) + 1 ) - 1 ];      % [samples]
h1                              = [ floor( h( 1 ) / binsize ) ceil( h( 2 ) / binsize ) ]; % [bins]
h2                              = [ floor( padSEC( 1 ) / binsizeSEC ) ceil( padSEC( 2 ) / binsizeSEC ) ]; % [bins]
halfwin                         = h1 + h2;
[ r, bins ]                     = get_rasters( clu, res, [], periods( :, 1 )...
    , 'binsize', binsize, 'halfwin', halfwin );
tim                             = bins / spkFs;

% (2) now, keep only the desired units and organize in the requested order
if ~isempty( shankclu )
    uidx                            = ismember( map( :, 2 : 3 ), shankclu( :, 1 : 2 ), 'rows' );
    if isempty( uidx )
        error( '' )
    end
    r                               = r( uidx );
    
    [ ~, bb ]                       = sortrows( shankclu( :, 1 : 2 ) );
    [ ~, ridx ]                     = sort( bb );
    r                               = r( ridx );
end

% (3) combine into one array
ntrials                         = size( periods, 1 );
nbins                           = size( r{ 1 }, 1 );
nunits                          = length( r );
if ~isequal( size( r{ 1 }, 2 ), ntrials )
    error( 'mismatch' )
end
if ~isequal( size( r{ 1 }, 1 ), nbins )
    error( 'mismatch' )
end
mat                             = sparse( nbins, nunits * ntrials );
for i                           = 1 : nunits
    cidx                        = ntrials * ( i - 1 ) + 1 : ntrials * i;
    mat( :, cidx )              = r{ i };
end

% (4) plot
if ~graphics
    return
end
newplot
% plot the raster
rh                              = plot_raster( mat, tim );
ah                              = alines( ntrials * ( 1 : ( nunits - 1 ) ) + 0.5, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel( 'Time [s]' )
ylabel( 'Trial number x unit' )

% add patches
if ~isempty( stimes )
    stimes                      = permute( stimes( 1, :, : ), [ 2 3 1 ] );
    nchans                      = size( stimes, 1 );
    for i                       = 1 : nchans
        xe                    	= stimes( i, : );
        ye                    	= ylim;
        ph                      = patch( xe( [ 1 2 2 1 1 ] ), ye( [ 1 1 2 2 1 ] ), pColor );
        set( ph, 'FaceAlpha', fAlpha, 'EdgeAlpha', eAlpha )
    end
end

% add ustr
if ~isempty( shankclu )
    ye                          = ntrials / 2 : ntrials : ntrials * ( nunits - 0.5 );
    unums                       = shankclu( :, 1 : 2 );
    ustr                        = cell( 1, nunits );
    for i                       = 1 : nunits
        ustr{ i }               = sprintf( '%d.%d', unums( i, 1 ), unums( i, 2 ) );
    end
    set( gca, 'ytick', ye, 'YTickLabel', ustr )
end

return

% EOF
