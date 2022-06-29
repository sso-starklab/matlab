% shirly_dtdt           2D place fields on DTDT task

% 21-mar-21 ES + SSo

function fig = shirly_dtdt( filebase )

% session specific parameters (specific to mS234_20)
tlim                = [ 1700 5370 ];                                        % temporal limits in s
xlims               = [ 70 inf ];
ylims               = [ 40 inf ];

% session specific parameters (specific to mS234_19)
% tlim                = [ 2450 5260 ];                                        % temporal limits in s
% xlims               = [ 80 inf ];
% ylims               = [ 50 130 ];
% 
% % session specific parameters (specific to mS234_04)
% tlim                = [ 12000 18000 ];                                        % temporal limits in s
% xlims               = [ 70 inf ];
% ylims               = [ 50 inf ];
% 
% % session specific parameters (specific to mS234_16)
% tlim                = [ 1 4300 ];                                        % temporal limits in s
% xlims               = [ 70 inf ];
% ylims               = [ 50 inf ];
% 
% 
% % session specific parameters (specific to mS234_18)
% tlim                = [ 6300 8966 ];                                        % temporal limits in s
% xlims               = [ 70 inf ];
% ylims               = [ 50 inf ];
% 
% % session specific parameters (specific to mS234_02)
% tlim                = [ 10000 14430 ];                                        % temporal limits in s
% xlims               = [ 70 inf ];
% ylims               = [ 50 inf ];
% 
% % session specific parameters (specific to mS234_17)
% tlim                = [ 1500 4500 ];                                        % temporal limits in s
% xlims               = [ 70 inf ];
% ylims               = [ 50 inf ];
% % unit selection constants
ilevel                          = 'B';

% PF calculation constants
minspd                          = 5; % [cm/s]
binSize                         = 2.5; % [cm]
smoothSD                        = 5; % [cm]

%----------------------------------------------------------------------
% (1) spatial information
load( [ filebase '.mov.mat' ], '-mat', 'mov' )
x                       = mov.pos( :, 1 );
y                       = mov.pos( :, 2 );
spd                     = mov.spd;
nx                      = length( x );
movFs                   = mov.Fs;
dt                      = 1 / movFs;
t                       = dt : dt : dt * nx;

% plot the position
fig( 1 )                = figure;
subplot( 2, 2, 1 )
plot( t, x, '-b' )

subplot( 2, 2, 2 )
plot( x, y, '.b', 'MarkerSize', 1 )

% temporal range
idx                 = ceil( tlim( 1 ) * mov.Fs ) : floor( tlim( 2 ) * mov.Fs );
tidx                = false( nx, 1 );
tidx( idx )         = 1;

subplot( 2, 2, 3 )
plot( x( tidx ), y( tidx ), '.b', 'MarkerSize', 1 )

% spatial range
xidx                = x >= xlims( 1 ) & x <= xlims( 2 );
yidx                = y >= ylims( 1 ) & y <= ylims( 2 );
idx                 = xidx & yidx & tidx;

subplot( 2, 2, 4 )
plot( x( idx ), y( idx ), '.b', 'MarkerSize', 1 )

if isinf( xlims( 2 ) )
    xlims               = [ xlims( 1 ) max( x( idx ) ) ];
end
if isinf( ylims( 2 ) )
    ylims               = [ ylims( 1 ) max( y( idx ) ) ];
end

%----------------------------------------------------------------------
% (2) spike times
spk                 = load_spikes( filebase, [], ilevel );
par                 = LoadXml( filebase );
spkFs               = par.SampleRate;

%----------------------------------------------------------------------
% (3) go over units and compute PF

% prepare edges for grids
edgesX                          = xlims( 1 ) - binSize / 2 : binSize : xlims( 2 ) + binSize / 2;
edgesY                          = ylims( 1 ) - binSize / 2 : binSize : ylims( 2 ) + binSize / 2;
binEdges                        = { edgesX, edgesY };

periods                         = tlim;

% prepare figure for plots
nunits                          = size( spk.shankclu, 1 );
nrows                           = floor( sqrt( nunits ) );
ncols                           = ceil( nunits / nrows );
[ ah, fig( 2 ) ]                = tilefig( nrows, ncols );
colormap( myjet ) 

for i                           = 1 : nunits
    fprintf( 1, 'Computing PF for %d.%d...\n', spk.shankclu( i, 1 ), spk.shankclu( i, 2 ) )
    clunum                      = spk.map( i, 1 );
    st                          = spk.res( spk.clu == clunum );
    
%     [ F, X ]                    = calc_field( st / spkFs, x, 5, movFs );
%     [ F, X, C, T ]              = calc_field( st / spkFs, x, binSize, movFs, periods, minspd, [], [], spd, smoothSD );
%    [ F, X, C, T ]              = calc_field( st / spkFs, [ x y ], binSize, movFs, periods, minspd, [], [], spd, smoothSD );
    [ F, X, C, T ]              = calc_field( st / spkFs, [ x y ], binEdges, movFs, periods, minspd, [], [], spd, smoothSD );
    
    %imagesc( X{ 1 }, X{ 2 }, T' ), colormap( myjet ), axis xy
    subplot( ah( i ) )
    ph                          = imagesc( X{ 1 }, X{ 2 }, F', 'CDataMapping','scaled' );
    set( ph, 'AlphaData', ~isnan( F' ) )
    axis xy, axis equal, axis tight
    set( gca, 'tickdir', 'out', 'box', 'off' )
    title( sprintf( '%d.%d (%d), %0.3g', spk.shankclu( i, 1 ) ...
        , spk.shankclu( i, 2 ), spk.shankclu( i, 3 ), max( F( : ) ) ) )
    axis off
    
end

for i                           = ( nunits + 1 ) : length( ah )
    subplot( ah( i ) )
    axis off
end

return 

% EOF

filebase = filebaseLookup( 'mS234', -20 );
fig = shirly_dtdt( filebase );

