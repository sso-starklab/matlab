% eeg2whl           convert eeg channels to movement channels
% 
% mov = eeg2whl( filebase )
%
% steps:
% (1) load eeg/dat file         default - eeg
% (2) determine channels        default - by par file
% (3) determine sampling rate   default - Fs/2^5 (39.062)
% (4) downsample the data       
% (5) identify jumps            default - 6 SDs of speed
% (6) interpolate linearly
% (7) smooth                    4 Hz
% 
% calls:
%   utility:            ParseArgPairs, LoadXml
%   downsampling:       lin_downsample
%   jumps:              derive, parse, uniteranges, enumerate, makefir, makegaussfir, linint
%   move direction:     car2pol
%   filtering:          circfilt, firfilt

% 25-feb-18 ES

% revisions
% 26-feb-18 downsampling
% 27-feb-18 jumps and interpolating; done!
% 30-aug-18 (1) modified to support only xy (without theta)
%           (2) modified parameters, e.g. pix2cm should be 0.45, not 0.357
%           (3) modified plotting - to work w/o theta; to include speed; to
%                   enable viewing all on one axis
% 31-aug-18 (1) also smooth speed, direction, and orientation (ang)
%           (2) modified back to 3.57 since otherwise the track scales to 200 cm..
% 02-sep-18 (1) modified Overwrite to default to -2
% 11-sep-18 (1) modified algorithm to support any combination of x/y/theta
%           (2) changed scaling frmo 3.57 to 4.27 (empircally to inter-sol distance of 130, was 108.5)
% 16-jul-19 (1) added NaN field for ang in case no theta signal
% 04-aug-19 (1) added check for 
% 12-sep-19 (1) figfname modified to support saving in PC
% 14-oct-19 (1) nvalid == 0 check added

% issues:
% (1) temporal discontinuities (jumps)              handled - done
% (2) block-wise analysis                           not necessary for eeg
% (3) if dat, concatenation over multiple files     not necessary for eeg
%
% (4) values not conforming to the expected. 
% should run blockwise over dat files

function mov = eeg2whl( filebase, varargin )

%--------------------------------------------------------------------%
% initialize output
figs                        = [];
mov                         = [];

% constants
deg2rad                     = pi / 180;             % spotter gives theta in degrees
DACs( 1, : )                = [ -3.8e-3 4.64 ];     % [V] corresponding to 0 and 639
DACs( 2, : )                = [ -3.8e-3 4.64 ];     % 
DACs( 3, : )                = [ -3.8e-3 4.64 ];     % 
maxBins                     = 640;

% parse arguments
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ Overwrite, graphics, preFixed, verbose ...
    , suffix, pad ...
    , LP, filtMode ...
    , pix2cm, xTH ]         = ParseArgPairs(...
    { 'Overwrite', 'graphics', 'preFixed', 'verbose' ...
    , 'suffix', 'pad' ...
    , 'LP', 'filtMode' ...
    , 'pix2cm', 'xTH' }...
    , { -2, 1, 0, 1 ...
    , 'eeg', [ -3 3 ] ...
    , 4, 'gauss' ...
    , 4.27 / 10, [ 1 6 ] ...
    } ...
    , varargin{ : } );
switch suffix
    case 'dat'
        dsf                 = 2^9;
    case 'eeg'
        dsf                 = 2^5;
end

%--------------------------------------------------------------------%
% get parameters
%--------------------------------------------------------------------%
mfname                      = upper( mfilename );
par                         = LoadXml( filebase );
nBits                       = par.nBits;
switch suffix
    case 'dat'
        rawFs               = par.SampleRate;
    case 'eeg'
        rawFs               = par.lfpSampleRate;
end
nchans                      = par.nChannels;
[ xchan, ~, xRange ]        = get_stimchans( par, [], 'x' );
[ ychan, ~, yRange ]        = get_stimchans( par, [], 'y' );
[ tchan, ~, tRange ]        = get_stimchans( par, [], 'theta' );
channels                    = NaN * ones( 1, 3 );
ranges                      = NaN * ones( 1, 3 );
if ~isempty( xchan )
    channels( 1 )           = xchan;
    ranges( 1 )             = xRange;
end
if ~isempty( ychan )
    channels( 2 )           = ychan;
    ranges( 2 )             = yRange;
end
if ~isempty( tchan )
    channels( 3 )           = tchan;
    ranges( 3 )             = tRange;
end
vchans                      = ~isnan( channels );
nvalid                      = sum( vchans );
if nvalid == 0
    return
end

% file handling
sourcefile                  = sprintf( '%s.%s', filebase, suffix );
whlfname                    = sprintf( '%s.whl', filebase );
movfname                    = sprintf( '%s.mov.mat', filebase );

[ ~, fname, extname ]       = fileparts( whlfname );
cpathi                      = strfind( filebase, 'dat' );
cpath                       = filebase( 1 : ( cpathi - 1 ) );
fpath                       = sprintf( '%sfigs/', cpath );
if ~exist( fpath, 'dir' )
    eval( sprintf( '!mkdir %s', fpath ) )
end
figfname                    = sprintf( '%s%s.mov', fpath, fname );
if ispc
    figfname                = replacetok( figfname, '\', '/' );
end

% load the mov file
if Overwrite < 0 && exist( movfname, 'file' )
    if verbose
        fprintf( '%s: Loading *.mov.mat file %s\n', mfname, movfname )
    end
    load( movfname, 'mov' )
    return
end

% load data
if verbose
    fprintf( 1, '%s: loading *%s data...\n', sourcefile, suffix );
end
a                           = memmapfile( sourcefile, 'Format', 'int16' );
n                           = length( a.data );
idx                         = ( channels( vchans )' * ones( 1, n / nchans ) + ones( nvalid, 1 ) * ( 0 : nchans : ( n - nchans ) ) )';
x                           = single( a.data( idx ) );
clear a

%--------------------------------------------------------------------%
% convert x/y/theta from eeg/dat to whl and mov
%--------------------------------------------------------------------%

% convert d2au to physical values (cm, radians)
if verbose
    fprintf( 1, 'converting and downsampling...\t' );
end
xhat                        = single( ones( size( x ) ) );
xchans                      = zeros( size( vchans ) );
tmp                         = 1 : sum( vchans );
xchans( vchans )            = tmp;
vDACs                       = DACs( vchans, : );
vranges                     = ranges( vchans );
for i                       = 1 : size( x, 2 )
    xhat( :, i )            = ( x( :, i ) / 2^nBits * vranges( i ) - vDACs( i, 1 )  )  / ( vDACs( i, 2 ) - vDACs( i, 1 ) ) * maxBins;
end
if vchans( 1 )
    xhat( :, xchans( 1 ) )  = xhat( :, xchans( 1 ) ) * pix2cm;                % pixels to cm
end
if vchans( 2 )
    xhat( :, xchans( 2 ) ) 	= xhat( :, xchans( 2 ) ) * pix2cm;                % pixels to cm
end
if vchans( 3 )
    xhat( :, xchans( 3 ) )  = xhat( :, xchans( 3 ) ) * deg2rad;               % angles to radians
end

% downsample the data
islin                       = logical( [ 1 1 0 ] );
xychans                     = islin( vchans );
tmp                         = lin_downsample( xhat( :, xychans ), dsf, true( 1, sum( xychans ) ) );
xy                          = single( zeros( size( tmp, 1 ), 2 ) );
xy( :, xychans )            = tmp;
Fs                          = rawFs / dsf;
n                           = size( xy, 1 );
if vchans( 3 )
    tmp                     = lin_downsample( xhat( :, xchans( 3 ) ), dsf, false( 1 ) );
else
    tmp                     = zeros( n, 1 );
end
whl                         = [ xy tmp ];

% identify segments of jumps
if verbose
    fprintf( 1, 'interpolating over discontinuities...\t' );
end
vel                         = derive( whl( :, 1 : 2 ), 1 / Fs );
spd                         = sqrt( sum( vel.^2, 2 ) );
mm                          = nanmean( spd );
ss                          = nanstd( spd );
idx                         = spd > ( mm + xTH( 2 ) * ss );
mat                         = parse( find( idx ) );
if ~isempty( mat )
    mat                     = bsxfun( @plus, mat, pad );
    if mat( 1 ) <= 1
        mat( 1, : )         = [];
    end
    if mat( end ) >= n
        mat( end, : )       = [];
    end
    mat                     = uniteranges( mat );
end
idx2                        = enumerate( mat );
whl0                        = whl;

% interpolate
switch filtMode
    case 'fir'
        mwin                = makefir( [ 0 LP ], Fs, [], 'low' );
    case 'gauss'
        mwin                = makegaussfir( 1/LP/3, Fs );
    case 'ma'
        MA                  = ceil( Fs / LP );
        mwin                = ones( MA, 1 ) / MA;
end
mwin                        = mwin( : ) / sum( mwin );
whl( idx2, : )              = -1000;
whl                         = linint( whl, -1000 );
whl                         = firfilt( whl, mwin );
mov.pos                     = whl( :, 1 : 2 );

% compute velocity
vel                         = derive( mov.pos, 1/Fs );
[ spd, dir ]                = car2pol( vel( :, 1 ), vel( :, 2 ) );
mov.spd                     = firfilt( spd, mwin );
mov.dir                     = mod( circfilt( dir, mwin ), 2 * pi );

% compute orientation
if vchans( 3 )
    mov.ang                 = mod( circfilt( whl( :, 3 ), mwin ), 2 * pi );
else
    mov.ang                 = NaN * ones( size( spd ) );
end

% fill the rest of the fields
mov.whl                 = [ fname extname ];
mov.Fs                  = Fs;

if verbose
    fprintf( 1, 'done!\n' );
end

%--------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%
if graphics
    
    figs                    = zeros( 1, 3 );
    fnamestr                = replacetok( [ fname extname ], '\_', '_' );
    t                       = ( 1 : size( whl, 1 ) ) / Fs;
    
    yoffset                 = 400;                  % [cm]
    soffset                 = -1200;                % cm/s
    sscale                  = 0.1;                   
   
    figs( 1 )               = figure;
    
    sp = axes( 'position', [ 0.13 0.4 0.775 0.5 ] );
    subplot( sp )
    plot( t, mov.pos( :, 1 ), '-b' )
    hold on
    plot( t, mov.pos( :, 2 ) + yoffset, '-b' )
    plot( t, mov.spd( :, 1 ) + soffset, '-b' )
    
    title( fnamestr )
    xlim( [ 0 n ] / Fs )
    alines( [ soffset 0 yoffset ], 'y', 'color', [ 1 1 1 ] * 0, 'linestyle', '--' );
    ylim( ylim + [ -1 1 ] * 100 )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    ylabel( sprintf( 'Speed%d [cm/s]; X, Y+%d [cm]', soffset, yoffset ) )
    xlabel( 'Time [s]' )
    
    sTH                     = 1;
    nbins                   = 200;
    sidx                    = mov.spd > sTH;
    
    if sum( sidx ) 
        subplot( 3, 3, 7 )
        vals                    = mov.spd( sidx );
        bins                    = logspace( log10( min( vals ) ), log10( max( vals ) ), nbins );
        hist( vals, bins )
        set( gca, 'xscale', 'log' )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        xlabel( 'Speed [cm/s]' )
        title( sprintf( '%0.3g%% > %0.3g cm/s', sum( sidx ) / length( mov.spd ) * 100, sTH ) )
        
        subplot( 3, 3, 8 )
        hist( mov.dir( sidx ), nbins )
        xlabel( 'Movement direction [rad]' )
        xlim( [ 0 2 * pi ] )
        set( gca, 'tickdir', 'out', 'box', 'off' )
    end
    
    if ~all( isnan( mov.ang ) )
        subplot( 3, 3, 9 )
        hist( mov.ang( sidx ), nbins )
        xlabel( 'Orientation [rad]' )
        xlim( [ 0 2 * pi ] )
        set( gca, 'tickdir', 'out', 'box', 'off' )
    end
    
    if preFixed
        
        figs( 2 )           = figure;
        
        subplot( 1, 1, 1 )
        plot( t, mov.pos( :, 1 ), '-b', t( idx2 ), mov.pos( idx2, 1 ), '.r' )
        hold on
        plot( t, mov.pos( :, 2 ) + yoffset, '-b', t( idx2 ), mov.pos( idx2, 2 ) + yoffset, '.r' )
        plot( t, mov.spd( :, 1 ) + soffset, '-b', t( idx2 ), mov.spd( idx2, 1 )  + soffset, '.r' )
        
        title( fnamestr )
        xlim( [ 0 n ] / Fs )
        alines( [ soffset 0 yoffset ], 'y', 'color', [ 1 1 1 ] * 0, 'linestyle', '--' );
        ylim( ylim + [ -1 1 ] * 100 )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        ylabel( sprintf( 'Speed+%d; X; Y+%d [cm]', soffset, yoffset ) )
        
        vel0                = derive( whl0, 1/Fs );
        spd0                = car2pol( vel0( :, 1 ), vel0( :, 2 ) );
        
        figs( 3 )           = figure;
        
        plot( t, whl0( :, 1 ), '-b', t( idx2 ), whl0( idx2, 1 ), '.r' )
        hold on
        plot( t, whl0( :, 2 ) + yoffset, '-b', t( idx2 ), whl0( idx2, 2 ) + yoffset, '.r' )
        plot( t, ( spd0( :, 1 ) * 0.1 + soffset ), '-b', t( idx2 ), ( spd0( idx2, 1 ) * 0.1 + soffset ), '.r' )
        
        title( fnamestr )
        xlim( [ 0 n ] / Fs )
        alines( [ soffset 0 yoffset ], 'y', 'color', [ 1 1 1 ] * 0, 'linestyle', '--' );
        ylim( ylim + [ -1 1 ] * 100 )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        ylabel( sprintf( '%0.3g*Speed+%d; X; Y+%d [cm]', sscale, soffset, yoffset ) )
        
    end
    
end

%--------------------------------------------------------------------%
% save the whl and mov files
%--------------------------------------------------------------------%
if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( whlfname, 'file' ) )
    
    % save the whl
    fprintf( '%s: Saving ASCII %s\n', mfname, whlfname )
    aformat         = sprintf( '%%d%s\\n', repmat( '\t%d', [ 1 size( whl, 2 ) - 1 ] ) );            % save as an ASCII file for backwards compatibility (more disk space though):
    fid             = fopen( whlfname, 'w' );
    fprintf( fid, aformat, whl' );
    fclose( fid );

    % save the mov
    fprintf( '%s: Saving MATLAB %s\n', mfname, movfname )
    save( movfname, 'mov' )
    
    % save the figure
    if graphics
        fprintf( '%s: Saving %s\n', mfname, figfname )
        print( figs( 1 ), '-dpng', [ figfname '.png' ] )
    end
    
end

return

% EOF

% 27-feb-18
filebase                    = datenum2filebase( { 'm637', 94 } );
mov                         = eeg2whl( filebase, 'Overwrite', 1 );
[ mat0 mat1 mov theta ]     = whl2states( filebase );

% 31-aug-18
filebase                    = datenum2filebase( { 'mC41', -24 } );
mov                         = eeg2whl( filebase, 'Overwrite', -2, 'graphics', 1 );


