% multipeth_plot            plot peths for multiple units and trigger events
%
% CALL          [ fig, ah ] = multipeth_plot( peth, bins, shankclu )
%
% GETS          peth          m x n x p array: bins x units x stimuli, [spikes/sec]
%               bins          bins: m x 1 array [sec]
%               shankclu      n x 3 array: [ shank clu cell_type ]; cell_type: 0-INT; 1-PYR
%
% OPTIONAL ARGUMENTS (name/value pairs):
% 
%               trigs         optional q x 1 array of labels
%               durs          optional q x 1 array; [sec]
%               vals          optional q x 1 array; [V]
%               str           optional string to be added to the figure
%               figname       optional full path and name to be saved as
%               savetype      {'png'}
%               ahandle       optional axis handle to be partitioned into 4 (single trig only)
%
% DOES          for a simple set of peth/bins/shankclu, plots all peths corresponding to the same
%               stimuli as one set of 4 subplots: PYR mean, PYR SUs (scaled 0-1), and the
%               same for INT.
%               if labels/durations/values/st are given, these are used to add details
%               if figname/savetype are given, the output is also saved
%
% CALLS         ParseArgPairs                   (general)
%               calc_sem, moutliers             (stats)
%               addpatch, alines, calibration, fig_out, imagescbar, patch_band, textf (graph)
%
% See also      get_spikes, get_triggers, multipeth, multipeth_make

% 26-feb-13 ES

% revisions
% 16-sep-13 do not plot unclassified (NaN) units
% 03-nov-13 (1) argument pairs supported
%           (2) scaled images/bands supported
%           (3) calibration proper
%           (4) enable plot single combination on externally-provided space
% 19-sep-14 minor modification in logic during saving
% 18-aug-19 cleaned up

% NOTE
% this presents the peth in both 0-max of each unit and the average
% firing rate of all units. ideally each unit should also have a bar at the
% side to indicate the max rate (by which it is scaled), like I did with
% the coherence data

function [ fig, ah ] = multipeth_plot( peth, bins, shankclu, varargin )

% constants
f                           = 0.85;                                     % determines the covered fig portion
COLORS                      = [ 0 0 0.7; 1 0 0 ];                       % INT/PYR coloring
TRIGCOLOR                   = [ 0.7 0.7 1 ];                            % patch
blackColor                  = [ 0 0 0 ];
sepStyle                    = '--';
CALIBCOLOR                  = [ 0 0 0 ];                                % calibration
YCALIB                      = [ 10 1 ];                                 % spks/s: 1 for pyr, 10 for IN
cLW                         = 2;                                        % calibration line width
LW                          = 0.5;                                      % separators between shanks
bColor                      = [ 1 0 0 ];                                % line color
FS                          = 10;
bLW                         = 2;                                        % line width for bar
tColor                      = [ 1 0 1 ];                                % text color
tFW                         = 'normal';
tFS                         = 12;
tY                          = 0.985;                                    % title

%----------------------------------------------------------------------%
% check the input and extract some parameters
nargs                       = nargin;
if nargs < 3
    error( 'missing arguments' )
end
if isempty( peth ) || isempty( bins ) || isempty( shankclu )
    return
end
[ trigs, durs, vals, str, figname, savetype...
    , ahandle, imageType, barType, avMode ] = ParseArgPairs(...
    { 'trigs', 'durs', 'vals', 'str', 'figname', 'savetype'...
    , 'ahandle', 'imageType', 'barType', 'avMode' }...
    , { [], [], [], '', '', 'png'...
    , [], 'imagesc', 'bar', 'mean' } ...
    , varargin{ : } );

% peth details
[ m, n, p ]                 = size( peth );
if n ~= size( shankclu, 1 )
    error( 'fix multipeth indexing...' )
end
if size( shankclu, 2 ) < 3
    error( 'fix shankclu format...' )
end
if m ~= length( bins )
    error( 'bins/peth mismatch..' )
end
bins                        = bins( : );
rmv                         = isnan( shankclu( :, 3 ) );
if sum( rmv )
    n                       = sum( ~rmv );
    shankclu( rmv, : )      = [];
    peth( :, rmv, : )       = [];
    
end
ctypes                      = unique( shankclu( :, 3 )' );

% stimulus details
if isempty( trigs ) || isempty( durs ) || isempty( vals )
    dflag                   = 0;
else
    dflag                   = 1;
    [ utrigs, ~, itrigs ]   = unique( trigs );
    ntrigs                  = length( utrigs );
    if ntrigs == 1
        stimreps            = length( trigs );
    else
        stimreps            = hist( trigs, utrigs );
    end
    chans                   = round( trigs / 1e9 );
    [ ~, ~, cats ]          = unique( rem( trigs, 1e9 ) );
    if isequaln( cats, 1 )
        cats                = NaN;
    end
    combs                   = [ chans cats ];
    ucombs                  = unique( combs, 'rows' );
    uchans                  = unique( chans );
    if 1 %collapse % temporary fix..: give all sim. stimuli the same channel allocation
        %        fprintf( '%s: This is a temporary fix!!\n', upper( mfilename ) )
        cidx                = uchans < 32;
        chans( ismember( chans, uchans( cidx ) ) ) = 1;
        ridx                = ismember( ucombs( :, 1 ), uchans( cidx ) );
        ucombs( ridx, 1 )   = 1;
        ucombs( ridx, 2 )   = 1 : sum( ridx );
        uchans( cidx )      = 1;
    end
    ncombs                  = uhist( ucombs( :, 1 ) )';
    uchans                  = unique( chans );
    nchans                  = length( uchans );
    if p ~= ntrigs
        error( 'fix trigger indexing...' )
    end
    stimdurs                = zeros( ntrigs, 1 );
    stimvals                = zeros( ntrigs, 1 );
    for i                   = 1 : length( utrigs )
        idx = itrigs == i;
        stimdurs( i )       = median( durs( idx ) );
        stimvals( i )       = mean( vals( idx ) );
    end
end
try
    gpos                    = get( ahandle, 'position' );
catch
    gpos                    = [];
end
if isempty( gpos )
    newfig                  = 1;
else
    newfig                  = 0;
end

%----------------------------------------------------------------------%
% tile a figure with the proper number of axes sets
% each set will show the multi-unit PETH in four subplots -
% PYR mean histogram (spikes/sec), PYR individual histograms (non averaged)
% and the same for INT

if dflag && nchans > 1
    nrows                   = nchans; % separate row for each channel
    ncols                   = max( [ ncombs; 2 ] );
else
    pp                      = max( p, 1 );
    nrows                   = ceil( sqrt( pp ) );
    ncols                   = ceil( pp / nrows );
end
if p <= 2
    alignmode               = 'center';
else
    alignmode               = 'edge';
end
pp                          = nrows * ncols;
wx                          = 1 / ncols * f;
wy                          = 1 / nrows * f;
switch alignmode
    case 'center'
        x0                  = ( 1 - f ) / 2 / ncols : 1 / ncols : ( 1 - ( 1 - f ) / 2 / ncols ); % align center
        y0                  = ( 1 - f ) / 2 / nrows : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align center
    case 'edge'
        y0                  = 0 : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align bottom
        x0                  = fliplr( 1 : -1 / ncols : ( 1 - f ) / 2 / ncols ) - wx; % align right
end
x0                          = ones( nrows, 1 ) * x0( 1 : ncols );
y0                          = flipud( y0( 1 : nrows )' * ones( 1, ncols ) );
wx                          = wx * ones( nrows, ncols );
wy                          = wy * ones( nrows, ncols );
x0                          = x0';
y0                          = y0';
wx                          = wx';
wy                          = wy';
ah                          = zeros( pp, 4 );
if newfig
    fig                     = figure;
end
for i                       = 1 : pp
    pos                     = [ x0( i ) y0( i ) wx( i ) wy( i ) ]; % split this into four: [ 1 3; 2 4 ]
    spos( 1, : )            = [ pos( 1 )                   pos( 2 ) + pos( 4 ) / 2     pos( 3 ) / 2    pos( 4 ) / 2 ];
    spos( 2, : )            = [ pos( 1 )                   pos( 2 )                    pos( 3 ) / 2    pos( 4 ) / 2 ];
    spos( 3, : )            = [ pos( 1 ) + pos( 3 ) / 2    pos( 2 ) + pos( 4 ) / 2     pos( 3 ) / 2    pos( 4 ) / 2 ];
    spos( 4, : )            = [ pos( 1 ) + pos( 3 ) / 2    pos( 2 )                    pos( 3 ) / 2    pos( 4 ) / 2 ];
    if newfig
        sposhat             = spos;
    else
        sposhat             = [ gpos( 1 ) + gpos( 3 ) * spos( :, 1 ) ...
            gpos( 2 ) + gpos( 4 ) * spos( :, 2 ) ...
            gpos( 3 ) * spos( :, 3 ) ...
            gpos( 4 ) * spos( :, 4 ) ];
    end
    for k = 1 : 4
        ah( i, k )          = axes( 'position', sposhat( k, : ) );
    end
end

%----------------------------------------------------------------------%
% plot one stimulus set at a time
uspi                        = zeros( p, 1 );
ahused                      = false( size( ah ) );
for k                       = 1 : length( ctypes )
    
    ctype                   = ctypes( k );
    idx                     = shankclu( :, 3 ) == ctype;
    
    for i                   = 1 : p
        % scale each unit 0-1; compute mean/SEM firing rate
        if dflag
            stimtime        = [ 0 stimdurs( i ) ];
        end
        peth0               = peth( :, idx, i );
        speth               = bsxfun( @rdivide, peth0, max( peth0 ) );
        speth( isnan( speth ) ) = 0;
        switch avMode
            case 'mean'
                mpeth       = nanmean( peth0, 2 );
                sem         = calc_sem( peth0, 2 );
            case 'moutliers'
                nSD = 3;
                mpeth       = moutliers( peth0, @nanmean, nSD, 0, 2, 2 );
                sem         = moutliers( peth0, @calc_sem, nSD, 0, 2, 2 );
            case 'median'
                mpeth       = nanmedian( peth0, 2 );
                sem         = mad( peth0, 1, 2 );
        end
        
        % subplot index
        if dflag && nchans > 1
            xi              = find( uchans == ucombs( i, 1 ) );
            yi              = find( ucombs( ucombs( :, 1 ) == ucombs( i, 1 ), 2 ) == ucombs( i, 2 ) );
            spi             = ncols * ( xi - 1 ) + yi;
        else
            spi             = i;
        end
        uspi( i )           = spi;
        
        % plot the mean peth
        subplot( ah( spi, 1 + 2 * ( ctype )  ) )
        ahused( spi, 1 + 2 * ( ctype ) ) = 1;
        if ismember( barType, { 'patch', 'patch_band' } )
            ylims           = [ 0 max( mpeth + sem ) * 1.1 ];
        else
            ylims           = [ 0 max( mpeth ) * 1.1 ];
        end
        if ylims( 2 ) == 0 || sum( isnan( ylims ) ) > 0
            axis off
            subplot( ah( spi, 2 + 2 * ( ctype )  ) )
            axis off
            continue;
        end
        if dflag && ~isnan( stimtime( 2 ) )
            addpatch( stimtime, ylims, TRIGCOLOR );
        end
        hold on
        switch barType
            case 'bar'
                bar( bins, mpeth, 'facecolor', COLORS( ctype + 1, : ), 'edgecolor', COLORS( ctype + 1, : ) )
            case { 'patch', 'patch_band' }
                patch_band( bins, mpeth, sem, COLORS( ctype + 1, : ) );
            case 'line'
                line( bins, mpeth, 'color', COLORS( ctype + 1, : ), 'linewidth', bLW )
                
        end
        xlims               = bins( [ 1 end ] );
        xlim( xlims )
        ylim( ylims )
        xlims               = xlim;
        ycalib              = YCALIB( k ) / ( ceil( YCALIB( k ) / max( ylims ) ) );
        if dflag && ~isnan( stimtime( 2 ) )
            calibration( [ stimtime( 2 ), ycalib ], [], [], [], 'color', CALIBCOLOR, 'linewidth', cLW );
        else
            calibration( [ 0 ycalib ], [], [], [], 'color', CALIBCOLOR, 'linewidth', cLW );
        end
        alines( 0, 'x', 'color', blackColor, 'linestyle', sepStyle );
        axis off
        
        % text info
        if dflag
            if ( k == 1 || length( ctypes ) == 1 ) && ~isnan( stimvals( i ) ) % value
                text( min( xlim ) + diff( xlim ) * 0.01 ...
                    , max( ylim ) + diff( ylim ) * 0.01 ...
                    , sprintf( '%0.2gV', stimvals( i ) )...
                    , 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'color', tColor );
            end
            if ( k == 2 || length( ctypes ) == 1 ) && ~isnan( stimreps( i ) ) && stimreps( i ) > 1 % n trigs
                text( max( xlim ) - diff( xlim ) * 0.01 ...
                    , max( ylim ) + diff( ylim ) * 0.01 ...
                    , sprintf( '%d', stimreps( i ) )...
                    , 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'color', tColor );
            end
        end
        
        % the scaled SU peths
        subplot( ah( spi, 2 + 2 * ( ctype ) ) )
        ahused( spi, 2 + 2 * ( ctype ) ) = 1;
        sh                  = shankclu( shankclu( :, 3 ) == ctype, 1 );
        borders             = find( diff( sh ) > 0 ) + 0.5;
        shnums              = sh( [ 1; ceil( borders ) ] );
        switch imageType
            case 'imagesc'
                imagesc( bins, 1 : size( speth, 2 ), speth' )
                axis ij
            case 'imagescbar'
                [ hhI, ahI ] = imagescbar( bins, 1 : size( speth, 2 ), peth0, [], 0 );
                bordersInv  = sum( shankclu( :, 3 ) == ctype ) - borders + 1;
                subplot( ahI( 2 ) )
                axis off
                alines( bordersInv, 'y', 'color', bColor, 'linewidth', LW );
                set( hhI( 2 ), 'facecolor', COLORS( ctype + 1, : ), 'edgecolor', COLORS( ctype + 1, : ) );
                curpos      = get( ahI( 1 ), 'position' );
                oldpos      = get( ah( spi, 1 + 2 * ( ctype )  ), 'position' );
                oldpos( 3 ) = curpos( 3 );
                set( ah( spi, 1 + 2 * ( ctype )  ), 'position', oldpos );
                subplot( ahI( 1 ) )
            case 'contourf'
                contourf( bins, 1 : size( speth, 2 ), speth', 64 );
                shading flat
                axis ij
        end
        alines( borders, 'y', 'color', bColor, 'linewidth', LW );
        if dflag && ~isnan( stimtime( 2 ) )
            alines( stimtime, 'x', 'color', blackColor, 'linestyle', sepStyle );
        else
            alines( 0, 'x', 'color', blackColor, 'linestyle', sepStyle );
        end
        
        % text info
        if k == 1 && dflag && ( nchans > 1 && yi == 1 || nchans == 1 && i == 1 ) && ~isnan( ucombs( i, 1 ) )% trigger channel
            text( min( xlim ) + diff( xlim ) * 0.1 ...
                , min( ylim ) + diff( ylim ) * 0.8 ... % 0.8 for ij, 0.2 for xy
                , sprintf( 'T%d', ucombs( i, 1 ) ) ...
                , 'HorizontalAlignment', 'left', 'color', tColor );
        end
        if i == 1 % elec groups
            if k == 1
                mult        = -0.05;
            else
                mult        = 1.05;
            end
            if strcmp( imageType, 'imagescbar' ) && ctype == 1
                subplot( ahI( 2 ) )
                for j = 1 : length( shnums )
                    th = text( diff( xlim ) * mult + min( xlim ), sum( shankclu( :, 3 ) == ctype ) - mean( find( sh == shnums( j ) ) ) + 1, sprintf( 'S%d', shnums( j ) ) );
                    set( th, 'color', bColor, 'fontsize', FS, 'HorizontalAlignment', 'center', 'FontWeight', 'normal' );
                end
                subplot( ahI( 1 ) )
            else
                for j = 1 : length( shnums )
                    th = text( diff( xlim ) * mult + min( xlim ), mean( find( sh == shnums( j ) ) ), sprintf( 'S%d', shnums( j ) ) );
                    set( th, 'color', bColor, 'fontsize', FS, 'HorizontalAlignment', 'center', 'FontWeight', 'normal' );
                end
            end
        end
        axis off
    end
    
end

th = textf( 0.5, tY, replacetok( str, '\_', '_' ) );
set( th, 'color', tColor, 'fontsize', tFS, 'fontweight', tFW )
cmap                        = flipud( gray );
colormap( cmap )

% clean up
ah0                         = ah( setdiff( 1 : pp, uspi ), : );
for i                       = 1 : length( ah0( : ) )
    subplot( ah0( i ) )
    axis off
end
for i                       = find( ~ahused( : ) ).'
    subplot( ah( i ) )
    axis off
end
if ~newfig
    subplot( ahandle )
    axis off
end

% save
if ~isempty( figname ) && ~any( isnan( savetype ) )
    fprintf( 'Saving %s..\n', figname )
    fig_out( fig, 1, [ figname '.' savetype ], savetype );
end

return

% EOF
