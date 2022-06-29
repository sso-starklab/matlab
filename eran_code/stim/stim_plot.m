% stim_plot         summarize graphically
%
% CALL              fig = stim_plot( stims, figname, savetype, Overwrite )
% 
% GETS              stims           structure
%                   figname         name (without full path; will be taken from stims( 1 ).filebase)
%                   savetype        {'png'}
%                   Overwrite       {0}
%
% CALLS             tilefig, PlotStims, replacetok, uhist, minmax, fig_out

% 25-jan-13 ES

% revisions
% 10-apr-13 modified to have only one legend on the last subplot
% 06-sep-13 1. excluded some outliers on x-axis
%           2. support of trigger
% 30-aug-18 changed call to legend to allow compatibility with R2018a
% 18-aug-19 cleaned up
% 07-sep-19 default figname created
% 10-sep-19 figdir directory issues handled; savetype not forced non-empty
% 12-sep-19 figfname modified to support saving in PC

function fig = stim_plot( stims, figname, savetype, Overwrite )

fig                         = [];

nargs                       = nargin;
if nargs < 1 || isempty( stims )
    return
end
if nargs < 2 || isempty( figname )
    figname                 = '';
end
if nargs < 3
    savetype                = 'png';
end
if nargs < 4 || isempty( Overwrite )
    Overwrite               = 0;
end
  
% determine filename
mfname                              = upper( mfilename );
[ pathname, filename, extname ]     = fileparts( stims( 1 ).filebase );
filename                            = [ filename extname ];

% determine where to save
if isempty( figname ) || ~exist( fileparts( figname ), 'dir' )
    homedir                 = fileparts( fileparts( pathname ) );
    figdir                  = sprintf( '%s/figs', homedir );
    if exist( homedir, 'dir' ) && ~exist( figdir, 'dir' )
        try
            mkdir( homedir, 'figs' )
        catch
            savetype        = '';
        end
    elseif ~exist( figdir, 'dir' )
        savetype            = '';
    end
    if ~isempty( savetype )
        figname             = sprintf( '%s/%s.stm.sim', figdir, [ filename extname ] );
    end
end

% accumulate
nstims                      = length( stims );
n                           = zeros( 1, nstims );
nc                          = zeros( 1, nstims );
types                       = [];
for i                       = 1 : nstims
    n( i )                  = length( stims( i ).types ); 
    nc( i )                 = length( stims( i ).chan ); 
    types                   = [ types; stims( i ).types ];
    if strcmp( stims( i ).source, 'trigger' )
        nc( i ) = -1;
    end
end
rmv                         = n == 0 | nc == 0;
stims( rmv )                = [];
nc( rmv )                   = [];
n( rmv )                    = [];
nstims                      = sum( ~rmv );
nchans                      = sum( nc == 1 );
utypes                      = unique( types );
colors                      = lines( length( utypes ) );
nsp                         = nstims + 1;

% plot
if nsp < 3
    alignment               = 'center';
else
    alignment               = 'topright';
end
[ ah, fig ]                 = tilefig( ceil( sqrt( nsp ) ), ceil( ( nsp ) / ceil( sqrt( nsp ) ) ), 1, 0.8, alignment );
xlims                       = zeros( nstims, 2 );
ylims                       = xlims;
vals                        = [];
durs                        = [];
adurs                       = [];
for i                       = 1 : nstims
    if isempty( stims( i ).durs )
        continue
    end
    vals                    = [ vals; max( stims( i ).vals ) ];
    durs                    = [ durs; max( stims( i ).durs ) ];
    adurs                   = [ adurs; stims( i ).durs ];
    if ischar( stims( i ).source )
        sourcetype{ i }     = stims( i ).source;
    else
        sourcetype{ i }     = 'sim';
    end
    subplot( ah( i ) );
    ph                      = PlotStims( stims( i ).durs, stims( i ).vals, stims( i ).types );
    [ ~,cidx ]              = ismember( unique( stims( i ).types ), utypes );
    for k                   = 1 : length( ph )
        set( ph( k ), 'color', colors( cidx( k ), : ) );
    end
    xlims( i, : )           = xlim;
    ylims( i, : )           = ylim;
end
if isempty( ylims )
    textf( 0.5, 0.975, replacetok( filename, '\_', '_' ) );
    return
end

% take care of axes - same for all subplots
[ nsources, usources ]      = uhist( sourcetype );
if sum( nsources > 1 ) == 1
    sourcetype( ismember( sourcetype, 'sim' ) ) = usources( find( nsources == max( nsources ), 1 ) );
end
usources                    = unique( sourcetype );
nusources                   = length( usources );
xlimsNew                    = zeros( nusources, 2 );
ylimsNew                    = zeros( nusources, 2 );
for k                       = 1 : nusources
    idx                     = ismember( sourcetype, usources{ k } );
    xlimsNew( k, : )        = minmax( xlims( idx, : ) );
    ylimsNew( k, : )        = minmax( ylims( idx, : ) );
    xlimsNew( k, 1 )        = max( 0, xlimsNew( k, 1 ) );
    xlimsNew( k, 2 )        = min( max( durs ) * 1.2, xlimsNew( k, 2 ) );
    ylimsNew( k, 1 )        = max( 0, ylimsNew( k, 1 ) );
    ylimsNew( k, 2 )        = min( max( vals ) * 1.2, ylimsNew( k, 2 ) );   
end
for i = 1 : nstims
    if isempty( stims( i ).durs )
        continue
    end
    subplot( ah( i ) );
    xlims                   = sort( xlimsNew( ismember( usources, sourcetype{ i } ), : ) );
    ylims                   = sort( ylimsNew( ismember( usources, sourcetype{ i } ), : ) );
    sdurs                   = sort( adurs ); 
    xlims( 2 )              = min( xlims( 2 ), max( floor( length( sdurs ) * 0.999 ), 1 ) );
    if xlims( 2 ) < xlims( 1 )
        continue
    end
    set( gca, 'xlim', xlims, 'ylim', ylims );   
    if i ~= nstims
        xlabel( '' )
        ylabel( '' )
    end
    nevents                 = length( stims( i ).types );
    channel                 = stims( i ).chan;
    tx                      = min( xlims ) + diff( xlim ) * 0.5; 
    ty                      = min( ylims ) + diff( ylim ) * 0.9; 
    if nc( i ) == 1
        th                  = text( tx, ty, sprintf( '%s; %d', num3str( channel ), nevents ) );
    elseif nc( i ) == -1
        th                  = text( tx, ty, sprintf( 'T%s; %d', num3str( channel ), nevents ) );
    else
        th                  = text( tx, ty, sprintf( 'sim; %d', nevents ) );
    end
    set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top' )
    legend( 'off' )
end

% keep only one legend on the last subplot
subplot( ah( nstims + 1 ) )
hold on
for k                   = 1 : length( utypes )
    ph( k )             = plot( 2 * xlims(2 ), 2 * ylims( 2 ), '.' );
    set( ph( k ), 'color', colors( k, : ) );
end
th                      = title( replacetok( filename, '\_', '_' ) );
set( th,'Fontsize', 12 )
lh                      = legend( ph, utypes );
set( lh, 'Box', 'off' )
if nstims > 1
    if xlims( 2 ) > xlims( 1 )
        set( ah( nstims + 1 ), 'xlim', xlims, 'ylim', ylims );
    end
    for i               = nstims + 1 : size( ah, 1 )
        axis( ah( i ), 'off' )
    end
end

% save
if exist( 'savetype', 'var' ) && ~isempty( savetype )
    if Overwrite  >= 1 || ( Overwrite ~= -1 && ~exist( [ figname '.' savetype ], 'file' ) )
        fprintf( 1, '%s: saving %s.%s\n', mfname, figname, savetype )
        figfname = [ figname '.' savetype ];
        if ispc
            figfname    = replacetok( figfname, '\', '/' );
        end
        print( fig, [ figfname '.png' ], [ '-d' savetype ] );
    end
end

return

% EOF
