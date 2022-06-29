% FIG_OUT       enlarge, save, or print a figure.
%
% CALL          fig_out( f, mode, fname, ftype, MAXFIG, verbose, MOD )
%
% GETS          f               figure handle (gcf)
%               mode            {0} view
%                                1  save
%                                2  print
%               fname           {'temp'}; full path and file name (if to be saved)
%               ftype           {'jpg'}; supported: 'jpg', 'jpeg', 'fig', 'eps', 'bmp', 'png', 'tif', 'pdf', 'print'
%               MAXFIG          { [ 1 54 1024 647 ] }; [pixels], resize to MAXFIG 
%               verbose         {0} 
%               MOD             {''}; argument to print.m
%
% DOES          resizes figure and optionally prints (to file to printer)

% 15-dec-02 ES

% revisions
% 29-dec-12 standalone
% 28-jun-18 modified f -> gcf; added MOD
% 30-aug-18 removed MOD to allow compatibility with R2018a
% 18-aug-19 cleaned up; added MOD as seventh argument
% 19-oct-19 modified default behavior of pdf

function fig_out( f, mode, fname, ftype, MAXFIG, verbose, MOD )

% constants
PAPTYPE                 = 'A4';

% arguments
nargs = nargin;
if nargs < 1 || isempty(f)
    f                   = gcf; 
end
if nargs < 2 || isempty(mode)
    mode                = 0; 
end
if nargs < 3 && mode==1
    disp( 'figure saved in pwd as temp.jpg' )
    fname               = 'temp'; 
end
if nargs < 4 || isempty(ftype)
    ftype               = 'jpg'; 
end
if nargs < 5 || isempty(MAXFIG)
    MAXFIG              = [ 1 54 1024 647 ]; % pixels
end
if nargs < 6 || isempty( verbose )
    verbose             = 0;
end
if nargs < 7 || isempty( MOD )
    MOD                 = '-bestfit';
end
if ~ismember( MOD, { '', '-bestfit' } )
    MOD                 = '';
end

[do_FIG, do_JPG, do_EPS, do_BMP, do_PRINT, do_PNG, do_TIF, do_PDF, do_NOTHING ] = local_parse_input( ftype );

f0                      = gcf;
figure( f );

switch mode
    case 0                                                              % enlarge
        set(f,'Position',MAXFIG);
    case 1                                                              % save
        if do_NOTHING
            return
        end
        if do_FIG                                                       % all data
            saveas(fname,'fig'),
            if verbose
                disp(['                 saved figure as ' fname '.fig'])
            end
        end
        if do_PDF                                                       % vector data
            set( gcf, 'Position', MAXFIG );
            print( '-dpdf', fname, MOD )
            if verbose
                disp(['                 saved figure as ' fname '.pdf'])
            end
        end
        if do_JPG                                                       % highest compression ratio
            set( gcf, 'Position', MAXFIG );
            print( '-djpeg', fname )
            if verbose
                disp(['                 saved figure as ' fname '.jpg'])
            end
        end
        if do_EPS                                                       % high resolution
            set( gcf, 'Position', MAXFIG );
            print( '-depsc', fname, MOD )
            if verbose
                disp(['                 saved figure as ' fname '.eps'])
            end
        end
        if do_BMP                                                       % enlarge first: ~650KB ( normal size - ~250KB )
            set( gcf, 'Position', MAXFIG );
            pause( 0.1 )
            print( '-dbitmap', fname )
            if verbose
                disp(['                 saved figure as ' fname '.bmp'])
            end
        end
        if do_PNG                                                       % enlarge first: ~650KB ( normal size - ~250KB )
            set( gcf, 'Position', MAXFIG );
            print( '-dpng', fname )
            if verbose
                disp(['                 saved figure as ' fname '.png'])
            end
        end
        if do_TIF                                                       % enlarge first: ~650KB ( normal size - ~250KB )
            set( gcf, 'Position', MAXFIG, MOD );
            print( '-dtiff', fname )
            if verbose
                disp(['                 saved figure as ' fname '.tif'])
            end
        end
        if do_PRINT
            set( 'PaperType', PAPTYPE, MOD );
            print -dwinc
            if verbose
                disp( '                 printed figure.' )
            end
        end
    case 2                                                              % print on A4
        set( gcf, 'PaperType', PAPTYPE );
        print -dwinc
end

figure( f0 );


return

%----------------------------------------------------------------%
% local_parse_input
function [ FIG, JPG, EPS, BMP, PRINT, PNG, TIF, PDF, NONE ] = local_parse_input( mode )

FIG                     = 0; 
JPG                     = 0; 
EPS                     = 0; 
BMP                     = 0; 
PRINT                   = 0; 
PNG                     = 0; 
TIF                     = 0; 
PDF                     = 0; 
NONE                    = 0;

parsed = local_parse_str( lower( mode ) );

for i                   = 1 : length( parsed )
    switch lower( parsed{ i } )
        case { 'jpg' 'jpeg' }
            JPG         = 1;
        case 'fig'
            FIG         = 1;
        case 'eps'
            EPS         = 1;
        case 'bmp'
            BMP         = 1;
        case 'png'
            PNG         = 1;
        case 'tif'
            TIF         = 1;
        case 'pdf'
            PDF         = 1;
        case 'print'
            PRINT       = 1;
        case 'none'
            NONE        = 1;
    end
end

return

%----------------------------------------------------------------%
% PARSE_STR     parses a string by underscores.
%
%               PSTR = PARES_STR(STR)
%
%               STR     a string
%               PSTR    cell array

function pstr = local_parse_str( str )

delims                  = find( str == '_' );
nstr                    = length( delims ) + 1;
j                       = 0;
pstr                    = cell( 1, nstr );
for i                   = 1 : nstr
    j                   = j + 1;
    if i == 1
        first           = 1; 
    else
        first           = delims( i - 1 ) + 1;
    end
    if i == nstr
        last            = length( str ); 
    else
        last            = delims( i ) - 1;
    end
    pstr{ j }           = str( first : last );
    if isempty( pstr{ j } )
        pstr( j )       = []; 
        j               = j - 1; 
    end
end
pstr                    = pstr( 1 : j );

return

% EOF
