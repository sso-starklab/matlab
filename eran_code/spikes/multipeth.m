% multipeth             compute PETH for multiple units and triggers
% 
% [ peth, bins, ntrig ] = multipeth( clu, res, trig, tim )
%
% clu, res      must be vectors of the same length
% trig, tim     must be vectors of the same length
% clu, trig     [labels], may overlap
% res, tim      [samples], same Fs
% 
% additional arguments, specified as parameter name/value pairs:
% binsize       {20}; [samples]
% halfwin       {50}; [bins]
% scale         {'hz'} (spike/sec/trigger) or 'count' (number of spikes)
% sdGauss       {0}; [samples] of Gaussian fir used for smoothing PETH
% graphics      {0}; 1 generates rudimentary graphical display
% clucat        {[]}; [ cluster category ]
% Fs            {20000}; Hz
%
% this routine simply computes all possible PETHs and organizes them in a
% 3D array: nbins x nclu x ntrig. the array is thus m x n x p, where
%   m = 2 * ceil( halfwin / binsize )
%   n = length( unique( clu ) )
%   p = length( unique( trig ) )
%
% calls         ParseArgPairs, uhist                        (general)
%               inranges                                    (sets)
%               CCG                                         (blab)
%               makegaussfir, firfilt                       (ssp)
%               imagescbar, myjet                           (graph)
%
% See also      get_spikes, get_triggers, get_rasters

% 14-feb-13 ES

% revisions
% 18-feb-13 modified graphics slightly to include sorting by cell-type +
%               shank number
% 04-mar-13 (1) implemented period selection properly
%           (2) proper argument parsing w/ ParseArgPairs
% 18-aug-19 cleaned up

function [ peth, bins, ntrig, hh, ah ] = multipeth( clu, res, trig, tim, varargin )

% constants
peth                        = [];
bins                        = [];
ntrig                       = [];
hh                          = [];
ah                          = [];

% arguments
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 4 || isempty( clu ) || isempty( res ) || isempty( tim )
    return
end
if isempty( trig )
    trig                    = ones( size( tim ) );
end
if length( clu ) ~= length( res ) || length( trig ) ~= length( tim )
    fprintf( '%s: input size mismatch\n', mfname );
    return
end
[ binsize, halfwin, scale, sdGauss, graphics, clucat, Fs ] = ParseArgPairs(...
    { 'binsize', 'halfwin', 'scale', 'sdGauss', 'graphics', 'clucat', 'Fs' }...
    , { 20, 50, 'count', 0, 0, [], 20000 }...
    , varargin{ : } );
if isempty( clucat ) || size( clucat, 2 ) ~= 2
    clunums                 = unique( clu );
    clucat                  = [ clunums NaN * ones( size( clunums ) ) ];
end
if nargout == 0
    graphics                = 1;
    sdGauss                 = 1;
end
if length( halfwin ) > 1
    halfwin0                = halfwin;
    halfwin                 = max( abs( halfwin0 ) );
else
    halfwin0                = [];
end

% offset the trigger indices to prevent overlaps
fac                         = ceil( max( clu ) / 1e3 ) * 1e3;
trig                        = trig + fac;
[ ntrig, utrig ]            = uhist( trig );
m                           = 2 * halfwin + 1;
p                           = length( utrig );

% keep only the relevant data to speed up computations
clucat0                     = clucat;
periods                     = tim * [ 1 1 ] + halfwin * binsize * ( ones( length( tim ), 1 ) * [ -1 1 ] );
kidx                        = inranges( res, periods );
res                         = res( kidx );
clu                         = clu( kidx );
kidx                        = ismember( clu, clucat( :, 1 ) );
res                         = res( kidx );
clu                         = clu( kidx );
uclu                        = unique( clu );
n                           = length( uclu );
clucat( ~ismember( clucat( :, 1 ), uclu ), : ) = [];

% merge spikes and triggers
Res                         = [ res; tim ];
Clu                         = [ clu; trig ];
uClu                        = unique( Clu );
uidx                        = ismember( uClu, uclu ); 
tidx                        = ismember( uClu, utrig );

% compute
fprintf( 1, '%s: called w/ %d/%d spikes/units, %d/%d triggers/caterories, %0.3g ms bins and %0.3g sec total window\n'...
    , mfname, length( clu ), n, length( trig ), p...
    ,  binsize / Fs * 1000, ( 2 * halfwin + 1 ) * binsize / Fs );
[ ccg, bins ]               = CCG( Res, Clu, binsize, halfwin, Fs, uClu, 'count' );
bins                        = bins / 1000; % ms -> s

% condition (extract, scale, smooth, expand)
peth                        = ccg( :, uidx, : );               % bins x clu x trig, but reversed
peth                        = peth( :, :, tidx );
peth                        = peth( m : -1 : 1, :, : );        % trigger on trig
switch lower( scale )
    case 'count'
    case { 'hz', 'spikes/sec', 'spk/s' }
        for i               = 1 : p
            peth( :, :, i ) = peth( :, :, i ) / ntrig( i ) / binsize * Fs;
        end
        scale               = 'Spikes/s';
end
if sdGauss 
    g                       = makegaussfir( sdGauss, 1 ); 
    for i                   = 1 : p
        peth( :, :, i )     = firfilt( peth( :, :, i ), g );
    end
end
cidx                        = ismember( clucat0( :, 1 ), clucat( :, 1 ) );
peth0                       = peth;
peth                        = zeros( m, size( clucat0, 1 ), p );
peth( :, cidx, : )          = peth0;
if ~isempty( halfwin0 )
    kidx                    = ismember( -halfwin : halfwin, halfwin0( 1 ) : halfwin0( 2 ) );
    peth( ~kidx, : )        = [];
    bins( ~kidx )           = [];
end

% graphics
if graphics
    clucat                  = clucat0;
    if all( isnan( clucat( :, 2 ) ) )
        clucat( :, 2 )      = 1;
    end
    [ ~, idx ]              = sort( clucat( :, 2 ) );
    if isempty( setdiff( unique( clucat( :, 2 ) ), [ 0 1 ] ) )
        b1                  = sum( clucat( :, 2 ) == 0 ) + 0.5;
    end
    for i                   = 1 : p
        if p == 1
            newplot
        else
            subplot( ceil( sqrt( p ) ), ceil( p / ceil( sqrt( p ) ) ), i )
        end
        mat                 = peth( :, idx, i );
        [ hh, ah ]          = imagescbar( bins, 1 : size( mat, 2 ), mat );
        subplot( ah( 1 ) )
        alines( 0, 'x', 'color', [ 1 0 0 ], 'linestyle', '--' );
        for j               = 1 : length( ah )
            subplot( ah( j ) )
            if exist( 'b1', 'var' )
                alines( b1, 'y', 'color', [ 1 0 0 ] );
            end
            set( gca, 'tickdir', 'out', 'box', 'off' )
        end
        subplot( ah( 1 ) )
        title( sprintf( 'T%s; n%d; m%0.3g', num2str( utrig( i ) ), ntrig( i ), max( mat( : ) ) ) )      
        if i == 1 
            ylabel( 'Cell #' )
            xlabel( sprintf( 'Time [s]' ) )
        end
    end
    colormap( myjet )
end

return

% EOF
