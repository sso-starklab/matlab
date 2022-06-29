% plot_s2s              plot a single CCH with bars and bands
%
% call                  ah = plot_s2s( filebase, n12 )
%
% gets                  filebase        can be a full-path filebase or an s2s structure
%                       n12             cluster numbers; either a 2-element vector or a two-row matrix,
%                                           [ shank1 clu1; shank2 clu2 ]
%
% optional arguments (given as name/value pairs):
% 
%                       plotmode         0  raw CCHs,            w/ patches, w/ lines
%                                        1  subtracts predictor, w/ patches, w/ lines
%                                      {-1} raw CCHs,            no patches, w/ lines
%                                       -2  raw CCHs,            no patches, no lines
%                                       -3  subtracts predictor, no patches, no lines
%                                       -4  subtracts predictor, no patches, w/ lines
%
%                       convType    { 'gauss' }
%                       suffix      { 's2s' }
% 
% returns               ah              handled to axes
%
% does                  plots in the current axes
% 
% calls                 ParseArgPairs, replacetok               (general)
%
% see also spikes2spikes

% 15-jan-12 ES

% revisions
% 24-jan-12 argument checking for cluster existence
% 30-jan-12 plotmode -1 for no patches
% 02-mar-12 FileBase can also be a 2-element cell array (datestr, fnum)
% 12-jan-13 slash_ -> replacetok
% 15-jan-13 in case plotmode is NaN, just return the [ bins cch pred ], no plot 
% 16-may-13 suffix added
% 23-feb-19 (1) bar modified from implicit 0.8 to 1
%           (2) added plotmode -4
% 16-sep-19 cleaned up; argument and file handling improved

function ah = plot_s2s( filebase, n12, varargin )

% constants
LS = '--';

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( n12 )
    n12                     = [ 1 2 ];
end
[ plotmode, convType, suffix ...
    ]                       = ParseArgPairs(...
    { 'plotmode', 'convType', 'suffix' ...
    }...
    , { -1, 'gauss', 's2s' ...
    }...
    , varargin{ : } );

% files
if isa( filebase, 'char' )
    s2sfname = [ filebase '.' suffix ];
    if strcmp( convType, 'jitter' )
        s2sfname            = [ s2sfname '.jit' ];
    end
    if exist( s2sfname, 'file' )
        load( s2sfname, 's2s', '-mat' )
    end
elseif isa( filebase, 'struct' )
    s2s                     = filebase;
end
if ~exist( 's2s', 'var' )
    error( 'input type mismatch' )
end
[ ~, fname, fsuffix ]        = fileparts( s2s.filebase );


% determine pair
if numel( n12 ) == 2
    n1                      = n12( 1 );
    n2                      = n12( 2 );
elseif isequal( [ 2 2 ], size( n12 ) )
    n1                      = find( s2s.shankclu( :, 1 ) == n12( 1, 1 ) & s2s.shankclu( :, 2 ) == n12( 1, 2 ) );
    n2                      = find( s2s.shankclu( :, 1 ) == n12( 2, 1 ) & s2s.shankclu( :, 2 ) == n12( 2, 2 ) );
else
    error( 'input size mismatch' )
end
if isempty( n1 ) || isempty( n2 ) || n1 > size( s2s.ccg, 2 ) || n2 > size( s2s.ccg, 3 )
    fprintf( 1, 'missing clusters\n' )
    ah                      = NaN;
    return
end

%---------------------------------------------------------------%
% precompute bands
cchbins                     = s2s.t;
cch                         = s2s.ccg( :, n1, n2 );
pred                        = s2s.pred( :, n1, n2 );
alpha                       = s2s.alpha;
nBonf                       = sum( s2s.t_ROI );
upperpt                     = poissinv( ( 1 - alpha ) * ones( size( pred ) ), pred );           % no correction
lowerpt                     = poissinv(   alpha * ones( size( pred ) ), pred );
uppergb                     = poissinv( ( 1 - alpha / nBonf ) * ones( size( pred ) ), pred );
lowergb                     = poissinv(  ( alpha / nBonf ) * ones( size( pred ) ), pred );      % bonferroni corrected for the number of bins of interest
gbUpper                     = poissinv( 1 - alpha / nBonf, max( pred(s2s.t_ROI), [], 1 ) );
gbLower                     = poissinv( alpha / nBonf, min( pred(s2s.t_ROI), [], 1 ) );
if plotmode == 1 || plotmode == -3 || plotmode == -4
    pred0                   = pred;                                                             % subtract the predictor
else
    pred0                   = zeros( size( pred ) );
end

if isnan( plotmode )
    ah                      = [ cchbins cch pred ];
    return
end
    
%---------------------------------------------------------------%
% plot
newplot
hold on
ah                          = gca;

% plot the CCH
bar( cchbins, cch - pred0, 1, 'facecolor', 'k', 'edgecolor', 'k' )
set( gca, 'XLim', [ min( cchbins ), max( cchbins ) ] )
set( gca, 'box', 'off', 'tickdir', 'out' )
xlabel( 'Time [ms]' )
ylabel( 'Counts' )
if plotmode == -2
    title( sprintf( '%d.%d x %d.%d'...
        , s2s.shankclu( n1, 1 ), s2s.shankclu( n1, 2 ), s2s.shankclu( n2, 1 ), s2s.shankclu( n2, 2 ) ) );
else
    title( sprintf( '%s%s: %d.%d x %d.%d', replacetok( fname, '\_', '_' ), fsuffix...
        , s2s.shankclu( n1, 1 ), s2s.shankclu( n1, 2 ), s2s.shankclu( n2, 1 ), s2s.shankclu( n2, 2 ) ) );
end
line( [ 0 0 ], ylim, 'color', [ 0 0 0], 'linestyle', '--' );

% add lines
if plotmode > -2 || plotmode == -4
    line( cchbins( s2s.t_ROI ), gbUpper * ones( sum( s2s.t_ROI ), 1 ) - pred0( s2s.t_ROI ), 'color', 'b' )
    line( cchbins( s2s.t_ROI ), gbLower * ones( sum( s2s.t_ROI ), 1 ) - pred0( s2s.t_ROI ), 'color', 'b' )
    line( cchbins, pred - pred0,    'linestyle', LS, 'color', 'b' )
    line( cchbins, upperpt - pred0, 'linestyle', LS, 'color', 'r' )
    line( cchbins, lowerpt - pred0, 'linestyle', LS, 'color', 'r' )
    line( cchbins, uppergb - pred0, 'linestyle', LS, 'color', 'm' )
    line( cchbins, lowergb - pred0, 'linestyle', LS, 'color', 'm' )
end
% highlight bars
if plotmode >= 0
    bh1 = bar( cchbins, s2s.hiBins( :, n1, n2 ) * max( ylim ), 1, 'edgecolor', 'b', 'facecolor', 'b' );
    bh2 = bar( cchbins, s2s.loBins( :, n1, n2 ) * max( ylim ), 1, 'edgecolor', 'r', 'facecolor', 'r' );
    set( get( bh1, 'Children' ), 'faceAlpha', 0.3, 'edgeAlpha', 0 )
    set( get( bh2, 'Children' ), 'faceAlpha', 0.3, 'edgeAlpha', 0 )
end
if plotmode == -1
    bh1 = bar( cchbins, s2s.hiBins( :, n1, n2 ) * min( ylim ), 1, 'edgecolor', 'b', 'facecolor', 'b' );
    bh2 = bar( cchbins, s2s.loBins( :, n1, n2 ) * min( ylim ), 1, 'edgecolor', 'r', 'facecolor', 'r' );
    set( get( bh1, 'Children' ), 'faceAlpha', 0.3, 'edgeAlpha', 0 )
    set( get( bh2, 'Children' ), 'faceAlpha', 0.3, 'edgeAlpha', 0 )
end
% add calibration
if plotmode == -2
    lh = line( cchbins( 1 ) / 2 + [ 10 10 20 ], 1.5 * mean( ylim ) * [ 1 1 1 ] + [ 20 0 0 ] ); % 20 counts, 10 ms calibration
    set( lh, 'color', [ 0 0.7 0 ], 'linewidth', 2 )
end

hold off

return

% EOF

% example:
filebase = '/Volumes/Data/phaser3/mouse428/26nov11/dat/es26nov11_1/es26nov11_1';
s2s = spikes2spikes( filebase );
figure, plot_s2s( s2s, [ 12 6 ], 'plotmode', -1 );
figure, plot_s2s( filebase, [ 1 13; 1 7 ], 'plotmode', 1 );

