% calc_cch_width    calculate width and lag of a waveform (CCH, EPSP, PSTH)
%
% call              [ wid, lag ] = calc_cch_width( cch, t, method, graphics )
% 
% gets              cch         conditional rate (gain) cch, [spk/s]
%                               will also work with count CCH
%                   t           vector of same length as cch
%                   methods     1   FWHM:   lag to peak, FWHM
%                               2   CDF:    lag to median, 95% of the CDF
%                               3   COM:    lag to COM, SD from COM
%                               {4} Span:   lag to peak, span >0.05 of peak
%                               all methods work on the scaled (0-1) waveform
%                               methods 1 and 4 remove mean of negative
%                               time lags (i.e., assume 'causality')
%                   graphics    {0} on the present plot
%                       
% calls             calc_com, scale             (general)
%                   alines                      (graph)
%                   calc_fwhh                   (ssp)

% 19-jan-22 ES

function [ wid, lag ] = calc_cch_width( cch, t, method, graphics )

% constants
f                               = 0.05;                                     % relevant for methods 2 and 4 only

% arguments
nargs                           = nargin;
if nargs < 2 || isempty( cch ) || isempty( t )
    error( 'missing arguments' )
end
if ~isvector( cch ) || ~isequal( size( cch ), size( t ) )
    error( 'input size mismatch - must be two equal-sized vectors' )
end
if nargs < 3 || isempty( method )
    method                      = 4;
end
if nargs < 4 || isempty( graphics )
    graphics                    = 0;
end

% initialize output
wid                             = NaN;
lag                             = NaN;

% preliminaries
cch                             = cch( : );
t                               = t( : );
dt                              = diff( t( 1 : 2 ) );
x                               = scale( cch );
xz                              = x - mean( x( t < 0 ) );
if length( unique( xz ) ) == 1
    method                      = 0;
end

% compute
switch method
    case 1
        % method #1: FWHM (lag to peak, FWHM)
        [ ww, ll ]              = calc_fwhh( xz );
        lag                     = ll * dt + t( 1 ) - dt;
        wid                     = ww * dt;
    case 2
        % method #2: cumulative distribution (lag to median, 95% of the probability)
        cx                      = cumsum( x ) / sum( x );
        ll                      = find( cx > 0.5, 1, 'first' );
        bR                      = find( cx > ( 1 - f / 2 ), 1, 'first' );
        bL                      = find( cx > ( f / 2 ), 1, 'first' );
        lag                     = ll * dt + t( 1 ) - dt;
        wid                     = ( bR - bL + 1 ) * dt;
    case 3
        % method #3: COM (lag to COM, SD relative to COM)
        w                       = x;
        w( w < 0.05 )           = NaN;                                      % modified weight vector
        [ cc, ss ]              = calc_com( ( 1 : length( x ) )', w );
        lag                     = cc * dt + t( 1 ) - dt;
        wid                     = ss * dt;
    case 4
        % method #4: span (lag to peak, span of values above 5% of the peak)
        [ maxval, midx ]        = max( xz( t > 0 ) );
        midx                    = midx + find( t > 0, 1, 'first' ) - 1;
        if ~isnan( maxval )
            bL                  = midx - find( xz( midx : -1 : 1 ) < maxval * f, 1, 'first' ) + 2;
            bR                  = midx + find( xz( midx : end ) < maxval * f, 1, 'first' ) - 2;
            if isempty(bR)
                wid = NaN;
            end
            lag                 = midx * dt + t( 1 ) - dt;
            wid                 = ( bR - bL + 1 ) * dt;
        end
end

% plot
if ~graphics
    return
end
newplot
if all( cch < 0 )
    plot( t, cch, 'k' );
else
    bar( t, cch, 1, 'k' );
end
hold on
if exist( 'bR', 'var' ) && exist( 'bL', 'var' )
    xx                          = t( [ bL bR ] ) + [ -1 1 ]' * dt / 2;
else
    xx                          = lag + wid / 2 * [ -1 1 ];
end
yy                              = ylim;
ph                              = patch( xx( [ 1 2 2 1 ] ), yy( [ 1 1 2 2 ] ), [ 1 0 0 ] );
set( ph, 'FaceAlpha', 0.25, 'FaceColor', [ 1 0 0 ], 'EdgeColor', [ 1 0 0 ], 'EdgeAlpha', 0.25 )
set( gca, 'tickdir', 'out', 'box', 'off' );
if all( cch < 0 )
    plot( t, cch, 'k' );
else
    bh                          = bar( t, cch, 1, 'k' );
	set( bh, 'FaceColor', [ 1 1 1 ] * 0.75, 'EdgeColor', [ 1 1 1 ] * 0.75 )
end
alines( lag, 'x', 'color', [ 1 0 0 ], 'linestyle', '-' );
ylim( yy )

return

% EOF