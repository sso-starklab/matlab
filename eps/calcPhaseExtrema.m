% calcPhaseExtrema              based on troughs and peaks 
%
% call                          [ phs, xf, pidxs ] = calcPhaseExtrema( x, BP, Fs, dt, dph, graphics )
%
% gets                          x           signal (vector or matrix)
%                               BP          {[ 1 60 ]}, [Hz] bandpass
%                               Fs          {1250}, [Hz], sampling rate
%                               dt          {Fs/20}, [samples], pruning
%                               dph         {2*pi/2}, [rad], smoothing
%                               graphics    {0}, plots i'th column
%
% returns                       phs         phases, [rad] (0 to 2*pi range)
%
% calls                         firfilt, local_max, makegausslpfir
%
% does                          (1) computes band-pass signal (2-pole Butterworth)
%                               (2) detects all peaks and troughs in the bandpass
%                               (3) prunes to dt between peaks, troughs
%                               (4) assigns 0 to peak, pi to troughs, and interpolates linearly between them
%                               (5) smooths in the time domain
%                               (6) wraps back to 0 2*pi range
%
% everything is done column-wise
%
% see also                      calcPhase

% 31-mar-21 ES

% revisions
% 10-apr-21 added intermediate products to output arguments
% 11-apr-21 offset cumsum by pi
% 29-jun-21 added smoothing of phase
% 04-Aug-21 updated help

function [ phs, xf, pidxs, mT, Fc ] = calcPhaseExtrema( x, BP, Fs, dt, dph, graphics )

% constants
npoles                          = 2;
imethod                         = 'linear';

% arguments
nargs                           = nargin;
if nargs < 1 || isempty( x )
    return
end
if nargs < 2 || isempty( BP )
    BP                          = [ 1 60 ];                                 % [Hz]
end
if nargs < 3 || isempty( Fs )
    Fs                          = 1250;                                     % [Hz]
end
if nargs < 4 || isempty( dt )
    dt                          = Fs / 20;                                  % [samples]
end
if nargs < 5 || isempty( dph )
    dph                         = ( 2 * pi ) / 2;                           % [rad]
end
if nargs < 5 || isempty( graphics )
    graphics                    = 0;
end


% build wideband filter
[ bwb, awb]                     = butter( npoles, BP / Fs * 2, 'bandpass' );

% filter wideband
[ m, n ]                        = size( x );
xf                              = filtfilt( bwb, awb, x );

% detect extrema, spaced by dt samples
[ pidxs, ~, etypes ]            = local_max( xf, 'ext', dt );
if n == 1
    pidxs                       = [ pidxs ones( length( pidxs ), 1 ) ];
end

% interpolate
phs                             = NaN( m, n );
phs_ns                          = NaN( m, n );


for i                           = 1 : n
    
    % column-specific values
    idx                         = pidxs( :, 2 ) == i;
    pidx                        = pidxs( idx, 1 );
    etype                       = etypes( idx );
    
    % interpolation grid
    npeaks                  	= length( etype );
    if etype( 1 ) == 1
        bias                    = 0;
    else
        bias                    = pi;
    end
    v                           = cumsum( [ 0; pi * ones( npeaks - 1, 1 ) ] );
    ti                          = ( pidx( 1 ) : pidx( npeaks ) )';
    
    % actually interpolate
    phi                         = interp1( pidx, v, ti, imethod );
    phi_ns                      = phi; %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % smooth
    if ~isnan( dph )
        mT                      = mean( diff( pidx( etype == 1 ) ) );       % cycle duration [samples]
        Fc                      = Fs / ( mT * dph / ( 2 * pi ) );           % cutoff frequency [Hz]
        win                     = makegausslpfir( Fc, Fs );
        phi                     = firfilt( phi, win );
    end
    
    % remove bias, modulus to [ 0 2*pi] range, pad edges, assign
    phi                         = mod( phi - bias, 2 * pi );
    prepad                      = NaN( pidx( 1 ) - 1, 1 );
    postpad                     = NaN( length( xf ) - pidx( npeaks ), 1 );
    phs( :, i )                 = [ prepad; phi; postpad ];
    
    phi_ns                      = mod( phi_ns - bias, 2 * pi ); %%%%%%%%%%%%%%%%%%%%%%%%%%%
    phs_ns( :, i )              = [ prepad; phi_ns; postpad ]; %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% output
if nargout >= 3
    pidxs( :, 3 )               = etypes;
end

% graphics
if ~graphics
    return
end

newplot
if graphics >= 1 && graphics < n
    i                           = round( n );
else
    i                           = 1;
end
tt                              = ( 1 : m ) / Fs;
maxx                            = max( x( :, i ) );
sphs                            = phs( :, i ) * maxx / ( 2 * pi ) - maxx / 2;
ph                              = plot( tt, x( :, i ), 'b', tt, xf( :, i ), 'r', tt, sphs, 'k' );
set( ph( 1 ), 'color', [ 1 1 1 ] * 0.5 );
hold on;
sphs_ns                         = phs_ns( :, i ) * maxx / ( 2 * pi ) - maxx / 2; %%%%%%%
plot( tt, sphs_ns, 'b' ); %%%%%%%%%
xlim( tt( [ 1 end ] ) )
ylim( [ min( x( :, i ) ) max( x( :, i ) ) ] )
idx                             = pidxs( :, 2 ) == i;
pks                             = pidxs( idx, 1 );
typs                            = etypes( idx );
alines( pks( typs == 1 ) / Fs, 'x', 'color', [ 1 0 0 ], 'linestyle', '-' );
alines( pks( typs == -1 ) / Fs, 'x', 'color', [ 1 0 0 ], 'linestyle', '--' );
set( gca, 'tickdir', 'out', 'box', 'off' )
ylabel( 'Amplitude/phase' )
xlabel( 'Time [s]' )

return

% EOF

