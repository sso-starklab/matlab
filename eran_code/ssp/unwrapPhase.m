% unwrapPhase           into discrete cycles
%
% call                  [ uphs, talign ] = unwrapPhase( phs, t0 )
% 
% gets                  phs         vector, phases 0-2*pi
%                       t0          center
%
% returns               uphs        unwrapped phases
%
% does
% -looks for the trough (phase pi) closest to t0. defines this as the
% trough of cycle0
% -goes forward through the data, define a peak as a 0-phase crossing
% preceeded by a pi-phase crossing (and vice versa for troughs)
% -do the same backwards
% -define cycles between peaks
%
% calls                 nothing
%
% note: 
% to solve wriggline, use monotonic
% 
% see also              calcPhase, monotonic

% 14-oct-13 ES

% revisions
% 23-mar-21 cleaned up

function [ uphs, talign ] = unwrapPhase( phs, t0 )

uphs                            = [];
talign                          = [];

nargs                           = nargin;
if nargs < 1 || isempty( phs )
    fprintf( '%s: phs missing\n', upper( mfilename )  )
    return
end
phs                             = phs( : );
n                               = length( phs );
if nargs < 2 || isempty( t0 )
    t0 = 1;
end
uphs                            = phs;
talign                          = t0;
if t0 > n || t0 < 1 || length( t0 ) > 1 || t0 ~= round( t0 ) 
    fprintf( '%s: t0 mismatch\n', upper( mfilename )  )
    return
end

% detect the cycle0 trough
zc                              = phs - pi;
trf                             = find( zc( 1 : n - 1 ) <= 0 & zc( 2 : n ) >= 0 );
if isempty( trf )
    return
end
[ ~, minidx ]                   = min( abs( t0 - trf ) );
trf0                            = trf( minidx );

% find first peak before it
uzi                             = unwrap( phs );
pk0                             = floor( uzi( trf0 ) / ( 2 * pi ) ) * 2 * pi;
uphs                            = uzi - pk0;
talign                          = t0 - trf0;                                % positive number - shift the trigger back

return

% EOF

