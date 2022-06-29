% makegaussfir          for 1D smoothing
%
% gwin = makegaussfir( sd, Fs )
% 
% sd        Gaussian SD [sec]
% Fs        sampling rate [Hz]

% 26-dec-12 ES

function gwin = makegaussfir( sd, Fs, support )

n = 2 * floor( sd * Fs / 2 ) + 1; % make odd
if ~exist( 'support', 'var' )
    support = 6;
end
switch length( sd )
    case 1
        win = gausskernel( n, 0, support * n( 1 ) + 1, 1 ); % make w/ support of 3 SDs
    case 2
        win = gausskernel( n( 1 ), n( 2 ), support * n( 1 ) + 1, support * n( 2 ) + 1 ); % make w/ support of 3 SDs
end
gwin = win / sum( win( : ) );

return

% EOF

% for uncorrelated noise, boxcar filtering reduces the variance by sqrt( n )
% and gaussian filtering by sqrt( n ) * 1.96 (actually slightly less reduction d.t. edge effects)
nw = 5; r = randn( 1e6, 1 ); 
w = ones( nw, 1 ) / nw; rf = firfilt( r, w ); 
g = makegaussfir( nw, 1, 18 ); rg = firfilt( r, g ); 
[ std( r ), std( r ) / sqrt( nw ) std( rf ) std( rg ) std( rg ) * 1.96  ]
