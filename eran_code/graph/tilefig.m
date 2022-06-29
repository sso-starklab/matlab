% tilefig           into subplots and sub-subplot
%
% [ ah, fig ] = tilefig( nrows, ncols, nsubs, f, align )
%
% nrows, ncols      any number is supported
% nsubs             1, 2 (right/left), -2 (top/bottom), and 4 ([1 3;2 4])
%                       are supported; default {1}
% f                 {0.85}; covered fraction of fig
% align             {'edge'} aligns to the bottom right. 
%                       Other options: 'topleft', 'center'
% 
% ah                handle to subplots; m x nsubs, where m = nrows x ncols
% fh                handle to figure

% 14-mar-13 ES

function [ ah, fh ] = tilefig( nrows, ncols, nsubs, f, align )

nargs = nargin;
if nargs < 1 || isempty( nrows )
    nrows = 1;
end
if nargs < 2 || isempty( ncols )
    ncols = 1;
end
if nargs < 3 || isempty( nsubs )
    nsubs = 1;
end
if ~ismember( nsubs, [ 1 -2 2 4 ] )
    error( 'not supported' )
end
if nargs < 4 || isempty( f )
    f = 0.85;                            % determines the covered fig portion
end
if ~isnumeric( f )
    error( 'erroneous format' )
end
f = f( 1 );
if nargs < 5 || isempty( align )
    align = '';
end

align = lower( align );
pp = nrows * ncols;
if isempty( align )
    if pp <= 9
        align = 'center';
    else
        align = 'edge';
    end
end

wx = 1 / ncols * f;
wy = 1 / nrows * f;
switch align
    case 'center'
        x0 = ( 1 - f ) / 2 / ncols : 1 / ncols : ( 1 - ( 1 - f ) / 2 / ncols ); % align center
        y0 = ( 1 - f ) / 2 / nrows : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align center
    case 'right'
        x0 = fliplr( 1 : -1 / ncols : ( 1 - f ) / 2 / ncols ) - wx; % align right
        y0 = ( 1 - f ) / 2 / nrows : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align center
    case { 'edge', 'bottomright' }
        x0 = fliplr( 1 : -1 / ncols : ( 1 - f ) / 2 / ncols ) - wx; % align right
        y0 = 0 : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align bottom
    case { 'bottomleft' }
        x0 = 0 : 1 / ncols : ( 1 - ( 1 - f ) / 2 / ncols ); % align left
        y0 = 0 : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align bottom
    case { 'topleft' }
        x0 = 0 : 1 / ncols : ( 1 - ( 1 - f ) / 2 / ncols ); % align left
        y0 = fliplr( 1 : -1 / nrows : ( 1 - f ) / 2 / nrows ) - wy; % align top
    case { 'topright' }
        x0 = fliplr( 1 : -1 / ncols : ( 1 - f ) / 2 / ncols ) - wx; % align right
        y0 = fliplr( 1 : -1 / nrows : ( 1 - f ) / 2 / nrows ) - wy; % align top
    case { 'centerright' }
        x0 = fliplr( 1 : -1 / ncols : ( 1 - f ) / 2 / ncols ) - wx; % align right
        y0 = ( 1 - f ) / 2 / nrows : 1 / nrows : ( 1 - ( 1 - f ) / 2 / nrows ); % align center
        
end
x0 = ones( nrows, 1 ) * x0( 1 : ncols );
y0 = flipud( y0( 1 : nrows )' * ones( 1, ncols ) );
wx = wx * ones( nrows, ncols );
wy = wy * ones( nrows, ncols );
x0 = x0';
y0 = y0';
wx = wx';
wy = wy';
ah = zeros( pp, abs( nsubs ) );

fh = figure;
for i = 1 : pp,
    pos = [ x0( i ) y0( i ) wx( i ) wy( i ) ];
    switch nsubs
        case 1
            spos( 1, : ) = pos;
        case 2 % split this into two: [ 1 2 ]
            spos( 1, : ) = [ pos( 1 )                   pos( 2 )                    pos( 3 ) / 2    pos( 4 ) ];
            spos( 2, : ) = [ pos( 1 ) + pos( 3 ) / 2    pos( 2 )                    pos( 3 ) / 2    pos( 4 ) ];
        case -2 % split this into two: [ 1; 2 ]
            spos( 1, : ) = [ pos( 1 )                   pos( 2 ) + pos( 4 ) / 2     pos( 3 )        pos( 4 ) / 2 ];
            spos( 2, : ) = [ pos( 1 )                   pos( 2 )                    pos( 3 )        pos( 4 ) / 2 ];
        case 4 % split this into four: [ 1 3; 2 4 ]
            spos( 1, : ) = [ pos( 1 )                   pos( 2 ) + pos( 4 ) / 2     pos( 3 ) / 2    pos( 4 ) / 2 ];
            spos( 2, : ) = [ pos( 1 )                   pos( 2 )                    pos( 3 ) / 2    pos( 4 ) / 2 ];
            spos( 3, : ) = [ pos( 1 ) + pos( 3 ) / 2    pos( 2 ) + pos( 4 ) / 2     pos( 3 ) / 2    pos( 4 ) / 2 ];
            spos( 4, : ) = [ pos( 1 ) + pos( 3 ) / 2    pos( 2 )                    pos( 3 ) / 2    pos( 4 ) / 2 ];
    end
    for k = 1 : abs( nsubs )
        ah( i, k ) = axes( 'position', spos( k, : ) );
    end
end

return

% EOF
