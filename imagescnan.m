% imagescnan        image in myjet colormap on current plot, NaNs in white
%
% call              [ h, ah ] = imagescnan( x )
%                   imagescnan( x, clim )
%                   imagescnan( xls, yls, x, clim, func )
%
% gets              arguments of imagesc + a colormap function
%
% returns           h               handle to image
%                   ah              handle to axes

% 04-aug-21 ES

function [ h, ah ]  = imagescnan( xls, yls, x, clim, map )

nargs                           = nargin;
switch nargs
    case 1
        x                       = xls;
        xls                     = [];
        yls                     = [];
        clim                    = [];
    case 2
        x                       = xls;
        clim                    = yls;
        xls                     = [];
        yls                     = [];
    case 3
        clim                    = [];
end
if nargs < 5 || isempty( func )
    map                         = myjet;
end

newplot
ah                              = gca;
if isempty( clim )
    h                           = imagesc( xls, yls, x, 'CDataMapping', 'scaled' );
else
    h                           = imagesc( xls, yls, x, clim,  'CDataMapping', 'scaled' );
end
set( h, 'AlphaData', ~isnan( x ) );

set( ah, 'tickdir', 'out', 'box', 'off' );
colormap( ah, map )
axis xy

return

% EOF