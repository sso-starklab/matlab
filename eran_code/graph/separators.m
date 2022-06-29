% SEPARATORS        plots seperators on current axis.
%
% call              LH = SEPARATORS( X )
%                   SEPARATORS( ..., AH, COLOR, MODE )
%
% gets              X               points to place separators
%                   AH              axis handle {gca}
%                   LINECOLOR       3 point vector {[ 1 0 0 ]}
%                   MODE            which axis to separate {'x'}
%                   LINEWIDTH       width {0.5}
%                   LINESTYLE       style {'-'}
%
% returns           LH      handle to the lines
%
% calls             nothing

% stability package
% Jan-2003 ES

% revisions
% 16-dec-03 lh returned
% 02-may-04 width and style added
% 11-sep-19 cleaned up

function lh = separators( x, ah, linecolor, mode, linewidth, linestyle )

nargs                   = nargin;
if nargs < 2 || isempty( ah )
    ah                  = gca; 
end
if nargs < 3 || isempty( linecolor )
    linecolor           = [ 1 0 0 ]; 
end
if nargs < 4 || isempty( mode )
    mode                = 'x'; 
end
if nargs < 5 || isempty( linewidth )
    linewidth           = 0.5; 
end
if nargs < 6 || isempty( linestyle )
    linestyle           = '-'; 
end

x                       = x(:);
[ nx, cols ]            = size(x);
if nx > 1 && cols > 1
    warning('input should be a vector'); 
end
YLIM                    = get(ah,'ylim');
XLIM                    = get(ah,'xlim');
if lower(mode) == 'x'
    X                   = [ x x NaN * ones( nx, 1 ) ]';
    X                   = X( : );
    Y                   = repmat( [ YLIM NaN ]', nx, 1 );
elseif lower(mode) == 'y'
    XLIM                = get( ah, 'xlim' );
    Y                   = [ x x NaN * ones( nx, 1 ) ]';
    Y                   = Y(:);
    X                   = repmat([XLIM NaN]',nx,1);
else
    warning('incorrect mode')
end
subplot(ah)
lh                      = line(X,Y);
set( lh, 'color', linecolor, 'LineWidth', linewidth, 'LineStyle', linestyle );
set( ah, 'ylim', YLIM, 'xlim', XLIM );

return

% EOF
