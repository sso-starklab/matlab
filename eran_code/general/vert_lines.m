% VERT_LINES    plot vertical lines.

% jan-03 ES

function lh = vert_lines(x,y,linelength,linecolor,linewidth)
lh = [];
if isempty(x), return, end
nx = length(x);
if length(y)==1
    Y = repmat([[y + [-1 1]*linelength/2] NaN]',nx,1);
elseif nx==length(y)
    Y = [y(:).'-linelength/2; y(:).'+linelength/2; NaN*ones(1,nx)];
else
    warning('input mismatch'), return
end
X = [x(:).'; x(:).'; NaN*ones(1,nx)];
lh = line(X(:),Y(:));
set(lh,'color',linecolor,'LineWidth',linewidth);
return
