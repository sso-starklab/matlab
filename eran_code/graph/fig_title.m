% FIG_TITLE             title for a figure.
%
% call                  th = fig_title( tstr, f, pos, TFS )
%
% gets                  tstr        text for the title
%                       f           { gcf }                         handle to figure
%                       pos         { [ 0.5 0.925 0.01 0.01 ] }     position of text
%                       TFS         { 8 }                           title font size
%
% 

% 22-dec-02 ES

% revisions
% 01-jan-02 revised
% 08-mar-03 f
% 02-apr-03 pos
% 12-dec-03 TFS
% 14-sep-19 cleaned up

function th = fig_title( tstr, f, pos, TFS )

nargs = nargin;
if nargs < 2 || isempty( f )
    f                   = gcf;
end
if nargs < 3 || isempty(pos)
    pos                 = [ 0.5 0.925 0.01 0.01 ]; 
end
if nargs < 4 || isempty( TFS )
    TFS                 = 8; 
end

th                      = textf( pos( 1 ), pos( 2 ), tstr, 'Fontsize', TFS );
set( th, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom' )

return

% EOF

