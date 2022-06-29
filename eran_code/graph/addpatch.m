% ADDPATCH      adds a rectangular patch to the current figure
%
% CALL          ph = addpatch( xx, yy, c, aa )
%
% CALLS         mypacth

% 08-may-11 ES

% revisions
% 18-aug-19 cleaned up

function ph = addpatch( xx, yy, c, aa )

% arguments
nargs                       = nargin;
if nargs < 3 || isempty( c )
    c = [ 0 0 1 ]; 
end
if nargs < 4 || isempty( aa )
    aa = 1; 
end

% preparations
xx                          = sort( xx );
yy                          = sort( yy );
xx                          = xx( 1 : 2 );
yy                          = yy( 1 : 2 );
c                           = c( 1 : 3 );
c( c < 0 )                  = 0;
c( c > 1 )                  = 1;
xx                          = [ xx xx( 2 ) xx( 1 ) ];
yy                          = [ yy( 1 ) yy( 1 ) yy( 2 ) yy( 2 ) ];

% plot
holdstate = ishold;
if ~holdstate
    hold on
end
xlims                       = xlim;
ylims                       = ylim;
ph                          = mypatch( xx, yy, c );
if ~holdstate
    hold off
end
set( gca, 'xlim', xlims, 'ylim', ylims );
if aa < 1 && aa > 0
    set( ph, 'EdgeAlpha', aa, 'FaceAlpha', aa )
end

return

% EOF
