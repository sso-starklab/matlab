% PlotStims             see GetTrigVals
%
% CALL                  ph = PlotStims( dur, amp, trig, durscale, ampscale, durlabel, amplabel, jitFac )
%
% GETS                  graphics  {0}
%                       durscale  {1/20000}
%                       ampscale  {1}
%                       durlabel  {'Time [s]'}
%                       amplabel  {'I [A]'}
%                       jitFac    {0.01}
%
% CALLS                 myjet
%
% see also LoadVals

% 12-jan-13 ES

% revisions
% 30-aug-18 changed call to legend to allow compatibility with R2018a
% 18-aug-19 cleaned up

function ph = PlotStims( dur, amp, trig, durscale, ampscale, durlabel, amplabel, jitFac )

% initialize output
ph                          = [];

% arguments
nargs = nargin;
if nargs < 3 || isempty( dur ) || isempty( amp ) || isempty( trig )
    return
end
if nargs < 4 || isempty( durscale )
    durscale                = 1;                        % samples->s (or any other time display units)
end
if nargs < 5 || isempty( ampscale )
    ampscale                = 1;                        % V->A (or any other amp display units)
end
if nargs < 6 || isempty( durlabel )
    durlabel                = 'Time [s]';
end
if nargs < 7 || isempty( amplabel )
    amplabel                = 'I [A]';
end
if nargs < 8 || isempty( jitFac )
    jitFac                  = 0.01;
end
if jitFac < 0 || jitFac > 1
    jitFac                  = 0.01;
end

% plot on current axes
trigchans                   = unique( trig( : ).' );
ntrigs                      = length( trigchans );
if ntrigs == 0
    return
end
map                         = myjet;
nmap                        = size( map, 1 );
if nmap < ntrigs
    map                     = repmat( map, [ ceil( ntrigs / nmap ) 1 ] );
    nmap                    = size( map, 1 );
end
colors                      = map( round( 1 : nmap / ntrigs : nmap ), : );
ph                          = zeros( size( trigchans ) );
legstr                      = cell( size( trigchans ) );
if ~ishold
    newplot
end
for i                       = 1 : length( trigchans )
    trigchan                = trigchans( i ); 
    idx                     = ismember( trig, trigchan );
    jit                     = 1 + ( rand( sum( idx ), 1 ) - 0.5 ) * jitFac; 
    ph( i )                 = plot( dur( idx ) .* jit * durscale, amp( idx ) * ampscale, '.' ); 
    set( ph( i ), 'color', colors( i, : ) );
    if isnumeric( trigchan )
        str                 = num2str( trigchan );
    elseif isa( trigchan, 'cell' ) && isa( trigchan{ 1 }, 'char' )
        str                 = trigchan{ 1 };
    else
        str                 = num2str( i );
    end
    legstr{ i }             = str;
    hold on; 
end
xlabel( durlabel )
ylabel( amplabel )
legend( ph, legstr )
grid on
set( gca, 'tickdir', 'out', 'box', 'off' )

return

% EOF
