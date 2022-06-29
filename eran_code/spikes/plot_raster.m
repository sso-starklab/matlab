% plot_raster           raster display for spike trains
%
% call                  [ sh, eh ] = plot_raster( st, timevec, ilevel, spike_length, spike_width, spike_color, et, ecolors, event_width )
%
% gets                  st              spike times. format is either:
%                                           1. columns of a full/sparse matrix (0,1)
%                                           2. vectors in a cell array (spike times)
%
% optional arguments (ordered list):
%
%                       timevec         {[]}        vector of all possible spike times (if st is sparse, then same length as number of rows in st)      
%                       ilevel          {1}         level (in axes) of first trial. may be either
%                                                   - an integer scalar (e.g. 5)
%                                                   - a vector of NT integers
%                       spike_length    {0.5}       "half height" of a pulse. if a character (e.g. '*', '.'), plots that character
%                       spike_width     {0.5}       width (thickness) of each tick
%                       spike_color     {[0 0 0]}   color of spike ticks
%                       et              {[]}        event times. each row - a different event; each column - different trial
%                       ecolors         {[1 0 0]}   color of event ticks; each row - a different event
%                       event_width     {spike_width}
%
% returns               sh              handle to spikes
%                       eh              handle to events
%
% calls                 vert_lines

% 18-mar-11 ES

% revisions
% 05-apr-11 correction of ytick scaling
% 02-feb-12 correction of cell array st, empty timevec case
% 24-aug-12 correction of cell array st, empty x case
% 18-feb-15 added symbol (e.g. '.' ) support in case st is cell array
% 06-dec-18 (1) modified event plotting to max between SPIKE_LENGTH and 1
%           (2) added argument EVENT_WIDTH to default to SPIKE_WIDTH

function [ sh, eh ] = plot_raster( st, timevec, ilevel, SPIKE_LENGTH, SPIKE_WIDTH, SPIKE_COLOR, et, ECOLORS, EVENT_WIDTH )

% arguments
nargs                   = nargin;
if nargs < 2 || isempty( timevec )
    timevec             = [];
end
if nargs < 3 || isempty( ilevel )
    ilevel              = 1;
end
if nargs < 4 || isempty( SPIKE_LENGTH )
    SPIKE_LENGTH        = 0.5;
end
if nargs < 5 || isempty( SPIKE_WIDTH )
    SPIKE_WIDTH         = 0.5;
end
if nargs < 6 || isempty( SPIKE_COLOR )
    SPIKE_COLOR         = [ 0 0 0 ];
end
if nargs < 7 || isempty( et )
    et                  = [];
end
if nargs < 8 || isempty( ECOLORS)
    ECOLORS             = [ 1 0 0 ];
end
if nargs < 9 || isempty( EVENT_WIDTH )
    EVENT_WIDTH         = SPIKE_WIDTH;
end

% preprocess
if iscell( st )
    [ ~, nt ]           = size( st );
elseif isa(st,'numeric') || isa( st, 'logical' )
    if isa( st, 'logical' )
        st              = double( st );
    end
    if ~issparse(st)
        st              = sparse( st );
    end
    [ win, nt ]         = size(st);
    if ~exist( 'timevec', 'var' ) || length( timevec ) ~= win
        timevec         = [];
    end
end
if numel(ilevel) ~= 1
    ilevel              = 1;
end
if isa( SPIKE_LENGTH, 'numeric' )
    d                   = SPIKE_LENGTH;
else
    d                   = 0.5;
end
LEVEL                   = ilevel : 2 * d : ( (nt -1 ) * 2 * d + ilevel );
ne                      = size( et, 1 );
if ne
    if size( et, 2 ) ~= nt
        fprintf( 1, 'should have an event for each trial\n' )
        et              = [];
    end
    if size( ECOLORS, 1 ) ~= ne
        ECOLORS         = repmat( ECOLORS( 1, : ), [ ne 1 ] );
    end
end

% plot spikes
sh                      = [];
eh                      = [];
hf0                     = ishold;
if ~hf0
    hold on
end
if nargs >= 1 && ~isempty( st )
    if issparse( st )
        [ x, y ]        = find( st );
        if isempty( timevec )
            t           = x;
        else
            t           = timevec( x );
        end
        if ~isempty( x )
            for i = 1 : nt
                levels( y == i ) = LEVEL( i );
            end
            if isa( SPIKE_LENGTH, 'char' )
                sh      = plot( t, levels, SPIKE_LENGTH );
            else
                sh      = vert_lines( t, levels, SPIKE_LENGTH, SPIKE_COLOR, SPIKE_WIDTH );
            end
            set( sh, 'color', SPIKE_COLOR );
        end
    else
        for trial       = 1 : nt
            x           = st{ trial };
            if isempty( x )
                continue
            end
            if isempty( timevec )
                t       = x;
            elseif x > length( timevec ) 
                continue
            else
                t       = timevec( x );
            end
            if isa( SPIKE_LENGTH, 'char' )
                sh( trial ) = plot( t, LEVEL( trial ) * ones( size( t ) ), SPIKE_LENGTH );
                set( sh( trial ), 'color', SPIKE_COLOR )
            else
                sh( trial ) = vert_lines(...
                    t, LEVEL( trial ), SPIKE_LENGTH, SPIKE_COLOR, SPIKE_WIDTH );
            end
        end
    end
end

% plot events
for ei = 1 : ne
    eh( ei )            = vert_lines( et( ei, : ), LEVEL, max( SPIKE_LENGTH, 1 ), ECOLORS( ei, : ), EVENT_WIDTH );
end

ylim( LEVEL( [ 1 end ] ) + [ -1 1 ] )
if ~isempty( timevec ) && length( unique( timevec ) ) > 1
    xlim( timevec( [ 1 end ] ) )
end
if hf0
    hold off
end

return

% EOF
