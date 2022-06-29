% parseSimEvents            parse simultaneous events from stims structure (multiple channels)
%   
% call                      [ STIM, STIMS ] = PARSESIMEVENTS( STIMS, CHANNELS, SPKFS )
%
% gets                      stims           vector of stim structures (see stim_make)
%                           channels        (optional), list of channels
%                           spkFs           (optional), otherwise from par
%
% returns                   stim            structure for simultaneous events
%                           stims           input with index field updated
%
% calls                     isoverlap, prunemat, stim_make, stim_get, LoadXml
%
% called by parseNchannels

% 01-aug-19 ES based on parseMultipleChannels

function [ stim, stims ] = parseSimEvents( stims, channels, spkFs )

stim            = stim_make;
mfname          = upper( mfilename );

nargs = nargin;
if nargs < 1 || isempty( stims )
    return
end
if nargs < 2 || isempty( channels )
    channels    = unique( [ stims.chan ] );
end
if nargs < 3 || isempty( spkFs )
    try 
        par     = LoadXml( stims( 1 ).filebase ); 
        spkFs   = par.SampleRate;
    catch
        spkFs   = 20000;
    end
end
nchans          = length( channels );

%------------------------------------------------------------------%
% detect simultaneous events, build a common structure, and update
% scheme is to check, for each channel, all other channels for overlap
%------------------------------------------------------------------%

fprintf( 1, '%s: Detecting simultaneous events...\n', mfname )
% do a first pass over the channels to know stimulus types
utypes          = [];
for i           = 1 : nchans
    utypes      = [ utypes; stims( i ).types ];
end
utypes          = upper( unique( utypes ) );
ntypes          = length( utypes );
x               = cell( 1, ntypes );
for k           = 1 : ntypes
    x{ k }      = zeros( 0, nchans );
end
y                   = x;
for i               = 1 : ( nchans - 1 )
    for k           = 1 : ntypes
        sidx0       = ismember( stims( i ).types, utypes{ k } );
        mat0        = stims( i ).times( sidx0, : );
        for j       = ( i + 1 ) : nchans
            sidx1               = ismember( stims( j ).types, utypes{ k } );
            mat1                = stims( j ).times( sidx1, : );
            [ idx0, idx1 ]      = isoverlap( mat0, mat1 );
            nidx = sum( idx0 );
            if nidx
                fprintf( 1, '\t\tdetected %d simultaneous %s events (%d x %d)\n'...
                    , nidx, utypes{ k }, i, j )
                block           = zeros( nidx, nchans );
                block( :, i )   = find( idx0 );
                block( :, j )   = idx1( idx0 );
                x{ k }          = [ x{ k }; block ];
                % the times:
                % [ mat0( idx0, : )  mat1( idx1( idx0 ), : ) ]
            end
        end
    end
end

% translate indices from relative to absolute
for i                   = 1 : nchans
    for k               = 1 : ntypes
        sidx            = find( ismember( stims( i ).types, utypes{ k } ) );
        idx             = x{ k }( :, i );
        idx( idx ~= 0 ) = sidx( idx( idx ~= 0 ) );
        x{ k }( :, i )  = idx;
    end
end

% prune pair-wise matrices
fprintf( 1, '%s: Pruning pair-wise matrices...\n', mfname )
for k           = 1 : ntypes
    mat0        = x{ k };
    if size( mat0, 1 ) == 0
        continue
    end
    y{ k }      = prunemat( mat0 );
end

% build the common structure
fprintf( 1, '%s: Building common structure...\n', mfname )
stim                = stim_make;
stim.filebase       = stims( 1 ).filebase;
stim.suffix         = stims( 1 ).suffix;
stim.duration       = stims( 1 ).duration;
stim.generator      = stims( 1 ).generator;
stim.chan           = channels;
for i                       = 1 : nchans
    stim.source{ i }        = stims( i ).source;
    stim.voltagerange( i )  = stims( i ).voltagerange;
    stim.median( i )        = stims( i ).median;
end

for k               = 1 : ntypes
    mat             = y{ k };
    if size( mat, 1 ) == 0
        continue
    end
    stim.index      = [ stim.index; mat ];
    for j           = 1 : size( mat, 1 )
        idx         = mat( j, : );
        cidx        = find( idx > 0 );
        idx         = idx( cidx );
        % collect
        nidx        = length( idx );
        types       = cell( nidx, 1 );
        times       = zeros( nidx, 2 );
        slopes      = zeros( nidx, 2 );
        vals        = zeros( nidx, 1 );
        stats       = zeros( nidx, 2 );
        franges     = zeros( nidx, 2 );
        plateaus    = zeros( nidx, 1 );
        h           = 0;
        for i = cidx
            h = h + 1;
            types{ h }          = stims( i ).types{ idx( h ) };
            times( h, : )       = stims( i ).times( idx( h ), : );
            slopes( h, : )      = stims( i ).slopes( idx( h ), : );
            vals( h, : )        = stims( i ).vals( idx( h ), : );
            stats( h, : )       = stims( i ).stats( idx( h ), : );
            franges( h, : )     = stims( i ).franges( idx( h ), : );
            plateaus( h, : )    = stims( i ).plateaus( idx( h ), : );
        end
        % average
        types                   = unique( types );
        [ minval, minidx ]      = min( times( :, 1 ) );
        [ maxval, maxidx ]      = max( times( :, 2 ) );
        times                   = [ minval maxval ];
        durs                    = ( diff( times ) + 1 ) / spkFs; % ( diff( round( times / spkFs * Fs ) ) + 1 ) / Fs;
        slopes                  = [ slopes( minidx, 1 ) slopes( maxidx, 2 ) ];
        vals                    = mean( vals );
        stats                   = [ mean( stats( :, 1 ) ) std( stats( :, 1 ) ) ];
        franges                 = mean( franges );
        plateaus                = mean( plateaus );
        % populate
        stim.types              = [ stim.types; types ];
        stim.durs               = [ stim.durs; durs ];
        stim.times              = [ stim.times; times ];
        stim.slopes             = [ stim.slopes; slopes ];
        stim.vals               = [ stim.vals; vals ];
        stim.stats              = [ stim.stats; stats ];
        stim.franges            = [ stim.franges; franges ];
        stim.plateaus = [ stim.plateaus; plateaus ];
    end
    fprintf( 1, '\t\tupdated with %d simultaneous %s events\n'...
        , size( mat, 1 ), utypes{ k } )
end

fprintf( 1, '%s: Updating individual structures...\n', mfname )
if ~isempty( stim.times )
    % sort by onset time
    [ ~, sidx ]         = sort( stim.times( :, 1 ) );
    stim                = stim_get( stim, sidx );
    stim.category       = zeros( size( stim.times, 1 ), 3 );
else
    fprintf( 1, '%s: No simultaneous events detected!\n', mfname )
end
for i                   = 1 : nchans
    if isempty( stims( i ).durs )
        continue
    end
    stims( i ).index( : )           = 1; % parse1channel default is not unique (0), so first tag all as unique
    if ~isempty( stim.times )
        idx                         = stim.index( :, i );
        idx( idx == 0 )             = [];
        stims( i ).index( idx )     = 0;
    end
end

return

% EOF

%------------------------------------------------------------------%
% Notes:
%------------------------------------------------------------------%
% assumption 1: simultaneous events are of similar type
% to release this assumption - set sameType flag to zero (the default).
% then all stimuli are checked. if a stimulus is erroneously classified
% erroneously and sameType is 1
%
% assumption 2: any overlap is considered simultaneity. thus partially
% overlapping stimuli are considered simultaneous. to distinguish, the
% individual stim structures have to be examined
%
% algorithm:
%
% if sameFlag
%       initialize a matrix with the number of channels
% end
% for channels:
%   if sameFlag
%       do the for loop only once
%   else
%       determine the stimulus types in this channel
%   end
%   for all stimulus types
%       if ~sameFlag & a new stimulus
%           initialize a matrix
%       end
%       for all subsequent channels (not all, saves half the time...)
%           -check for overlap
%           -keep overlapping stimulus indices in a matrix: columns are channels
%           and rows are indices. indices are zero if a channel does not have a
%           an overlapping stimulus
% for all stimulus matrices
%   -prune rows (e.g. channels 1 4 6 have overlapping stimuli, these will be
%   detected twice - when 1>4,1>6 are checked and again when 4>6; in
%   principle overlapping stimuli could be excluded from the list to
%   prevent this behavior and speed things up)
%
%   -get stimulus parameters:
%   'times' is 2 col, from the onset of the earliest to the offset of the latest
%   'types' is the stimulus type. if sameType is 0, will be a cell array matrix
%   'slopes' is like time, the onset slope of the earliest and the offset slope of the latest
%   'vals' is the mean over channels
%   'stats' is again the mean (for the mean) and the SD (over the means,
%   i.e. it has a different meaning than that of the indvidual channels)
%   for the SDs the individual channels structs have to be accessed
%   'franges' is the mean of the ranges
%   'plateaus' is again the mean value

%------------------------------------------------------------------%
% end of notes
%------------------------------------------------------------------%