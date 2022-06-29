% stim_get          a new structure with the given indices/type
%
% out = stim_get( stim, sidx )
%
% stim      has to conform to stim_check
% sidx      can be either a vector of indices, a char type, or a cell array
%               of types
% 
% see also          stim_make, stim_check
%
% example
%   out = stim_get( stim, 'ZAP' );
%
% note
% in principle should support retrieval by any field/combination of
% fields... this can evolve to a full DB functionality... so we keep it
% simple right now, either by type or by indices
%
% calls             nothing
%
% see also          stim_select (enables higher level functionality)
%                   stim_get_channels (enables selection by chan number)

% 18-jan-13

% revisions
% 28-jan-13 extended to support char arguments


function out = stim_get( stim, sidx )

out                 = stim;

% determine input type
if ~stim_check( stim ) || length( stim ) ~= 1
    return
end
n                   = size( stim.times, 1 );
ok                  = true( 1 );
nidx                = length( sidx );
if isa( sidx, 'logical' ) && nidx == n
elseif isa( sidx, 'numeric' ) && max( sidx ) <= n
elseif isa( sidx, 'char' )
    sidx            = ismember( upper( stim.types ), upper( sidx ) );
elseif isa( sidx, 'cell' )
    for i           = 1 : length( sidx )
        if ~isa( sidx{ i }, 'char' )
            ok      = 0;
            break;
        end
    end
    if ok
        sidx        = ismember( upper( stim.types ), upper( sidx ) );
    end
else
    ok              = 0;
end
if ~ok
    return
end

% fill output
if size( stim.index, 1 ) == n
    out.index       = stim.index( sidx, : );
else
    out.index       = [];
end
if size( stim.types, 1 ) == n
    out.types       = stim.types( sidx, : );
else
    out.types       = {};
end
if size( stim.category, 1 ) == n
    out.category    = stim.category( sidx, : );
else
    out.category    = [];
end
if size( stim.durs, 1 ) == n
    out.durs        = stim.durs( sidx, : );
else
    out.durs        = [];
end
if size( stim.times, 1 ) == n
    out.times       = stim.times( sidx, : );
else
    out.times       = [];
end
if size( stim.slopes, 1 ) == n
    out.slopes      = stim.slopes( sidx, : );
else
    out.slopes      = [];
end
if size( stim.vals, 1 ) == n
    out.vals        = stim.vals( sidx, : );
else
    out.vals        = [];
end
if size( stim.stats, 1 ) == n
    out.stats       = stim.stats( sidx, : );
else
    out.stats       = [];
end
if size( stim.franges, 1 ) == n
    out.franges     = stim.franges( sidx, : );
else
    out.franges     = [];
end
if size( stim.plateaus, 1 ) == n
    out.plateaus    = stim.plateaus( sidx, : );
else
    out.plateaus    = [];
end

return

% EOF
