% stim_get_channels          a new structure with the given channels
%
% out = stim_get_channels( stim, channels, eq )
%
% stim          has to conform to stim_check
% channels      a list of channels to include
% eq            operation (relevant for simultaneous structures only):
%                   {'le'} indicates any subset of channels
%                   'eq' requires an exact match (out of the available channels)
%                   'ge' allows additional channels to be active at the same time 
%
% see also          stim_make, stim_check, stim_get
%
% for example, calling 
%   >> out = stim_get_channels( stim, [ 33 35 ], cmp );
% with a simultaneous strucutre while setting cmp to 
%   'le', will include any sim events in which either 33 or 35 participated;
%   'ge', will include events in which both 33 and 35 (and possibly other) participated; and
%   'eq', will include only events in which exactly 33 and 35 participated
%
% see also          stim_select

% 25-feb-13 ES

function out = stim_get_channels( stim, channels, cmp )

out = [];
nargs = nargin;
if nargs < 2 || isempty( stim ) || ~stim_check( stim ) || length( stim ) ~= 1
    return
end
if nargs < 2 || isempty( channels )
    return
end
if nargs < 3 || isempty( cmp )
    cmp = 'le'; % le
end
if strcmp( cmp, 'eq' ) || strcmp( cmp, 'ge' ) || strcmp( cmp, 'le' )...
        || isequal( cmp, @eq ) || isequal( cmp, @ge ) || isequal( cmp, @le )...
        || strcmp( cmp, '==' ) || strcmp( cmp, '>=' ) || strcmp( cmp, '<=' )
    switch cmp
        case '==', cmp = 'eq';
        case '>=', cmp = 'ge';
        case '<=', cmp = 'le';
    end
else
    return
end
if isempty( stim.types )
    return
end

cidx = ismember( stim.chan, channels );
if sum( cidx ) == 0
    return
elseif length( stim.chan ) > 1
    nmatches = sum( stim.index( :, cidx ) > 0, 2 );
    nchans = sum( stim.index > 0, 2 );
    ridx = feval( cmp, nmatches, sum( cidx ) ) & nmatches >= 1;
    if strcmp( cmp, 'eq' ) || isequal( cmp, @eq )
        ridx = ridx & nchans == length( channels );
    end
    if sum( ridx ) > 0
        stim = stim_get( stim, ridx );
        stim.chan = stim.chan( cidx );
        stim.source = stim.source( cidx );
        stim.voltagerange = stim.voltagerange( cidx );
        stim.median = stim.median( cidx );
        stim.index = stim.index( :, cidx );
        out = stim;
    end
else
    out = stim;
end

return

% EOF