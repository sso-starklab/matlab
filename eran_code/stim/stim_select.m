% stim_select          a new structure with the given indices/type
%
% CALL                  out = stim_select( stim, fieldname1, fieldvals1, fieldname2, fieldvals2, ... )
%
% GETS
%                       stim              has to conform to stim_check
%                       fieldname1          e.g. 'types' 
%                       fieldvals1          e.g. { 'PULSE', 'SINE' }
% 
% CALLS                 LoadStims, isoverlap, inranges, stim_get, stim_check
%
% examples:
% 10-20 ms pulses of ~40 mA current:
%   out = stim_select( stims( 1 ), 'types', { 'PULSE', 'SINE' }, 'durs', [ 0.01 0.02 ] );
% descending chirps, 0:40 Hz:
%   out = stim_select( stim, 'types', 'ZAP', 'franges', [ 0 40 ], 'dfranges',  [ -inf 0 ] );
%
% 14-feb-13

% revisions
% 28-feb-13 extended to support char arguments
% 04-mar-13 extended to support chirp direction selection
% 02-jan-14 extended to support channel-specific value selection in case of
%               multi-channel simultaneous selection
% 09-mar-15 if empty, mirrored out
% 17-aug-19 cleaned up
% 31-aug-19 replaced call to inrange with call to inranges

function out = stim_select( stim, varargin )

out                     = stim;

% parse input
if ~stim_check( stim ) || length( stim ) ~= 1
    return
end
params                  = cell( 1 );
if length( varargin ) == 1
    if length( varargin{ 1 } ) >= 2
        params          = varargin{ : };
    end
elseif length( varargin ) >= 2
    params              = varargin;
end
if isempty( params )
    return
end
if isempty( stim.types )
    out                 = stim;
    return
end

% determine whether sim. and to load individual structs
nchans                  = length( stim.chan );

% determine relevant entries
sidx                        = true( length( stim.types ), 1 );
fields                      = fieldnames( stim_make );
fields                      = [ fields; 'dfranges' ];
stim.dfranges               = diff( stim.franges, [], 2 );
for i                       = 1 : length( params ) / 2
    field                   = params{ 2 * i - 1 };
    value                   = params{ 2 * i };
    simvals                 = 0;
    if strcmp( field, 'vals' ) && nchans > 1 && size( value, 1 ) == nchans
        [ ~, ~, mstims ]    = LoadStims( stim.filebase );
        chans               = stim.chan;
        valshat             = NaN * ones( size( sidx, 1 ), nchans );
        for j               = 1 : nchans
            chan0           = chans( j );
            for k           = 1 : length( mstims )
                chan1       = mstims( k ).chan;
                if isequal( chan0, chan1 )
                    idx     = stim.index( :, j );
                    valshat( :, j ) = mstims( k ).vals( idx );
                end
            end
        end
        simvals             = 1;
    elseif ~isequal( size( value ), [ 2 2 ] )
        value               = value( : )';
    end
    if ~ismember( field, fields ) || isempty( value )
        continue
    end
    vals                    = getfield( stim, field );
    if isa( vals, 'cell' ) && ischar( vals{ 1 } ) % string field: types
        sidx                = sidx & ismember( vals, value );
    elseif isnumeric( vals )
        if simvals % special case: simultaneous events, filtered by value
            for j           = 1 : nchans
                sidx        = sidx & inranges( valshat( :, i ), value( i, : ) );
            end
        elseif size( vals, 2 ) == 1 % integer field: durs, vals
            if length( value ) == 1 % integer specified
                sidx        = sidx & ismember( vals, value );
            elseif length( value ) == 2 % range specified, e.g. durs, vals, dfranges
                sidx        = sidx & inrange( vals, value );
            end
        elseif strcmp( field, 'times' ) %strcmp( field, 'franges' ) % times
            % any event within the range specified
            sidx = sidx & isoverlap( vals, value, 0 );
        elseif strcmp( field, 'franges' ) % frequency
            if isequal( size( value ), [ 2 2 ] )
                sidx = sidx & inrange( vals( :, 1 ), value( 1, : ) ) & inrange( vals( :, 2 ), value( 2, : ) );
            else
                sidx = sidx & inrange( vals( :, 1 ), value ) & inrange( vals( :, 2 ), value );
            end
        end
    end
end

% assign output
stim                        = rmfield( stim, 'dfranges' );
out                         = stim_get( stim, sidx );

return

% EOF
