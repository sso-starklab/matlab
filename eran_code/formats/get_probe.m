% get_probe             get probe layout in matrix format
%
% CALL                  probe = get_probe( par, spkGrp )
%
% GETS                  par             par structure
%                       spkGrp          {1}: generate the layout based on the spike groups
%                                           0: generate based on the anatomical groups (xml file)
%                                               or the neuronal channels (prm.xml file)
% OPTIONAL (name/value pairs)
% 
%                       linearize       {0}
%                       keepGrps        {[]}
%                       tosort          {0}
%                       toflip          {0}
%                       tospace         {1}
%
% tospace               place missing channels in their anatomical position.
%                           can be done without external knowledge only if
%                           channels are sorted (i.e. probe is 1,2,..10 on each
%                           shank)
%
% RETURNS               probe             matrix, groups in columns
%
% CALLS                 ParseArgPairs
%
% 03-apr-13 ES

% revisions
% 13-oct-13 post-hoc manipulations added
% 17-aug-19 cleaned up

function probe = get_probe( par, spkGrp, varargin )

probe                       = [];
nargs                       = nargin;
if ~isa( par, 'struct' ) || ~isfield( par, 'AnatGrps' ) || ~isfield( par, 'ElecGp' )
    return;
end
if nargs < 2 || isempty( spkGrp )
    spkGrp                  = 1;
end

% post-hoc manipulations (optional)
[ linearize, keepGrps, tosort, toflip, tospace...
    ] = ParseArgPairs(...
    { 'linearize', 'keepGrps', 'tosort', 'toflip', 'tospace'...
    }...
    , { 0, [], 0, 0, 1 ...
    }...
    , varargin{ : } );

% get the channels
if spkGrp
    
    nSpkGps                 = length( par.SpkGrps );
    nchansGps               = zeros( nSpkGps, 1 );
    for i                   = 1 : nSpkGps
        nchansGps( i )      = length( par.SpkGrps( i ).Channels );
    end
    probe                   = NaN * ones( max( nchansGps ), nSpkGps );
    for i                   = 1 : nSpkGps
        probe( 1 : nchansGps( i ), i ) = par.SpkGrps( i ).Channels + 1;
    end
    
else
    % assume a multi-shank probe, each shank numbered from bottom to top
    if isfield( par.AnatGrps( 1 ), 'Type' ) % prm.xml file
        
        nAnatGps            = length( par.AnatGrps );
        nchansGps           = zeros( nAnatGps, 1 );
        for i               = 1 : nAnatGps
            nidx            = ismember( par.AnatGrps( i ).Type, 'neuronal' );
            nchansGps( i )  = sum( nidx );
        end
        kidx                = nchansGps > 0;
        AnatGps             = find( kidx );
        nchansGps           = nchansGps( kidx );
        nAnatGps            = length( AnatGps );
        probe               = NaN * ones( max( nchansGps ), nAnatGps );
        for i               = 1 : nAnatGps
            nidx            = ismember( par.AnatGrps( AnatGps( i ) ).Type, 'neuronal' );
            probe( 1 : nchansGps( i ), i ) = par.AnatGrps( AnatGps( i ) ).Channels( nidx ) + 1;
        end
        
    else % xml file
        
        nElecGps            = par.nElecGps;
        nchansGps           = zeros( nElecGps, 1 );
        for i               = 1 : nElecGps
            nchansGps( i )  = length( par.ElecGp{ i } );
        end
        probe               = NaN * ones( max( nchansGps ), nElecGps );
        for i               = 1 : nElecGps
            probe( 1 : nchansGps( i ), i ) = par.ElecGp{ i } + 1;
        end
        
    end
    
end


% make post-hoc manipulations

% remove unnecessary groups
if isempty( keepGrps )
    keepGrps                = 1 : size( probe, 2 );
end
keepGrps                    = intersect( keepGrps, 1 : size( probe, 2 ) );
probe                       = probe( :, keepGrps );

% sort then flip
if tosort
    probe                   = sort( probe );
end
if toflip
    probe                   = flipud( probe );
end

% space if sorted
if tospace
    probe( all( isnan( probe ), 2 ), : ) = [];
    nans                    = isnan( probe );
    nrows                   = size( probe, 1 );
    for col                 = find( sum( nans ) > 0 )
        vals                = probe( ~nans( :, col ), col );
        if issorted( vals )
            fvals           = ( ( nrows * ( col - 1 ) + 1 ) : ( nrows * col ) )';
            [ ~, idx ]      = ismember( vals, fvals );
            if length( idx ) ~= nrows
                continue
            end
            newcol          = NaN * ones( nrows, 1 );
            newcol( idx )   = vals;
            probe( :, col ) = newcol;
        end
    end
end

% linearize 
if linearize == 1
    probe                   = probe( : );
end

return

% EOF