% get_egroup            from par/filebase and channel number
%
% call                  [ egroup, chans, achans, par ] = get_egroup( par, chan, bySpkGrp )
%
% gets                  par         *.xml
%                       chan        {[]}; if given, gives chans for the
%                                       single spike group that contains chan
%                       bySpkGrp    {1}; if 0, bases the selection on the
%                                       single electrode group that contains chan
% 
% returns               egroup        electrode group
%                       chans         spike group channels
%                       achans        anatomical group channels
%
% calls                 LoadXml

% 05-nov-12 ES

% revisions
% 14-jul-20 (1) added bySpkGrp
%           (2) changed LoadPar to LoadXml

function [ egroup, chans, achans, par ] = get_egroup( par, chan, bySpkGrp )

% initialize output
egroup              = [];
chans               = [];
achans              = [];

% argument handling
nargs                   = nargin;
if nargs < 1 || isempty( par )
    return
end
if nargs < 2 || isempty( chan )
    chan                = [];
else
    chan                = chan( 1 );
end
if nargs < 3 || isempty( bySpkGrp )
    bySpkGrp            = 1;
end

% load and check the par, populate all relevant fields
if ~isa( par, 'struct' ) && isa( par, 'char' ) && exist( fileparts( par ), 'dir' )
    par                 = LoadXml( par );
end
if ~isa( par, 'struct' )
    return
end 
if ~isfield( par, 'SpkGrps' ) && isfield( par, 'AnatGrps' )
    par.SpkGrps         = par.AnatGrps;
end
if ~isfield( par, 'ElecGp' )
    for i               = 1 : length( par.SpkGrps )
        par.ElecGp{ i } = par.SpkGrps( i ).Channels;
    end
end

% go over egroups and accumulate channels
nSpkGrps                = length( par.SpkGrps );
nAnatGrps               = length( par.AnatGrps );
for i                   = 1 : length( par.ElecGp )
    if isempty( chan )
        if i <= nSpkGrps
            chans       = [ chans par.SpkGrps( i ).Channels + 1 ];
        end
        if i <= nAnatGrps
            achans      = [ achans par.AnatGrps( i ).Channels + 1 ];
        end
    else
        if bySpkGrp && sum( chan == ( par.SpkGrps( i ).Channels + 1 ) )
            egroup          = i; 
        elseif ~bySpkGrp && sum( chan == ( par.AnatGrps( i ).Channels + 1 ) )
            egroup          = i;
        end
        if ~isempty( egroup )
            if i <= nSpkGrps
                chans       = par.SpkGrps( egroup ).Channels + 1;
            else
                chans       = NaN;
            end
            if i <= nAnatGrps
                achans      = par.AnatGrps( egroup ).Channels + 1;
            else
                achans      = NaN;
            end
            break
        end     % after detection of chan
    end         % detection attempt
end             % egroups

return

% EOF
