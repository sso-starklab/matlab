% get_spikes        times and lables, of the specified file and units
%
% CALL              [ clu, res, map ] = get_spikes( filebase, shankclu, suffix )
%
% GETS              filebase          full base or par structure
%                   shankclu          two-column matrix, see determine_units
%                                         {[]}, then take spikes of all units from all shanks
%                   suffix            {'clu'}; spike files of the type 
%                                         filebase.res.ShankNum
%                                         filebase.suffix.ShankNum
%                                     will be used
%
% RETURNS           clu               labels (n-element vector)
%                   res               spike times (same size as clu)
%                   map               3-column matrix mapping the clu to the original shankclu
%
% CALLS             LoadClu, LoadRes, LoadXml, verb
%
% NOTE              repetitive calls with different shankclu will result 
%                   in the same unit given different labels in Clu; 
%                   check the Map for the actual source
% 
% See also          determine_units, get_triggers


% 25-feb-13 ES

% revisions
% 07-mar-13 for single-shank loading, map is the identity map
% 14-mar-13 return the modified shankclu
% 22-jun-15 case of empty Map (no units) handled
% 18-aug-19 cleaned up

function [ Clu, Res, Map, shankclu ] = get_spikes( filebase, shankclu, suffix, minclu )

% constants
verbose                     = 1;
clustr0                     = '.clu.';          % backup source..

% initialize
Clu                         = [];
Res                         = [];
Map                         = [];

%----------------------------------------------------------------------%
% parse input
nargs                       = nargin;
if isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
    par                     = LoadXml( filebase );
elseif isa( filebase, 'struct' )
    if isfield( filebase, 'FileName' )
        par                 = filebase;
        filebase            = filebase.FileName;
    else
        fprintf( '%s: missing filebase\n', upper( mfilename ) );
        return
    end
else
    fprintf( '%s: missing filebase\n', upper( mfilename ) );
    return
end

if nargs < 2 || isempty( shankclu )
    shankclu                = [];
end
if isvector( shankclu ) && ( numel( shankclu ) > 3 || size( shankclu, 2 ) == 1 )
    shanknums               = shankclu( : );
    shankclu                = [];
elseif ~isempty( shankclu ) && size( shankclu, 2 ) <= 3
    shanknums               = unique( shankclu( :, 1 ) );
else
    shankclu                = [];
    shanknums               = [];
end
if isempty( shanknums )
    shanknums               = 1 : par.nElecGps;
end

if nargs < 3 || isempty( suffix )
    suffix                  = 'clu';
end
clustr                      = sprintf( '.%s.', suffix );

if nargs < 4 || isempty( minclu )
    if ~isempty( shankclu )
        minclu              = min( shankclu( :, 2 ) );
    else
        minclu              = 2;
    end
end

%----------------------------------------------------------------------%
% load all spikes of the relevant shanks
[ ~, filename, extname ]    = fileparts( filebase );
verb( sprintf( '%s: Loading %s%s spike data... ', upper( mfilename ), filename, extname ), -verbose )

maxClu                      = 0;
for i                       = shanknums( : ).'
    % load single-shank data
    clufile                 = [ filebase clustr num2str( i ) ];
    if ~exist( clufile, 'file' )
        clufile             = [ filebase clustr0 num2str( i ) ];
        if ~exist( clufile, 'file' )
            continue
        end
    end
    resfile                 = [ filebase '.res.' num2str( i ) ]; 
    clu                     = LoadClu( clufile );
    res                     = LoadRes( resfile );   
    % remove artifacts/noise
    indx                    = clu >= minclu;
    if ~sum( indx )
        continue
    end
    clu                     = clu( indx );
    res                     = res( indx );
    % renumber
    if length( shanknums ) == 1
        uniclu              = unique( clu );
        umap                = uniclu;
    else
        [ uniclu, ~, clu ]  = unique( clu );
        clu                 = maxClu + clu;
        umap                = unique( clu );
    end
    % accumulate and update map
    Clu                     = [ Clu; clu ];
    Res                     = [ Res; res ];
    Map                     = [ Map; [ umap i * ones( size( umap ) ) uniclu ] ];
    maxClu                  = max( Clu );
end
[ Res, sidx ]               = sort( Res );
Clu                         = Clu( sidx );

%----------------------------------------------------------------------%
% keep only the desired units
if isempty( shankclu )
    shankclu                = Map( :, 2 : 3 );
end
if isempty( Map )
    return
end
ridx0                       = ~ismember( Map( :, 2 : 3 ), shankclu( :, 1 : 2 ), 'rows' );
ridx                        = ismember( Clu, Map( ridx0, 1 ) );
Map( ridx0, : )             = [];
Res( ridx )                 = [];
Clu( ridx )                 = [];
shankclu( ~ismember( shankclu( :, 1 : 2 ), Map( :, 2 : 3 ), 'rows' ), : ) = [];

verb( sprintf( 'A total of %d/%d spikes from %d/%d units (%d shanks) were loaded.'...
    , length( Res ), length( ridx ), size( shankclu, 1 )...
    , length( ridx0 ), length( shanknums ) ), verbose )

return

% EOF
