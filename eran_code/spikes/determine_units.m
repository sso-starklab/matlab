% determine_units       to use by filebase and class
%
% shankclu = determine_units( filebase, shankclus, ilevel, params )
% 
% DOES:
% loads (or computes) the sst for the merged directory and 
% selects units based on the class/parameters
%
% ARGUMENTS:
% filebase          full path
% shankclus         anatomical selection: 
%                       {[]}: all valid units will be loaded
%                       vector (if two element, column vector) - only valid
%                           units from these shanks/elec groups will be used
%                       matrix (or 2-element row vector) - only the
%                           specified units will be used (validity not checked,
%                           just existence and cell type) 
% ilevel            {'B'}; class. Additional arguments can be specified 
%                           as in check_cluster_quality
% 
% OUTPUT:
% shankclu          [ elec_group cluster_number cell_type ]
% 
% calls             get_merged_filenum, check_cluster_quality, spikes_stats
%
% see also          get_spikes

% 25-feb-13 ES

% revisions
% 20-mar-13 secondary output sst (loaded anyhow)
% 04-apr-13 clean up sst format bug (from April. 2012)
% 16-apr-13 enable processing by local sst; output the classification
% 01-aug-19 modified call to spikes_stats to match field/value pairs

function [ shankclu, sst, cls ] = determine_units( filebase, shankclus, varargin )

% initialize
shankclu            = [];
sst                 = [];
nargs               = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( shankclus )
    shankclus       = [];
end

% parse the input argument shankclus
if isvector( shankclus ) && ( numel( shankclus ) > 2 || size( shankclus, 2 ) == 1 )
    shanks          = shankclus( : );
    shankclus       = [];
elseif size( shankclus, 2 ) == 2 || size( shankclus, 2 ) == 3
    shanks          = unique( shankclus( :, 1 ) );
else
    shanks          = [];
    shankclus       = [];
end

% get the sst structure (from the merged file unless explicitly given an sst)
if isa( filebase, 'char' ) && exist( fileparts( filebase ), 'dir' )
    mfilebase       = get_merged_filenum( filebase );
    %sst = spikes_stats( mfilebase, [], [], [], [], 0 );
    sst = spikes_stats( mfilebase, 'Overwrite', -2 );
elseif isa( filebase, 'struct' ) && isfield( filebase, 'shankclu' ) && isfield( filebase, 'pyr' )
    sst             = filebase;
else
    return
end
try
    shankclu        = [ sst.shankclu sst.pyr ];
catch
    fprintf( '%s: cleaning up file %s\n', upper( mfilename ), [ mfilebase '.sst' ] )
    sst             = sst_cleanup( filebase, 1 );
    shankclu        = [ sst.shankclu sst.pyr ];
end

% get the list of good units and their classification:
if ~isempty( shankclus )
    shankclu        = shankclu( ismember( shankclu( :, 1 : 2 ), shankclus( :, 1 : 2 ), 'rows' ), : );
else
    gidx            = check_cluster_quality( sst, varargin{ : } );
    shankclu        = shankclu( gidx, : );
end
if ~isempty( shanks ) && ~isempty( shankclu )
    shankclu        = shankclu( ismember( shankclu( :, 1 ), shanks ), : );
end

% get the general classification
if nargout > 2
    classes         = { 'A', 'B', 'C', 'D' }; 
    n               = length( classes );
    m               = size( sst.shankclu, 1 );
    lg              = zeros( m, n );
    for i           = 1 : n
        lg( :, i )  = check_cluster_quality( sst, classes{ i } ); 
    end
    cls             = zeros( m, 1 ); 
    for i           = n : -1 : 1
        cls( lg( :, i ) == 1 ) = i; 
    end
end

return

% EOF
