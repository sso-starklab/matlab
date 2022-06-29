% load_spikes           according to geometry and quality
% 
% CALL                  [ s, sst ] = load_spikes( filebase, shanknums )
%
% RETURNS               s           structure w/ fields: filebase, shankclu, clu, res, and map
%
% CALLS                 determine_units, get_spikes

% 03-mar-13 ES

% revisions
% 20-nov-13 extended to support loading junk units ('E' includes noise cluster, 'F' includes artifact too. 
% 18-aug-19 cleaned up
% 05-jul-20 organized second output to match first

function [ s, sst ] = load_spikes( filebase, shanknums, varargin )

% constants
minclu                      = 2;

% arguments
nargs                       = nargin;
if nargs < 2 || isempty( shanknums )
    shanknums               = [];
end
if ~isempty( varargin )
    if isequal( varargin{ : }, 'E' )
        minclu              = 1;
        varargin{ : }       = 'D';
    elseif isequal( varargin{ : }, 'F' )
        fprintf( '%s: finish writing this...\n', upper( mfilename ) )
        varargin{ : }       = 'D';
        minclu              = 0;
    end
end

% determine units
[ shankclu, sst ]           = determine_units( filebase, shanknums, varargin{ : } );

% get the spikes
if isempty( shankclu )
   clu                      = [];
   res                      = [];
   map                      = [];
else
    [ clu, res, map ]       = get_spikes( filebase, shankclu, 'clu', minclu );
    [ ~, ~, i2 ]            = intersect( map( :, 2 : 3 ), shankclu( :, 1 : 2 ), 'rows' ); 
    shankclu                = shankclu( i2, : );
    idx                     = ismember( sst.shankclu, shankclu( :, 1 : 2 ), 'rows' ); 
    sst                     = struct_select( sst, idx );
end

% generate the output structure
s.filebase                  = filebase;
if isempty( varargin )
    s.ilevel                = 'null';
else
    s.ilevel                = varargin{ 1 };
end
s.shankclu                  = shankclu;
s.clu                       = clu;
s.res                       = res;
s.map                       = map;

return

% EOF
