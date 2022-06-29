% rips_select           a subset of events
%
% CALL                  rips = rips_select( rips, sidx )
% 
% CALLS                 struct_select
% 
% See also              find_hfos

% 10-jul-13 ES

% revisions
% 18-aug-19 cleaned up

function rips = rips_select( rips, sidx )

% arguments
if nargin < 2 || isempty( rips ) || isempty( sidx ) || ~isa( rips, 'struct' ) || ~isa( sidx, 'logical' )
    return
end
sidx                        = sidx( : );
ridx                        = ~sidx;
if sum( ridx ) == 0
    return
end
if size( sidx, 1 ) ~= length( rips.t )
    fprintf( '%s: check input size!!\n', upper( mfilename ) )
end
if size( sidx, 1 ) <= 2
    if sum( sidx ) ~= size( sidx, 1 )
        rips.seg            = [];
        rips.t              = [];
        rips.pow            = [];
        rips.sd             = [];
        rips.f              = [];
        rips.edges          = [];
        rips.trigs          = [];
        rips.peaks          = [];
    end
    return
end

% remove all irrelevant events
rips                        = struct_select( rips, sidx );
if isempty( rips.t )
    return
end

% handle the peaks
rips.peaks( ismember( rips.peaks( :, 1 ), find( ridx ) ), : ) = [];
allpeaks                    = rips.peaks;
tmat                        = [ unique( allpeaks( :, 1 ) ) ( 1 : length( rips.t ) )' ];
for i                       = 1 : size( tmat, 1 )
    allpeaks( ismember( allpeaks( :, 1 ), tmat( i, 1 ) ), 1 ) = tmat( i, 2 );
end
rips.peaks                  = allpeaks;

return

% EOF
