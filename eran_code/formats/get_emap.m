% get_emap              map channels to electrode groups/sites (matrix)
%
% call                  [ map, pmap ] = get_emap( filebase, tosort )
%
% gets                  filebase        either full path + base or a par structure
%                       tosort          {1}
%
% returns               map             [ egroup site_in_egroup channel ]
%
% does                  for conventional arrays, this routine uses the par
%                       when there are *ind* files, they are used instead
%
% calls                 LoadXml                                 (blab)
%                       LoadInd                                 (formats)

% 04-apr-13 ES

% revisions
% 17-sep-19 (1) extended to support *ind* files
%           (2) channel numbers sorted by default
%           (3) argument checks added

function [ map, pmap ] = get_emap( filebase, tosort )

map                     = [];
pmap                    = [];

nargs = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( tosort )
    tosort              = 1;
end

if isa( filebase, 'struct' ) && isfield( filebase, 'SpkGrps' ) && isfield( filebase, 'nChannels' )
    par                 = filebase;
    filebase            = par.FileName;
else
    par                 = LoadXml( filebase );
end

% make a channel mapping matrix
ntot                    = par.nChannels;
map                     = NaN * ones( ntot, 3 );
si                      = 1;
for i                   = 1 : length( par.SpkGrps )
    indfname            = sprintf( '%s.ind.%d', filebase, i );
    if exist( indfname, 'file' )
        ind             = LoadInd( indfname );
        chans           = setdiff( unique( ind ), 0 );
    else
        chans           = par.SpkGrps( i ).Channels + 1;
    end
    if tosort
        chans           = sort( chans );
    end
    nchans              = length( chans );
    ei                  = si + nchans - 1;
    map( si : ei, : )   = [ i * ones( nchans, 1 ) ( 1 : nchans )' chans( : ) ];
    si                  = ei + 1;
end
if ( ei + 1 ) < ntot
    map( ei + 1 : ntot, : ) = [];
end

% make an anatomical matrix (assume dense):
if nargout > 1
    ugroups             = unique( map( :, 1 ) );
    usites              = unique( map( :, 2 ) );
    pmap                = NaN * ones( length( usites ), length( ugroups ) );
    for i               = 1 : length( ugroups )
        idx                         = map( :, 1 ) == ugroups( i );
        pmap( map( idx, 2 ), i )    = map( idx, 3 );
    end
end

return

% EOF

% usage example: get the anatomical location of waveforms

% get the list of intra-group sites for each unit
load( [ filebase '.sst' ], '-mat' );
loc                         = [ sst.shankclu( sidx, 1 ) round( sst.geo_com( sidx ) ) ];
% get the map - for each spike group, list the intra-group channels, and the overall channel
par                         = LoadXml( filebase );
emap                        = get_emap( par );
% get the overall channels for each unit
sidx                        = ismember( sst.shankclu, shankclu( :, 1 : 2 ), 'rows' );
[ ~, idx ]                  = ismember( loc, emap( :, 1 : 2 ), 'rows' );
loc_chans                   = emap( idx, 3 )

