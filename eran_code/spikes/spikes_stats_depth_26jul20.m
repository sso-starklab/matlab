% spikes_stats_depth            adds a field 'depth' to the sst structure
%
% call                          [ sst, depth ] = spikes_stats_depth( filebase )
%
% gets                          filebase 
%
% optional arguments (given as name/value pairs):
%                               Overwrite       {-2} determinez compute/Overwrite policy:  
%                                               1   always computes and overwrites
%                                               0   always computes; overwrites only if no depth field
%                                               -1  recomputes only if no depth field, never overwrites
%                                               -2  recomputes only if no depth field, overwrites only if no depth field 
%                                               -3  always computes; never overwrites
%                               graphics        {0} visualization
%
%                               bySpkGrp        {0} 0 uses Anatomical groups (see get_probe)
%                               flipLFP         {0} in case of discrepancy between LFP and SPK order, will flip the LFP
%                               moveNaNs        {-1} move NaNs to end;
%                                               1   moves NaNs to start
%                                               0   keeps where they are
%
%                               sdTH            {5} ripple SD 
%                               freqRange       {[110 220] } ripple bandwidth
%
% returns                       sst             updated sst structure
%                               depth           the new field
%
% files required:               basename.xml
%                               basename.sst
%                               basename.sps
%
% does                          determines the optimal HFO channel from the *sps, 
%                               and the soma location by the sst.geo_com,
%                               then computes the difference
%
% conventions                   (1) positive is above (in probe geometry)
%                               (2) depth is measured in sites (i.e. 20 um for dense, 50 or 100 um for linear)
%                               (3) probe should be sorted by ascending order of channels (i.e. [0 1 2 3 ...]), 
%                                   otherwise should be flipped using flipLFP
%
% calls                         LoadXml                     (blab) 
%                               get_probe                   (formats)
%                               ParseArgPairs, replacetok   (general)
%                               alines, textf               (graph)
%                               inrange                     (stats)

% 31-oct-13 ES

% revisions
% 15-apr-20 cleaned up
% 16-jun-20 added moveNaNs
% 26-jul-20 (1) added argument bySpkGrp, used during call to get_probe
%           (2) updated help to explain Overwrite -3
%           (3) modified dashed line to median
%           (4) modified colors to new conventions and added to constants

function [ sst, depth ] = spikes_stats_depth( filebase, varargin )

% constants
colors                      = [ 46 125 50; 106 27 154 ] / 255; % INT, PYR
colors_light                = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR

% argument and file handling
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    error( 'missing argument filebase' ); 
end
[ Overwrite, graphics...
    , bySpkGrp, flipLFP, moveNaNs ...
    , sdTH, freqRange ...
    ] = ParseArgPairs(...
    { 'Overwrite', 'graphics'...
    , 'bySpkGrp', 'flipLFP', 'moveNaNs'...
    , 'sdTH', 'freqRange'
    }...
    , { -2, 0 ...
    , 0, 0, -1 ...
    , 5, [ 110 220 ] ...
    }...
    , varargin{ : } );
if ~exist( [ filebase '.xml' ], 'file' )
    msg                     = sprintf( '%s: Missing %s file!!', upper( mfilename ), filebase );
    throw( MException( 'spikes_stats_depth:MissingXMLfile', msg ) );
end
if ~exist( [ filebase '.sst' ], 'file' )
    msg                     = sprintf( '%s: Missing %s file!!', upper( mfilename ), filebase );
    throw( MException( 'spikes_stats_depth:MissingSSTfile', msg ) );
end
if ~exist( [ filebase '.sps' ], 'file' )
    msg                     = sprintf( '%s: Missing %s file!!', upper( mfilename ), filebase );
    throw( MException( 'spikes_stats_depth:MissingSPSfile', msg ) );
end
[ ~, filename, extname ]    = fileparts( filebase );
filename                    = [ filename extname ];

% load everything
par                         = LoadXml( filebase );
probe                       = get_probe( par, bySpkGrp );
sps                         = load( [ filebase '.sps' ], '-mat' );
load( [ filebase '.sst' ], 'sst', '-mat' )

% check if needed to re-compute
if Overwrite < 0 && isfield( sst, 'depth' ) && Overwrite ~= -3
    depth                   = sst.depth;
    return
end

% unit shanks
shank                       = sst.shankclu( :, 1 ); 
nclu                        = length( shank );

% threshold HFOs
sd                          = sps.stats( :, 7 );
freq                        = sps.stats( :, 5 );
vidx                        = sd > sdTH & inrange( freq, freqRange );
layer                       = NaN * ones( size( vidx ) );
shanks                      = sps.stats( :, 1 );
shanks( shanks == 0 )       = NaN;
layer( vidx )               = sps.stats( vidx, 2 );
vidx( isnan( layer ) )      = 0;

% error checking
shank( ~ismember( shank, shanks ) ) = NaN;
if length( layer ) ~= size( probe, 2 )
    error( 'mismatch between probe and sps - check!' )
end
if sum( ~ismember( shank( ~isnan( shank ) ), shanks ) )
    error( 'additional shanks in unit list - modify code to ignore those!' )
end

% assign depths
if flipLFP || ismember( filename, { 'm100_1', 'm100_2', 'm101_1', 'm101_2', 'm101_3' } )
    probe                   = flipud( probe );
end
if moveNaNs
    for i                       = 1 : size( probe, 2 )
        v = probe( :, i ); nans = isnan( v ); 
        if any( nans )
            if moveNaNs == -1 % move NaNs to start
                probe( :, i ) = [ v( nans ); v( ~nans ) ]; 
            elseif moveNaNs == 1 % move NaNs to end
                probe( :, i ) = [ v( ~nans ); v( nans ) ]; 
            end
        end
    end
end
[ ii, ~ ]                   = find( ismember( probe, layer ) );
layer( shanks( vidx ) )     = ii;
elayer                      = NaN * ones( nclu, 1 );
for i = 1 : length( layer )
    if ~vidx( i )
        continue
    end
    elayer( shank == shanks( i ) ) = layer( i );
end
depth                       = sst.geo_com - elayer;
sst0                        = sst;
sst.depth                   = depth;

% overwrite
if Overwrite == 1 || ( Overwrite ~= -1 && Overwrite ~= -3 && ~isfield( sst0, 'depth' ) )
    save( [ filebase '.sst' ], '-v6' )
end

if graphics
    fig                     = figure;
    for ct                  = 0 : 1
        subplot( 2, 1, ct + 1 ),
        col                 = colors( ct + 1, : );
        col_light           = colors_light( ct + 1, : );
        [ h, b ]            = hist( sst.depth( sst.pyr == ct ), -7 : 7 );
        bh                  = bar( b, h, 'EdgeColor', col, 'FaceColor', col );
        xlim( [ -7 7 ] ),
        ylims = ylim;
        ylim( [ 0 ylims( 2 ) + 1 ] )
        alines( 0, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
        alines( nanmedian( sst.depth( sst.pyr == ct ) ), 'x', 'color', col_light, 'linestyle', '--' );
        set( gca, 'tickdir', 'out', 'box', 'off' )
        axis square
    end
    textf( 0.5, 0.975, replacetok( filename, '\_', '_' ) );
end

return

% EOF

% to call with plotting, without saving, with flipping (as we usually do in
% postprocess_spikes):
sst                 = spikes_stats_depth( filebase, 'graphics', 1,  'Overwrite', -3, 'flipLFP', 1 );

