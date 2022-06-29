% spikes_stats_supplement   add geometric dispersion and signed extremum to sst structure
%
% call                      sst = spikes_stats_supplement_bpi( sst )
%
% gets and returns          sst structure
%
% optional arguments (given as name/value pairs):
%                           nsites0         {32} 
%                           PEAKSAMPLE      {16}
%
% calls                     calc_com, fft_upsample, ParseArgPairs

% 30-mar-21 SSo


function sst = spikes_stats_supplement_bpi( sst, varargin )

% default values
nsites0              = 32;                   % maximal number of sites per shank
PEAKSAMPLE0          = 16;                   % fixed given par.SpkGrps.PeakSample
% argument handling
nargs                       = nargin;
if nargs < 1 || isempty( sst )
    return
end
[ nsites, PEAKSAMPLE ...
    ]                       = ParseArgPairs(...
    { 'USF', 'PEAKSAMPLE' ...
    }...
    , { nsites0, PEAKSAMPLE0 ...
    }...
    , varargin{ : } );

% allocate space
n                   = size( sst.shankclu, 1 );
nsites              = length(sst.mean{1,1});
z0                  = NaN( n, nsites );
vA                  = z0;
vT                  = z0;
bpi                 = NaN( n, 1 );
tV                  = z0;
vB                  = z0;
pA                  = z0;
nA                  = z0;
nT                  = z0;

% go over units
for i               = 1 : n
    w              = sst.mean{ i };
    [ vA, vT, bpi, tV, vB, pA, nA, pT, nT ] = calc_spatial_waveform_features( w, [], [], 0 );

    % update the structure
    sst.vA (i, 1 : length( vA ) )                        = vA';
    sst.vT (i, 1 : length( vT ))                         = vT';
    sst.bpi (i)                                          = bpi';
    sst.tV (i, 1 : length( tV ))                         = tV';
    sst.vB (i, 1 : length( vB ))                         = vB';
    sst.pA (i, 1 : length( pA ))                         = pA';
    sst.nA (i, 1 : length( nA ))                         = nA';
    sst.pT (i, 1 : length( pT ))                         = pT';
    sst.nT (i, 1 : length( nT ))                         = nT';
end


return

% EOF


% use example #1 (look at a specific dataset, do not save anything):

% load sst for B and above units
[ shankclu, sst ]   = determine_units( filebase, [], 'B' );
[ imat, i1, i2 ]    = intersect( sst.shankclu, shankclu( :, 1 : 2 ), 'rows' );
idx1                = false( size( sst.shankclu, 1 ), 1 );
idx1( i1 )          = 1;
sst                 = struct_select( sst, idx1 );

% add the new fields
sst2                = spikes_stats_supplement( sst );

% calculate some statistics...

% use example #2 (update the sst on the disk):
sstfname            = [ filebase '.sst' ];
load( sstfname, '-mat' )
sst                 = spikes_stats_supplement_bpi( sst );
save( sstfname, 'sst' )

