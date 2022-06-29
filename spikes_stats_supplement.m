% spikes_stats_supplement   add geometric dispersion and signed extremum to sst structure
%
% call                      sst = spikes_stats_supplement( sst )
%
% gets and returns          sst structure
%
% optional arguments (given as name/value pairs):
%                           USF             {32} 
%                           PEAKSAMPLE      {16}
%
% calls                     calc_com, fft_upsample, ParseArgPairs

% 11-sep-19 ES

% revisions
% 09-jul-20 (1) changed spatial metric from peak-to-peak to sum-of-squares
%               (influences all three metrics)
%           (2) changed algorithm for FWHM computation

function sst = spikes_stats_supplement( sst, varargin )

% default values
USF0                 = 32;                   % usampling factor - increase to increase spatial resolution
PEAKSAMPLE0          = 16;                   % fixed given par.SpkGrps.PeakSample

% argument handling
nargs                       = nargin;
if nargs < 1 || isempty( sst )
    return
end
[ USF, PEAKSAMPLE ...
    ]                       = ParseArgPairs(...
    { 'USF', 'PEAKSAMPLE' ...
    }...
    , { USF0, PEAKSAMPLE0 ...
    }...
    , varargin{ : } );

% allocate space
n                   = size( sst.shankclu, 1 );
z0                  = zeros( n, 1 );
geo_sd              = z0;
geo_fwhm            = z0;
extremum            = z0;

% go over units
for i               = 1 : n
    
    % p2p by recording site
    sw              = sum(  sst.mean{ i }.^2, 2 );                                                % spatial signature of each waveform

    % determine if positive/negative
    [ ~, midx ]     = max( sw );
    ext             = sst.mean{ i }( midx, PEAKSAMPLE );                            % [microVolts]
    extremum( i )   = ext;
    
    % geometrical dispersion by SD of COM
    nchans          = length( sw );                                                 % number of channels for this unit
    [ ~, geo_sd( i ) ] = calc_com( ( 1 : nchans )', sw );                           % [sites]
    
    % geometrical dispersion by FWHM
    zu              = fft_upsample( sw( : ), USF, 1 );
    hm              = ( max( zu ) - min( zu ) ) / 2 + min( zu );
    fwhm            = sum( zu >= hm ) / USF;                                        % [sites]
    geo_fwhm( i )   = fwhm;
    
end

% update the structure
sst.geo_sd          = geo_sd;
sst.geo_fwhm        = geo_fwhm;
sst.extremum        = extremum;

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
sst                 = spikes_stats_supplement( sst );
save( sstfname, 'sst' )

