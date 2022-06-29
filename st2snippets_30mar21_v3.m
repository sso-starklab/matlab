% st2snippets                   extracts waveforms for two set of channels during one clu
%
% CALL                          [ mat2, vec, t, mat1 ] = st2snippets( filebase, shankclu1, shankclu2 )
%
% GETS                          filebase
%                               shankclu1
%                               shankclu2
%
% OPTIONAL ARGUMENTS            halfWidthMS         {1.5}
%                               graphics            {1}
%                               nplot               {100}
%
% CALLS                         LoadXml
%                               load_spikes, mydetrend
%                               intersectranges
%
% does: 
% (1) determines times of shankclu1 spikes
% (2) determines a single channel for COM of shankclu2 (call this channel2)
% (3) extracts wide-band data from channel2 during shankclu1 spikes (requires halfWidthMS)
% (4) determines for each event whether there is a shankclu2 spike (requires 'causality')
% 
% returns                       mat1            [microV]; m x n1 matrix of n1 events of shankclu1 in channel1, 
%                                                           each of m samples
%                               mat2            [microV]; m x n2 matrix of n2 events in channel2, in shackclu1 times, 
%                                                           each of m samples
%                               vec             logical; 1 x n indicator of whether the event occurred simultaneouly 
%                                                with a shankclu1 spike
%                               t               [ms]; m x 1 vector of time, from -halfWidthMS : halfWidthMS, at spkFs

% 29-mar-21 ES + SSo

% revisions
% 29-mar-21 (1) added documentation
%           (2) scaling to microV
% 30-mar-21 (1) added file checks 
%           (2) modified disk handling to block structure
%           (3) explicit type casting
%           (4) computed cISI and tagged waveforms properly
%           (5) updated graphics

% 30-mar-21 to do:
% (1) update documentation
% (2) in future - may want to allow external determination of channel1, 2

function [ mat2, vec, t, mat1 ] = st2snippets( filebase, shankclu1, shankclu2, varargin)

% argument handling
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( shankclu1 )
    return
end
if nargs < 3 || isempty( shankclu2 )
    return
end

[ halfWidthMS, blocksize, occupancyMS ...
    , graphics, nplot ] = ParseArgPairs(...
    { 'halfWidthMS', 'blocksize', 'occupancyMS' ...
    , 'graphics', 'nplot' }...
    , { 1.5, 100, [] ...
    , 1, 100 }...
    , varargin{ : } );
if isempty( occupancyMS )
    occupancyMS                 = [ 0 halfWidthMS ];
end
occupancyMS                     = intersectranges( [ -1 1 ] * halfWidthMS, occupancyMS );

% file handling
par                             = LoadXml( filebase );
datfname                        = [ filebase '.dat' ];
sstfname                        = [ filebase '.sst' ];
if ~exist( datfname, 'file' ) || ~exist( sstfname, 'file' )
    error( 'missing source files' )
end

% parameter from *xml (amplitude scaling)
spkFs                           = par.SampleRate;
nchans                          = par.nChannels;
p2pvoltage                      = par.VoltageRange;
nBits                           = par.nBits;
Amplification                   = par.Amplification;
scalefactor                     = 1 / 2^nBits * p2pvoltage * 1e6 / Amplification;           % a2d units -> microVolts

% (1) get spike times
s1                              = load_spikes( filebase, shankclu1 );
s2                              = load_spikes( filebase, shankclu2 );
st1                             = s1.res;
st2                             = s2.res;
if ~issorted( st1 )
    st1                         = sort( st1 );
end
if ~issorted( st2 )
    st2                         = sort( st2 );
end
nspikes1                        = length( st1 );
nspikes2                        = length( st2 );
vec                             = false( 1, nspikes1 );

% (2) determine channels for shankclu1 and shankclu2
load( sstfname, '-mat', 'sst' )
idx1                            = ismember( sst.shankclu, shankclu1, 'rows' );
gcom1                        	= sst.geo_com( idx1 );
schannels1                    	= sort( par.SpkGrps( shankclu1( 1 ) ).Channels );
channel1                        = schannels1( round( gcom1 ) );
idx2                          	= ismember( sst.shankclu, shankclu2, 'rows' );
gcom2                           = sst.geo_com( idx2 );
schannels2                      = sort( par.SpkGrps( shankclu2( 1 ) ).Channels );
channel2                        = schannels2( round( gcom2 ) );

% (3) prepare to extract wide-band data during st1 events
halfWidth                       = round( halfWidthMS * spkFs / 1000 );      % [samples]
win                             = ( -halfWidth : halfWidth )';
t                               = win / spkFs * 1000;
nwin                            = length( win );
owin                            = round( occupancyMS * spkFs / 1000 );      % [samples]

% make blocks
bs                              = ( 1 : blocksize : nspikes1 )'; 
nblocks                         = length( bs ); 
be                              = [ bs( 2 : nblocks ) - 1; nspikes1 ];
blocks                          = [ bs unique( be ) ];

% (4) go over st1 spikes and extract channel1 and channel2 data
a                               = memmapfile( datfname, 'Format', 'int16' );
nsamples                        = length( a.Data ) / nchans;
if nsamples ~= round( nsamples )
    error( 'problem with dat file' )
end
mat1                            = NaN( nwin, nspikes1 );
mat2                            = NaN( nwin, nspikes1 );
fprintf( 1, 'Extracting waveforms for %d spikes from channels %d and %d ', nspikes1, channel1, channel2 )
for bidx                        = 1 : nblocks
    if ~mod( bidx, 10 )
        fprintf( 1, '. ' )
    end
    nspk                        = diff( blocks( bidx, : ) ) + 1;
    spkidx                      = blocks( bidx, 1 ) : blocks( bidx, 2 );
    sidx1                       = ( ones( nwin, 1 ) * st1( spkidx )' - 1 + win * ones( 1, nspk ) ) * nchans + channel1;
    sidx2                       = ( ones( nwin, 1 ) * st1( spkidx )' - 1 + win * ones( 1, nspk ) ) * nchans + channel2;
    ev1                         = double( a.Data( sidx1 ) );
    ev2                         = double( a.Data( sidx2 ) );
    mat1( :, spkidx )           = ev1;
    mat2( :, spkidx )           = ev2;
end
fprintf( 1, 'done!\n' )
clear a

% scale and detrend each waveform
mat1                            = mydetrend( mat1 * scalefactor );
mat2                            = mydetrend( mat2 * scalefactor );

% (5) determines for each event whether there is a shankclu2 spike
% combine spike trains
res                             = [ st1; st2 ];
clu                             = [ ones( nspikes1, 1 ); 2 * ones( nspikes2, 1 ) ];
[ res, sidx ]                   = sort( res );
clu                             = clu( sidx );

% compute cross-intervals
pidx                            = find( diff( clu ) == 1 );                 % indices of st1 right before st2
pcisi                           = res( pidx + 1 ) - res( pidx );            % positive C-ISIs
nidx                            = find( diff( clu ) == -1 );                % indices of st1 right after st2
ncisi                           = res( nidx ) - res( nidx + 1 );
mat                             = [ [ pidx pcisi ]; [ nidx ncisi ] ];       % combined c-ISIs
cisi                            = sortrows( mat, 1 );

% determine indices of st1 with an st2 spike within owin
ovlp                            = cisi( :, 2 ) >= owin( 1 ) & cisi( :, 2 ) <= owin( 2 );
oidx                            = cisi( ovlp, 1 );
vidx                            = ismember( find( clu == 1 ), oidx );
vec( vidx )                     = 1;

% graphics
if ~graphics 
    return
end

figure
colormap( myjet )

maxval1                         = max( abs( mean( mat1, 2 ) ) );
maxval2                         = max( abs( mean( mat2, 2 ) ) );
maxval                          = max( [ maxval1 maxval2 ] );
clim                            = [ -1 1 ] * 2 * maxval;

subplot( 2, 2, 1 )
imagesc( t, 1 : nspikes1, mat1', clim )
title( sprintf( 'ch%d (%d.%d; %d spikes)', channel1, shankclu1( 1 ), shankclu1( 2 ), nspikes1 ) )
axis xy

subplot( 2, 2, 2 )
imagesc( t, 1 : nspikes1, mat2', clim )
title( sprintf( 'ch%d (%d.%d); %d %d \\muV', channel2, shankclu2( 1 ), shankclu2( 2 ), round( clim( 1 ) ), round( clim( 2 ) ) ) )
axis xy

pidx                            = round( 1 : nspikes1 / nplot : nspikes1 );

subplot( 2, 2, 3 )
hold on
ph                              = plot( t, mat1( :, pidx ) , 'b' );
set( ph, 'color', [ 0.5 0.5 1 ] )
ylims                           = [ ylim; zeros( 1, 2 ) ];
ph                              = plot( t, mean( mat1, 2 ), 'b' );
set( ph, 'LineWidth', 2 )
xlabel( 'Time [ms]' )
ylabel( 'Amplitude [{\mu}V]' )

subplot( 2, 2, 4 )
hold on
if sum( pidx & ~vec( pidx ) ) > 0
    ph                          = plot( t, mat2( :, pidx( ~vec( pidx ) ) ), 'k' );
    set( ph, 'color', [ 1 1 1 ] * 0.5 )
end
if sum( pidx & vec( pidx ) ) > 0
    ph                          = plot( t, mat2( :, pidx( vec( pidx ) ) ) , 'b' );
    set( ph, 'color', [ 0.5 0.5 1 ] )
end
ph                              = plot( t, mean( mat2( :, ~vec ), 2 ), 'k' ...
    , t, mean( mat2( :, vec ), 2 ), 'b' );
set( ph, 'LineWidth', 2 )
ylims( 2, : )                   = ylim;
ylims                           = [ min( ylims( :, 1 ) ) max( ylims( :, 2 ) ) ];
xlabel( 'Time [ms]' )

for i                           = 1 : 4
    subplot( 2, 2, i )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    if i == 3 || i == 4
        ylim( ylims )
    end
end

return

% EOF

filebase = '/Volumes/slab1/mF108/dat/mF108_02/mF108_02';
shankclu1 = [ 4 6 ];                                                        % 'intermediate dendrite'
shankclu2 = [ 3 34 ];                                                       % 'axon'
[ mat2, vec, t, mat1 ] = st2snippets( filebase, shankclu1, shankclu2 );
