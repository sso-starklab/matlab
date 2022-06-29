% does: 
% (1) determines times of shankclu1 spikes
% (2) determines a single channel for COM of shankclu2 (call this channel2)
% (3) extracts wide-band data from channel2 during shankclu1 spikes (requires halfWidthMS)
% (4) determines for each event whether there is a shankclu2 spike (requires 'causality')
% 
% returns           mat         [microV]; m x n matrix of n events
%                   vec         logical; 1 x n indicator of whether the
%                                   event occurred simultaneouly with a shankclu1 spike
%                   t           [ms]; m x 1 vector of time, from -halfWidthMS : halfWidthMS, at spkFs

% to do:
% scaling to microV
% intersection
% documentation

% 29-mar-21 ES + SSo

function [ mat2, vec, t, mat1 ] = st2snippets( filebase, shankclu1, shankclu2, halfWidthMS )


halfWidthMS                     = 1.5; 
nplot                           = 100;

par                             = LoadXml( filebase );
spkFs                           = par.SampleRate;
nchans                          = par.nChannels;

% (1) get spike times
s1                              = load_spikes( filebase, shankclu1 );
s2                              = load_spikes( filebase, shankclu2 );
st1                             = s1.res;
st2                             = s2.res;
nspikes1                        = length( st1 );
nspikes2                        = length( st2 );

vec                             = false( 1, nspikes1 );

% (2) determine channel for shankclu2
load( [ filebase '.sst' ], '-mat', 'sst' )
idx2                          	= ismember( sst.shankclu, shankclu2, 'rows' );
gcom2                           = sst.geo_com( idx2 );
schannels2                      = sort( par.SpkGrps( shankclu2( 1 ) ).Channels );
channel2                        = schannels2( round( gcom2 ) );
idx1                            = ismember( sst.shankclu, shankclu1, 'rows' );
gcom1                        	= sst.geo_com( idx1 );
schannels1                    	= sort( par.SpkGrps( shankclu1( 1 ) ).Channels );
channel1                        = schannels1( round( gcom1 ) );

% (3) get wide-band data of channel2 during st1 events
halfWidth                       = round( halfWidthMS * spkFs / 1000 );      % [samples]
win                             = ( -halfWidth : halfWidth )';
t                               = win / spkFs * 1000;
nwin                            = length( win );
% % indices of all samples
% midx1                           = ones( nwin, 1 ) * st1';
% midx2                           = win * ones( 1, nspikes1 );
% midx                            = midx1 + midx2;

a                               = memmapfile( [ filebase '.dat' ], 'Format', 'int16' );
nsamples                        = length( a.Data ) / nchans;
if nsamples ~= round( nsamples )
    error( 'problem with dat file' )
end
% go over st1 spikes and extract channel2 data
mat1                            = NaN( nwin, nspikes1 );
mat2                            = NaN( nwin, nspikes1 );
fprintf( 1, 'Extracting waveforms for %d spikes from channels %d and %d ', nspikes1, channel1, channel2 )
for i                           = 1 : nspikes1
    if ~mod( i, 1000 )
        fprintf( 1, '. ' )
    end
    sidx1                       = ( st1( i ) - 1 + win ) * nchans + channel1;
    sidx2                       = ( st1( i ) - 1 + win ) * nchans + channel2;
    ev1                         = a.Data( sidx1 );
    ev2                         = a.Data( sidx2 );
    mat1( :, i )                = ev1;
    mat2( :, i )                = ev2;
end
fprintf( 1, 'done!\n' )

% detrend each waveform
mat1                            = mydetrend( mat1 );
mat2                            = mydetrend( mat2 );

%
figure, 
colormap( myjet )
subplot( 2, 2, 1 ), imagesc( t, 1 : nspikes1, mat1' ), axis xy
subplot( 2, 2, 2 ), imagesc( t, 1 : nspikes1, mat2' ), axis xy
pidx                            = round( 1 : nspikes1 / nplot : nspikes1 );
subplot( 2, 2, 3 ), plot( t, mat1( :, pidx ), 'b' )
subplot( 2, 2, 4 ), plot( t, mat2( :, pidx ), 'b' )
xlabel( 'Time [ms]' )

return

% EOF

filebase = '/Volumes/slab1/mF108/dat/mF108_02/mF108_02';
shankclu1 = [ 4 6 ]; 
shankclu2 = [ 3 34 ];
[ mat2, vec, t, mat1 ] = st2snippets( filebase, shankclu1, shankclu2 );