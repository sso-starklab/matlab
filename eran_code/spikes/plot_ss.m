% plot_ss                   plot statistics for a single cluster
%
% call                      sst = plot_ss( filebase, n1 )
%
% gets                      filebase    can be a filebase or an sst structure. if a filebase, an s2s
%                                           file is also loaded (if existing) and all sig. CCHs are plotted
%                           n1          cluster numbers; either a scalar or a 2-element vector
%                                           [ shank1 clu1 ]
%
% additional arguments (given as name/value pairs):
%
%                           ilevel      {B}         inclusion level (cluster quality): A, B+, {B}, C
%                           plotmode    {-1}        -1 raw CCHs without patches; 0 with; 1 subtracts predictor
%                           suffix      {'s2s'}     suffix of spikes2spikes file (see spikes2spikes.m)
%                           s2s         { [] }      structure with proper fields (see spikes2spikes.m)
% 
% returns                   the sst structure used
%
% calls                     LoadXml                                                 (blab)
%                           LoadInd                                                 (formats)
%                           ParseArgPairs                                           (general)
%                           separators                                              (graph)
%                           check_cluster_quality, check_mono, plot_spk_waveforms, 
%                                   plot_s2s, wfeatures                             (spikes)
%                           my_spectrum                                             (ssp)
%                           ctsts                                                   (stats)
%                           struct_select                                           (structs)
%
% see also                  spikes2spikes, spikes_stats

% 30-jan-12 ES

% revisions
% 12-jan-13 (1) limited gaussian WN samples to 1000 (otherwise can be
%           untractable, or must compute blockwise, but that is unnecerssary)
%           (2) added trough-peak time, stability p-value
%           (3) should also add partitioning to files (depending on srslen)
% 16-may-13 (1) modified for potentially different list of units in sst/mono
%           (2) suffix added to enable flexible use of *s2s* files
% 16-sep-19 cleaned up; argument and file handling improved
% 17-sep-19 ind supported
% 24-dec-19 (1) s2s input argument option created
%           (2) if *xml file not available, default values are used
%           (3) default values for *xml fields can be user modified
% 16-jan-20 (1) toFlip added, default is 1 since neuroscope origin is top
%               left, MATLAB bottom left
% 05-jul-20 (1) toFlip defaulted to 0

function s = plot_ss( filebase, n1, varargin )

% constants
CCHBS                       = 0.001;            % [s] - assumed (see spikes2spikes.m)
maxRandSpks                 = 1e3;              % for monte carlo WN spectrum 

BARCOLOR                    = [ 0 0 0 ];        % ACH/CCH plots
LW                          = 2;                % width of all line plots
COLOR                       = [ 0.7 0.7 1 ];    % waveform plot
COLORRAND                   = [ 1 0 0 ];        % white noise spectrum plot
CALIBWIDTH                  = 2;                % width of calibration 
CALIBCOLOR                  = [ 0 0.7 0 ];      % color of calibaration 

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
if nargs < 2 || isempty( n1 )
    n1                      = 1;
end
[ ilevel, plotmode, suffix, s2s ...
    , Fs, p2pvoltage, nBits, Amplification ...
    , toFlip ...
    ]                       = ParseArgPairs(...
    { 'ilevel', 'plotmode', 'suffix', 's2s' ...
    , 'Fs', 'p2pvoltage', 'nBits', 'Amplification' ...
    , 'toFlip' ...
    }...
    , { 'B', -4, 's2s', [] ...
    , 20000, 2.45, 16, 192 ...
    , 0 ...
    }...
    , varargin{ : } );


% sst and s2s files
if isa( filebase, 'char' )
    sstfname                = [ filebase '.sst' ];
    if exist( sstfname, 'file' )
        load( sstfname, 'sst', '-mat' )
    end
    s2sfname                = [ filebase '.' suffix ];
    if exist( s2sfname, 'file' )
        load( s2sfname, 's2s', '-mat' )
    end
elseif isa( filebase, 'struct' )
    sst                     = filebase;
    filebase                = sst.filebase;
    if isa( s2s, 'struct' )
        s2sLoaded           = 1;
    else
        s2sLoaded           = 0;
    end
    if ~s2sLoaded && exist( fileparts( filebase ), 'dir' )
        s2sfname                = [ filebase '.' suffix ];
        if exist( s2sfname, 'file' )
            load( s2sfname, 's2s', '-mat' )
            s2sLoaded       = 1;
        end
    end
    if ~s2sLoaded
        fprintf( 'Cannot load s2s file for %s\n', filebase )
    end
end
if ~exist( 'sst', 'var' )
    error( 'input type mismatch' )
end

% xml file
try
    par                     = LoadXml( filebase );
    Fs                      = par.SampleRate;                       % for firing rate stationarity analysis
    p2pvoltage              = par.VoltageRange;                     % amplitude scaling
    nBits                   = par.nBits;                            % amplitude scaling
    Amplification           = par.Amplification;                    % amplitude scaling
catch
    fprintf( 'Cannot load par file for %s, using default values for Fs (%d), p2p (%0.3g), nBits (%d), Amplification (%0.3g)\n' ...
        , filebase, Fs,p2pvoltage, nBits, Amplification  )
end

% determine shank and cluster to focus on
if numel( n1 ) == 1                                                 % single number: index of sst
    shanknum                = sst.shankclu( n1, 1 );
    clunum                  = sst.shankclu( n1, 2 );
elseif numel( n1 ) == 2                                             % two numbers: [ shank clu ]
    shanknum                = n1( 1 );
    clunum                  = n1( 2 );
    n1                      = find( sst.shankclu( :, 1 ) == n1( 1, 1 ) & sst.shankclu( :, 2 ) == n1( 1, 2 ) );
else
    error( 'input size mismatch' )
end
if isempty( n1 ) || n1 > size( sst.nspks, 1 )
    fprintf( 1, 'missing cluster\n' )
    return
end

% determine cluster-specific channels
indfname                    = sprintf( '%s.ind.%d', filebase, shanknum );
if exist( indfname, 'file' )
    ind                     = LoadInd( indfname );
    cluchans                = ind( clunum, : );
else
    cluchans                = [];
end
% determine ordering
nChannels                   = size( sst.mean{ n1 }, 1 );
%ridx                        = 1 : length( par.SpkGrps( shanknum ).Channels );
ridx                        = 1 : nChannels;

%---------------------------------------------------------------%
% find sharp CCHs 

% align s2s and sst 
if exist( 's2s', 'var' ) && ~isempty( s2s ) && ~isequal( sst.shankclu, s2s.shankclu )
    [ ~, i1, i2 ]           = intersect( sst.shankclu, s2s.shankclu, 'rows' );
    i1e                     = false( size( sst.shankclu, 1 ), 1 );
    i1e( i1 )               = 1;
    i2e                     = false( size( s2s.shankclu, 1 ), 1 );
    i2e( i2 )               = 1;
    sst                     = struct_select( sst, i1e, 1 );
    s2s                     = struct_select( s2s, i2e, 1 );
end

% determine cluster quality
gidx                        = check_cluster_quality( sst, ilevel );

% check for mono-synaptic connections
if exist( 's2s', 'var' ) && ~isempty( s2s )
    mono                    = check_mono( s2s, gidx, 'suffix', suffix );
    % find the participating clusters
    a                       = [ mono.pairsExc; mono.pairsInh; mono.pairsSync; mono.pairsDesync ];
    b                       = [ ones( sum( gidx ), 1 ) * n1 find( gidx ) ];
    ahat                    = [ mono.shankclu( a( :, 1 ),  1 : 2 ) mono.shankclu( a( :, 2 ),  1 : 2 ) ];
    bhat                    = [ sst.shankclu( b( :, 1 ),  1 : 2 ) sst.shankclu( b( :, 2 ),  1 : 2 ) ];
    pidx                    = ismember( ahat, bhat, 'rows' ) | ismember( ahat, bhat( :, [ 3 4 1 2 ] ), 'rows' );
    tmp                     = ahat( pidx, : )';
    uidx                    = unique( reshape( tmp( : ), [ 2 sum( pidx ) * 2 ] )', 'rows' );
    uidx                    = setdiff( uidx, sst.shankclu( n1, 1 : 2 ), 'rows' );
    uidx                    = find( ismember( sst.shankclu( :, 1 : 2 ), uidx, 'rows' ) );
end

%---------------------------------------------------------------%
% precompute
% waveform
spk( ridx, :, 1 )           = sst.mean{ n1 };
if isequal( size( sst.sd{ n1 } ), size( sst.mean{ n1 } ) )
    spk( ridx, :, 2 )       = sst.sd{ n1 };
end
scalefactor                 = 1 / 2^nBits * p2pvoltage * 1e6 / Amplification;           % a2d units -> microVolts
[ ~, chan ]                 = min( sum( spk( :, :, 1 )' - sst.max( :, n1 ) * ones( 1, size( spk, 1 ) ) ) );
if ~isempty( cluchans ) && length( cluchans ) >= chan
    chan                    = cluchans( chan );
end
% spectrum
nsamps                      = size( sst.mean{ n1 }, 2 );
win                         = hanning( nsamps );
nfft                        = 2 ^ nextpow2( 2 * length( sst.freqs ) );
gnoise                      = normrnd( 0, 1, nsamps, min( sst.nspks( n1 ), maxRandSpks ) );
powrand                     = mean( my_spectrum( gnoise, nfft, Fs, win, 0, 0 ), 2 );
% ach
cch                         = sst.ach( :, n1 );
MAXLAG                      = floor( length( cch ) / 2 ) * CCHBS;
cchbins                     = [ fliplr( -CCHBS : -CCHBS : -MAXLAG )  0 : CCHBS : MAXLAG ];
cchbins                     = cchbins / CCHBS;                                          % ms
% rate/time
binsize                     = sst.Tsec / size( sst.frateb, 1 ); 
tsec                        = binsize / 2 : binsize : sst.Tsec;

%---------------------------------------------------------------%
% plot

mcolor                      = zeros( 1, 3 );
mcolor( COLOR == max( COLOR ) ) = 1;
newplot

% plot all waveforms (scaled)
subplot( 1, 3, 1 )
plot_spk_waveforms( spk, 0, COLOR, 1, scalefactor, [ 0 1 1 ], toFlip );
title( sprintf( '%d.%d (%d): LR=%0.3g; ID=%0.3g (G:%d)'...
    , shanknum, clunum, sst.nspks( n1 ), sst.Lratio( n1 ), sst.ID( n1 ), gidx( n1 ) ) )
ylabel( 'Amplitude [{\mu}V]' )
xlabel( 'Time [ms]' )
lh                          = line( [ -0.5 -0.5 -0.25 ], [ -50 -100 -100 ] );   % 0.25 ms/50 uV calibration
set( lh, 'COLOR', CALIBCOLOR, 'linewidth', CALIBWIDTH )
axis off

% plot mean waveform for peak channel
subplot( 4, 3, 2 )
[ ~, ~, tpNew ]             = wfeatures( sst.max( :, n1 ), 1 );
xlabel( 'Time [samples]' )
lh                          = line( nsamps * [ 0.125 0.125 0.375 ], [ -25 -75 -75 ] );  % 0.25 ms/50 uV
set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
ylabel( 'Amplitude [{\mu}V]' ),
title( sprintf( 'Max: %d; Geo: %0.3g', chan, sst.geo_com( n1 ) ) )
axis off

% plot spectrum for that channel
subplot( 4, 3, 5 )
line( sst.freqs, sst.spec( :, n1 ), 'color', mcolor, 'Linewidth', LW );
set( gca, 'xscale', 'log', 'xlim', [ 10 Fs/2 ] );
ylabel( 'Power ({\mu}V^2)' )
xlabel( 'Frequency [Hz]' )
title( sprintf( '%d Hz; %0.2g/%0.2g ms; PYR:%d'...
    , round( sst.fmax( n1 ) ), sst.tp( n1 ), tpNew/Fs*1000, sst.pyr( n1 ) ) )

% overlay spectrum of white noise
line( sst.freqs, powrand( 2 : end ) / max( powrand ) * diff( ylim )...
    , 'color', COLORRAND, 'LineStyle', '--', 'LineWidth', LW );
axis off

% plot ACH
subplot( 4, 3, 8 );
bh                          = bar( cchbins, cch, 1 );
xlim( cchbins( [ 1 end ] ) - diff( cchbins( 1 : 2 ) ) * [ -1 1 ] )
ylims                       = [ 0 max( max( cch ), 1 ) * 1.2 ];
ylim( ylims )
set( bh, 'EdgeColor', BARCOLOR, 'FaceColor', BARCOLOR )
ylabel( 'spks/s' )
xlabel( 'Time [ms]' )
set( gca, 'box', 'off', 'tickdir', 'out' )
lh = line( [ -25 -15 ], mean( ylims ) * [ 1 1 ] );                                  % 10 ms calibration
set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
axis off
title( sprintf( '%0.3g spk/s; ISI: %0.3g; COM: %0.3g ms'...
    , sst.frate( n1 ), sst.ISIindex( n1 ), sst.ach_com( n1 ) ) )

% plot frate vs. time
subplot( 4, 3, 11 );
ph                          = line( tsec, sst.frateb( :, n1 ), 'linewidth', LW );
[ ~, ~, stab_pval ]         = ctsts( sst.frateb( :, n1 ), 100 );                    % low p-value - not stable
title( sprintf( 'p(zero slope): %0.2g', stab_pval ) )
xlim( [ 0 sst.Tsec ] )
ylims = ylim;
lh                          = line( [ 10 10 70 ], mean( ylims ) + [ 2 0 0 ] );      % 2 Hz/60 sec calibration
set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
axis off

% plot CCHs
if exist( 's2s', 'var' ) && ~isempty( s2s )
    ncchs                   = length( uidx );
else
    ncchs                   = -1;
end
if ncchs >= 0
    if ncchs > 4
        plotmode            = -3;
    end
    
    for i = 1 : ncchs
        subplot( ncchs, 3, 3 * i ),
        plot_s2s( s2s, [ n1 uidx( i ) ], 'plotmode', plotmode, 'suffix', suffix ); % n1 is the trigger
        xlim( [ -1 1 ] * MAXLAG / 2 * 1000 );
        ylims               = ylim;
        lh                  = line( -MAXLAG / 2 * 1000 + [ 5 5 10 ], mean( ylims ) + [ 50 0 0 ] ); % 5 ms/50 counts calibration
        set( lh, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH )
        separators( 0 );
        tstr                = sprintf( '%d.%d', s2s.shankclu( uidx( i ), 1 ), s2s.shankclu( uidx( i ), 2 ) );
        th                  = text( MAXLAG * 0.25 * 1000, ylims( 1 ) + diff( ylims ) * 0.25, tstr );
        set( th, 'HorizontalAlignment', 'center', 'color', [ 1 0 1 ] )
        title( '' )
        ylabel( '' );
        if i ~= ncchs
            xlabel( '' );
        end
        axis off
    end
    
end

if nargout > 0
    s = sst;
end

return

% EOF

% example:
filebase                    = '/Volumes/Data/phaser3/mouse365/25nov11/dat/es25nov11_4/es25nov11_4';
s2s                         = spikes2spikes( filebase );
sst                         = spikes_stats( filebase );
gidx                        = check_cluster_quality( sst, 'B' ); 
sum( gidx )
figure, plot_ss( filebase, [ 4 9 ] )
figure, plot_ss( filebase, [ 3 13 ] )



