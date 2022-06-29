% widebandPlotter
%
% widebandPlotter( filebase, tims, neurochans )
%
% filebase - full filebase
% tims - [ time1 time2 dur ].   time1: MM.SS
%                               time2: TTTT
%                               dur:   TTTT
%           i.e. 13.02 340 400 is minute 13, second 2, 340 ms. duration 400 ms.
% neurochans - a selection of neuronal channels, intended to be e.g. one/shank
%                   (up to six supported with different colors)
%
% optional: (argument name/value pairs):
% 
%
% calls         
% ParseArgPairs, LoadXml
% get_stimchans, get_merged_filenum, get_probe, load_spikes
% makegausslpfir, makefir, min2sec, firfilt, scale, uprobeV2P
% parseSchmitt, dehighpass, plotTraces, calibration, clipmat, plot_raster, alines
% fig_out
  
% 02-dec-20 ES + HS based on hfpAnalysisPlotter

% revisions
% 03-dec-20 modified chcolors (weird bug, ad-hoc solution)

function [ fig, splot ] = widebandPlotter( filebase, tims, neurochans, varargin )

% initialize output
fig                             = [];
splot                           = [];

% constants
maxChannels                     = 12;
cellcolors                      = [ 0 0 0.7; 1 0 0 ];                       % INT, PYR
chcolors                        = [ 1 0 0; 0 0.5 0; 0 0 1; 0 0 0; 1 0 1; 0 0.5 1 ]; % 6 shanks
spikeSamples                    = -15 : 16;                                 % for plotting wide-band spike traces

calib                           = [ 0.05 0.5 ];                             % [s, mV]
calibSprobe                     = [ 0.05 -0.5 ];                            % [s, mA]
calibSuprobe                    = [ 0.05 2 ];                               % [s, uW]
%TH2 = 1; % 1 mA, for extra channels

SPIKE_LENGTH                    = 1;
SPIKE_WIDTH                     = 2;
basedur                         = 0.005;                                    % [s]

% arguments
nargs                           = nargin;
if nargs < 3 || isempty( filebase ) || isempty( tims ) || isempty( neurochans )
    return
end
[ stimchans, exchans, savetype, figdir, digitize, addspikes, addspikewaveforms, ilevel...
    , fixedRaster, plotStim, nSD, LPF, LPFstim, LPFex, RemoveDC, RemovePoly, scaling, TH, TH2, toSmoothStim, suffix ...
    , uprobeStim, uprobeStimRescale, uprobeStimDeslope, spikeshanknums ...
    , chColors, imat ...
    , BPF, npoles, filtSelect ...
    , HPF, calcCSD, padBuffer, filtMode, toFlipTraces, flipStimPolarity, verbose ...
    ] = ParseArgPairs(...
    { 'stimchans', 'exchans', 'savetype', 'figdir', 'digitize', 'addspikes', 'addspikewaveforms', 'ilevel'...
    , 'fixedRaster', 'plotStim', 'nSD', 'LPF', 'LPFstim', 'LPFex', 'RemoveDC', 'RemovePoly', 'scaling', 'TH', 'TH2', 'toSmoothStim', 'suffix' ...
    , 'uprobeStim', 'uprobeStimRescale ', 'uprobeStimDeslope', 'spikeshanknums' ...
    , 'chColors', 'imat' ...
    , 'BPF', 'npoles', 'filtSelect' ...
    , 'HPF', 'calcCSD', 'padBuffer', 'filtMode', 'toFlipTraces', 'flipStimPolarity', 'verbose' }...
    , { [], [], 'pdf', pwd, 0, 1, 0, 'B'...
    , 1, 1, 3, 6000, [], [], 0, 0, 1, NaN, 1, 1, 'dat' ...
    , 0, 0, 1, [] ...
    , [], [] ...
    , [], 2, [] ...
    , [], 0, [ -0.01 0.01 ], 'median', 1, 0, 0 }...
    , varargin{ : } );

if length( neurochans ) > maxChannels || length( stimchans ) > maxChannels
    error( 'Max channels exceeded' )
end

if isempty( stimchans )
    plotStim                        = 0;
end

if ~isempty( chColors ) && size( chColors, 1 ) >= length( neurochans ) && size( chColors, 2 ) == 3
    chcolors                        = chColors;
end
chcolors                            = repmat( chcolors, [ 2 1 ] );
stmcolors                           = chcolors; % 6 stim channels
stmcolors( stmcolors == 0 )         = 0.3;                                  % lighter colors
excolors                            = [ 0 0 0 ];

if uprobeStim
    calibS                          = calibSuprobe;
    calib_XY                        = complex( 0, [ 0.4 0.1 ] );
else
    calibS                          = calibSprobe;
    calib_XY                        = [];
end

% start working
chans                               = [ neurochans stimchans exchans ];
nidx                                = find( ismember( chans, neurochans ) );
sidx                                = find( ismember( chans, stimchans ) );
xidx                                = find( ismember( chans, exchans ) );

% check filtSelect
if length( filtSelect ) ~= length( nidx ) 
    filtSelect                  = [];
end
if isempty( filtSelect )
    filtSelect                  = true( 1, length( nidx ) );
end
filtSelect                      = logical( filtSelect );

% load parameters from par
[ ~, filename, extname ]        = fileparts( filebase );
filename                        = [ filename extname ];
figname                         = [ figdir '/' filename '_examples' ];
par                             = LoadXml( filebase );
fwinsize                        = par.SampleRate * 0.001;                   % for stimulus digitization
fwin                            = ones( fwinsize, 1 ) / fwinsize;
medwinsize                      = par.SampleRate * 0.002;                   % for stimulus smoothing
scalefactor                     = 1 / 2^par.nBits * par.VoltageRange * 1e3 / par.Amplification; % a2d -> mV
try
    [ ~, ~, stimvranges ]       = get_stimchans( par );
    scalefactor2                = 1 / 2^par.nBits * stimvranges( 1 ) * 1e3; % [mA]
catch
    scalefactor2                = scalefactor;
end
if flipStimPolarity
    scalefactor2                = scalefactor2 * -1;
end

datfname                        = [ filebase '.' suffix ];
mfilebase                       = get_merged_filenum( filebase );
if isempty( chans )
    spsfile                     = [ mfilebase '.sps' ];
    sps                         = load( spsfile, '-mat' );
    chans                       = sps.stats;
end
nchans                          = par.nChannels;

% load the spikes
probe                           = get_probe( par );
[ ~, shanknums ]                = find( ismember( probe, neurochans ) );
if isempty( spikeshanknums )
    spikeshanknums              = shanknums;
end
if addspikes || addspikewaveforms
    s                           = load_spikes( filebase, spikeshanknums, ilevel );
end

% set colors
chcolors                        = chcolors( [ shanknums; 1 : length( sidx ) ], : );
stmcolors                       = stmcolors( [ shanknums; 1 : length( sidx ) ], : );

% compute filter coefficients
if isnan( LPF )
    gwin                        = 1;
else
    gwin                        = makegausslpfir( LPF, par.SampleRate );
end
if isnan( LPFstim )
    gwinStim                    = 1;
else
    gwinStim                    = makegausslpfir( LPFstim, par.SampleRate );
end
if isnan( LPFex )
    gwinEx                      = 1;
else
    gwinEx                      = makegausslpfir( LPFex, par.SampleRate );
end
Fs                              = par.SampleRate;
if ~isnan( BPF )
    [ bcoefs, acoefs ]          = butter( npoles, BPF / Fs * 2, 'bandpass' );
end
if ~isnan( HPF )
    switch filtMode
        case 'median'
            hwin                = floor( Fs ./ HPF( 1 ) / 2 ) * 2 + 1;
        case 'fir'
            hwin                = makefir( [ HPF( 1 ) NaN ], Fs, [], 'high' );
        case 'dog'
            hwin                = makegausslpfir( HPF( 1 ), Fs, 6 );           
        case 'iir'
    end
else
    hwin                        = [];
end
if ~isempty( HPF ) || ~isempty( BPF )
    padBuffer                   = padBuffer * par.SampleRate;
end

% get a pointer to the data on disk
a                               = memmapfile( datfname, 'Format', 'int16' );
n                               = length( a.data );

% plot a figure for each row of tims
fig                             = zeros( size( tims, 1 ), 1 );
for i                           = 1 : size( tims, 1 )
    
    %---------------------------------------------------------------------
    % (1) get all wide-band data
    tim                         = tims( i, : );
    t0                          = min2sec( tim( 1 ) ) + tim( 2 ) / 1000;
    win                         = round( ( t0 + [ 0 tim( 3 ) / 1000 ] ) * par.SampleRate + [ 1 0 ] );
    tx                          = ( win( 1 ) : win( 2 ) )' / Fs;
    if ~isempty( HPF ) || ~isempty( BPF )
        win                     = win + padBuffer;
    end
    
    idx                         = ones( diff( win ) + 1, 1 ) * ( ( win( 1 ) - 1 ) * nchans + chans ) + ( 0 : nchans : ( diff( win ) * nchans ) )' * ones( 1, length( chans ) );
    x                           = double( a.data( idx ) );
    
    %---------------------------------------------------------------------
    % (2) process the neuronal continuous data
    x( :, nidx )                = x( :, nidx ) * scalefactor;              % scales
    
    % (2.1) iir BPF
    if ~isempty( BPF )               
        nfidx                   = nidx( filtSelect );
        xtmp                    = x( :, nfidx );
        xtmpf                   = filtfilt( bcoefs, acoefs, xtmp );
        x( :, nfidx )          	= xtmpf;
    end
    
    % (2.2) fir HPF
    if ~isempty( HPF )
        nfidx                   = nidx( filtSelect );
        tmp                     = x( :, nidx );
        switch filtMode
            case 'median'
                lo              = medfilt1( tmp, hwin );
                x( :, nfidx )   = tmp - lo;
            case { 'dog', 'fir' }
                lo              = firfilt( tmp, hwin );
                x( :, nfidx )   = tmp - lo;
        end
    end
    
    % remove the padding
    if ~isempty( HPF ) || ~isempty( BPF )
        x( 1 : abs( padBuffer( 1 ) ), : ) = [];
        x( end - padBuffer( 2 ) + 1 : end, : ) = [];
        win                     = win - padBuffer;
    end
    
    % (2.3) remove DC
    if RemoveDC
        x( :, nidx )            = bsxfun( @minus, x( :, nidx ), mean( x( :, nidx ), 1 ) );
    end
    
    % (2.4) remove polynomial fit
    if RemovePoly
        for j                   = 1 : length( nidx )
            xx                  = x( :, nidx( j ) );
            tt                  = ( 1 : length( xx ) )'; 
            pp                  = polyfit( tt, xx, RemovePoly ); 
            rr                  = polyval( pp, tt );
            xxhat               = xx - rr;
            x( :, nidx( j ) )   = xxhat;
        end
    end
    
    % (2.5) compute the CSD
    if calcCSD
        xtmp                    = x( :, nidx );
        xtmp                    = [ xtmp( :, 1 ) xtmp xtmp( :, end ) ];
        xcsd                    = 2 * xtmp( :, 2 : end - 1 ) - xtmp( :, 1 : end - 2 ) - xtmp( :, 3 : end  );
        x( :, nidx )            = xcsd;
    end
    
    %---------------------------------------------------------------------
    % (3) process the non-neuronal continuous data
    
    % (3.1) process the stim data
    x( :, sidx )                = x( :, sidx ) * scalefactor2;
    x( :, sidx )                = firfilt( x( :, sidx ), gwinStim );           % smooth the stim part
    
    if ~isempty( sidx ) && uprobeStim
        x0                      = x( :, sidx ) / 1000; % mV->V
        if all( sum( x0 < 0 ) / size( x0 , 1 ) ) > 0.3 && uprobeStimRescale % acute - phase reconstruction was -1:1...
            x0                  = scale( x0 ) * 4;
        end
        if uprobeStimDeslope
            for col             = 1 : size( x0, 2 )
                if x0( 1, col ) ~= x0( end, col )
                    slp         = x0( 1, col ) : diff( x0( [ 1 end ], col ) ) / ( size( x0, 1 ) - 1 ) : x0( end, col );
                    x0( :, col ) = x0( :, col ) - slp( : );
                end
            end
        end
        p0                      = uprobeV2P( x0 * 0.83 ) / 1000; % [micro-Watts]
        x( :, sidx )            = p0;
    end
    
    % (3.2) process the external channels
    x( :, xidx )                = x( :, xidx ) * scalefactor;
    x( :, xidx )                = firfilt( x( :, xidx ), gwinEx );
    
    % (3.3) smooth (+digitize) the stim channels (+determine onset/offset, maxval)
    sval                        = NaN * ones( 1, length( sidx ) ); % max value of the stim
    tstim                       = NaN * ones( length( sidx ), 2 ); % stim times (onset/offset)
    txstim                      = NaN * ones( length( xidx ), 2 );
    xstim                       = x( :, [ sidx xidx ] );
    if toSmoothStim
        xstimf                  = firfilt( medfilt1( fliplr( xstim ), medwinsize ), fwin ); % median filter  +smooth
    else
        xstimf                  = xstim;
    end
    if ~isempty( imat )
        xstimf                  = interp1( imat( :, 1 ), imat( :, 2 ), xstimf, 'linear', 'extrap' );    % apply reverse non-linear transform to P [mW]
    end
    bidx                        = 1 : ( basedur * par.SampleRate ); 
    xstimf                      = bsxfun( @minus, xstimf, mean( xstimf( bidx, : ) ) ); % remove mean of beginning
    if isnan( TH )
        sTH                     = std( xstimf( bidx, : ) ) * nSD; % detect using a 3SD TH
    else
        sTH                     = TH * ones( size( sidx ) );
    end
    for j                       = 1 : length( sidx )
        stime                   = parseSchmitt( xstimf( :, j ), sTH( j ) * 2, sTH( j ) );
        if size( stime, 1 ) > 1
            fprintf( '%s: IGNORING %d events\n', upper( mfilename ), size( stime, 1 ) - 1 )
            stime               = stime( 1, : );
        end
        if ~isempty( stime )
            tstim( j, : )       = tx( stime )';
            idx                 = false( size( xstimf, 1 ), 1 );
            idx( stime( 1 ) : stime( 2 ) ) = 1;
            sval( j )           = max( xstimf( idx, j ) );
            if digitize
                xstimf( idx, j ) = sval( j );
                xstimf( ~idx, j ) = 0;
            end
        end
    end
    
    % (3.4) digitize the extra channels
    if digitize
        xtmp2                   = x( :, xidx );
        xtmp2                   = firfilt( dehighpass( xtmp2 ), fwin ); % in case this is a repetitive trigger
        xtmp2( xtmp2 < TH2 )    = 0; 
        for j                   = 1 : length( exchans ) 
            idx                 = xtmp2( :, j ) >= ( TH2 );
            xtmp2( idx, j )     = mean( xtmp2( idx, j ) ); 
        end
        for j                   = 1 : length( xidx )
            stmidx              = find( diff( xtmp2( :, j ) ) );
            if length( stmidx ) <= 2
                continue
            end
            txstim( j, : )      = tx( stmidx )';
        end
        
    end
    x( :, [ sidx xidx ] )       = [];

    % (2.6) smooth the neuronal part
    x                           = firfilt( x, gwin );
    
    %---------------------------------------------------------------------
    % (4) process spikes
    
    if addspikes || addspikewaveforms
        % (4.1) get the corresponding spikes
        idx                     = inrange( s.res, win );
        res                     = round( s.res( idx )- win( 1 ) + 1 );
        clu                     = s.clu( idx );
        uclu                    = unique( clu );
        rast                    = sparse( res, clu, 1, round( diff( win ) + 1 ), max( s.map( :, 1 ) ) );
        if fixedRaster
            uidx                = find( ismember( 1 : max( s.map( :, 1 ) ), s.map( :, 1 ) ) );
        else
            uidx                = uclu;
        end
        rast                    = rast( :, uidx );
        celltypes               = s.shankclu( ismember( s.map( :, 1 ), uidx ), 3 );
        cellshanks              = s.shankclu( ismember( s.map( :, 1 ), uidx ), 1 );
        borders                 = find( diff( cellshanks ) ) + 0.5;
        
        rast                    = fliplr( rast );
        borders                 = length( uidx ) - borders + 1;
        
        % (4.2) divide into PYR and INT
        celltypes               = flipud( celltypes );
        irast                   = rast;
        prast                   = rast;
        irast( :, celltypes == 1 ) = 0;
        prast( :, celltypes == 0 ) = 0;
        borders                 = ( borders - 0.5 ) / 0.5 * SPIKE_LENGTH + 0.5;
        
        splot( i ).clu          = clu;
        splot( i ).res          = res;
        splot( i ).irast        = irast;
        splot( i ).prast        = prast;
    else
        splot( i ).clu          = [];
        splot( i ).res          = [];
        splot( i ).irast        = [];
        splot( i ).prast        = [];
    end
    
    %---------------------------------------------------------------------
    % (5) plot
    fig( i )                    = figure;
    ah( 1 )                     = gca;
    
    % (5.1) plot wide band
    if toFlipTraces
        xplot                   = fliplr( x ); % first on top
    else
        xplot                   = x;
    end
    [ pp, ph ]                  = plotTraces( tx, xplot, size( x, 2 ), scaling, [ 0 0  ] );
    calibration( calib, calib_XY );
    for j                       = nidx
        set( ph( j ), 'color', chcolors( length( nidx ) - j + 1, : ) );
    end
    splot( i ).time             = tx;
    splot( i ).neuronal         = pp;
    splot( i ).neuroCalib       = { calib, calib_XY };
    
    % (5.2) plot spike waveforms
    if addspikewaveforms
        hold on
        iidx                    = [];
        pidx                    = [];
        ucellshanks             = unique( cellshanks );
        for k                   = 1 : length( ucellshanks )
            ashank              = ucellshanks( k );
            cidx                = flipud( cellshanks ) == ashank; % cell indices
            if toFlipTraces
                chidx           = flipud( shanknums ) == ashank; % channel indices
            else
                chidx           = shanknums == ashank;
            end
            ppK                 = pp( :, chidx );
            xlims               = [ 1 length( tx ) ];
            
            [ spk_times, ~ ]    = find( irast( :, cidx ) );
            for j               = 1 : length( spk_times )
                idx             = clipmat( spk_times( j ) + spikeSamples, xlims );
                iidx            = [ iidx; idx ];
                phj             = plot( tx( idx ), ppK( idx, : ) );
                set( phj, 'color', cellcolors( 1, : ) );
            end
            [ spk_times, ~ ]    = find( prast( :, cidx ) );
            for j               = 1 : length( spk_times )
                idx             = clipmat( spk_times( j ) + spikeSamples, xlims );
                pidx            = [ pidx; idx ];
                phj             = plot( tx( idx ), ppK( idx, : ) );
                set( phj, 'color', cellcolors( 2, : ) );
            end
        end
        splot( i ).iidx         = iidx;
        splot( i ).pidx         = pidx;
    else
        splot( i ).iidx         = [];
        splot( i ).pidx         = [];
    end
    
    % (5.3) plot spike times
    if addspikes
        pos                     = get( ah( 1 ), 'position' ); 
        set( ah( 1 ), 'position', [ pos( 1 ) pos( 2 ) + pos( 4 ) / 2 pos( 3 ) pos( 4 ) / 2 ] ); 
        axis off
        ah( 2 )                 = axes( 'position', [ pos( 1 : 3 ) pos( 4 ) / 2 ] );
        plot_raster( irast, tx, 1, SPIKE_LENGTH, SPIKE_WIDTH, cellcolors( 1, : ) );
        plot_raster( prast, tx, 1, SPIKE_LENGTH, SPIKE_WIDTH, cellcolors( 2, : ) );
        alines( borders, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
        axis off
    end
    
    % (5.4) plot stim times
    if digitize
        for k                   = 1 : length( ah )
            subplot( ah( k ) )
            for j               = 1 : length( sidx )
                alines( tstim( length( sidx ) - j + 1, : ), 'x', 'color', stmcolors( j, : ), 'linestyle', '--' );
            end
        end
        for j                   = 1 : length( xidx )
            alines( txstim( length( sidx ) - j + 1, : ), 'x', 'color', excolors( j, : ), 'linestyle', '--' );
        end
    end
    
    % (5.5) plot stim waveforms
    if plotStim
        if addspikes
            pos                 = get( ah( 2 ), 'position' );
            set( ah( 2 ), 'position', [ pos( 1 : 3 ) pos( 4 ) / 2 ] )
            ah( 3 )             = axes( 'position', [ pos( 1 ) pos( 2 ) + pos( 4 ) / 2 pos( 3 ) pos( 4 ) / 2 ] );
        else
            pos                 = get( ah( 1 ), 'position' );
            set( ah( 1 ), 'position', [ pos( 1 ) pos( 2 ) + pos( 4 ) / 4 pos( 3 ) pos( 4 ) * 3 / 4 ] );
            axis off
            ah( 3 )             = axes( 'position', [ pos( 1 : 3 ) pos( 4 ) / 4 ] );
        end
        if toFlipTraces
            xplot               = xstimf;
        else
            xplot               = fliplr( xstimf );
        end
        [ ~, ph ]               = plotTraces( tx, xplot, size( xstim, 2 ), 1, [ 0 0  ] );
        for j                   = 1 : length( sidx )
            for k               = 1 : length( ah )
                if ah( k ) == 0 || k == 2 && ~addspikes
                    continue
                end
                subplot( ah( k ) )
                alines( tstim( length( sidx ) - j + 1 , : ), 'x', 'color', stmcolors( j, : ), 'linestyle', '--' );
                set( ph( length( sidx ) - j + 1 ), 'color', stmcolors( j, : ) );
            end
        end
        calibration( calibS, calib_XY );
        axis off
        splot( i ).stim         = xplot;
        splot( i ).stimCalib    = { calibS, calib_XY };
    else
        splot( i ).stim         = [];
        splot( i ).stimCalib    = { [], [] };
    end
    splot( i ).axes             = ah;
    
    % add title
    subplot( ah( 1 ) )
    tstr                        = sprintf( '%s, (%03g.%d; %d ms)%s', filename...
        , tim( 1 ), round( tim( 2 ) ), round( tim( 3 ) )...
        , num2str( sval( end : -1 : 1 ) ) );
    title( tstr )
    splot( i ).titlestr         = tstr;
    
    if verbose
        for ii                  = 1 : length( nidx )
            fprintf( 'Ch%d: range: %0.3g mV\n', chans( nidx( ii ) ), range( x( :, nidx( ii ) ) ) )
        end
    end
    
end

clear a

% save
if ~isempty( savetype ) && ( isa( savetype, 'cell' ) || ~all( isnan( savetype ) ) ) && ~isempty( figname )
    if ~isa( savetype, 'cell' )
        savetype                = { savetype };
    end
    for j                       = 1 : length( savetype )
        for i                   = 1 : length( fig )
            if fig( i ) == 0
                continue
            end
            fignameI            = [ figname '.part' num2str( i ) '_' suffix '.' savetype{ j } ];
            fig_out( fig( i ), 1, fignameI, savetype{ j } );
            fprintf( 'Saving %s...\n', fignameI )
        end
    end
end

return

% EOF

