% make_dualcolor_figures

% 27-may-18 ES

% revision
% 28-jul-18 (1) added P-I curves
%           (2) made a standard make_XXX_figures routine
% 23-feb-19 (1) wide-band data from mP23
% 04-nov-19 (1) organized calls
%           (2) tested temperature - seems nothing was recorded on mP23_03
% 25-dec-19 (1) added temperature (figure )
% 15-jan-20 (1) modified figure 2 to include more wideband and stim channels, less spike channels, 
%               and equal durations for blue and red pulses
%           (2) expanded figure 3 with raster plots
% 28-jan-20 (1) added figure 7 - population PSTHs
% 09-feb-20 (1) modified figure 7 to run from full data structures
%           (2) started adding index / significance information 
% 10-feb-20 (1) figure 7 - added plots with "effect size" (index)
%           (2) figure 7 - removed units with very low baseline rates from
%           the analysis
% 01-may-20 (1) added stats of eta for each set (8 sets)
%           (2) etimated voltages used for each unit

% figure 1: basic device measurements
% figure 2: wide band data (mP23_03)
% figure 3: light responses - PSTHs (mP23_03)
% figure 4: ACH and CCH of selected units (mP23_03, 1.80,86,85)
% figure 5: model simulations
% figure 6: Temperature and spike rate vs. pulse number
% figure 7: summary of population PSTHs


function make_dualcolor_figures( fignums, savef, savetype, outdir )

Nfigs = 10;

nargs = nargin;
if nargs < 1 || isempty( fignums )
    fignums = 0;
end
if ~isempty( setdiff( fignums, [ 0 1 ] ) ) || length( fignums ) < Nfigs
    idx = round( fignums );
    fignums = zeros( Nfigs, 1 );
    fignums( idx ) = 1;
end

if nargs < 2 || isempty( savef )
    savef = 1;
end
if nargs < 3 || isempty( savetype )
    savetype = 'pdf';
end
if nargs < 4 || isempty( outdir )
    outdir = '';
end

% set plotting parameters:
if isempty( outdir ) || ~exist( outdir, 'dir' )
    if ismac
        outdir = '/Users/eranstark/Documents/graphics/dualcolor/data_figures';
        datadir = '/Users/eranstark/Documents/da/dualcolor/';
    end
end

colors = [ 0 0 0.7; 1 0 0 ];
% ----------------------------------------------------------------------
% Figure 1
% Basic device measurements
% ----------------------------------------------------------------------

if fignums( 1 )
    
    % spectrum
    L = load( '~/Documents/slab/projects/ori_noked_DUALCOLOR/part7_paper/measurements_27may18/mat.mat', 'mat' );
    mat = L.mat;
    
    % temperature
    load( '~/Documents/slab/projects/ori_noked_DUALCOLOR/part7_paper/measurements_27may18/mat.mat', 'y' )
    a =	3.91e-3;
    b =	-5.78e-7;
    R0 = 100;
    R = y;
    %R = R0 * ( 1 + a * y + b * y.^2 );
    temp = ((-a/b) - sqrt( (a/b)^2 - 4*(R0-R)/(R0*b) ))/2;
    
    % light power
    % hand-copied from IPCurve DB Rev.1.xlsm:
    mat0 = [
        0.3169      0.16842       0.2552      0.18278      0.21893
        0.62977      0.34038      0.49529       0.3714      0.41408
        1.0345      0.55736      0.84279      0.57761      0.72229
        1.6577      0.92137        1.346      0.86303       1.1378
        2.7748       1.5932       2.1946        1.283       1.8415
        5.4476       3.4875       4.0203       1.8648       3.3499
        39.657       127.76       11.951       2.7045       10.739
        975.77       574.56       609.74       3.9921       526.04
        1903.3       1300.2       1412.2       6.0558       1259.1
        2873.8       1851.4       2255.2       9.4896       1998.2
        3830.6       2397.3       3103.7       15.941       2734.7
        4774.5       2961.5       3936.3       32.351       3459.2
        5729.7         3443         4699       97.307       4181.2
        6696.1         3988       5533.2       183.27       4853.5
        7653.2       4477.9       6288.7       231.51         5465
        8625.4         5068       7049.9        246.2         6212
        9656.9       5560.5       7778.6       245.19       6788.8
        10690       6008.7       8410.7       230.68       7120.1
        11646       6396.4       8981.8       229.27       7628.2
        ]; % red
    I = ( 10 : 5 : 100 )';
    mat0( :, 4 )  = []; % really low
    mat0 = mat0 / 1000; % uW->mW
    mat0ori = mat0;
    
    
    mat1 = [
        0     0.033812     0.054711     0.075614     0.047789
        1.3268     0.081084      0.17872      0.23897      0.10234
        48.983       6.9779       21.673       31.066       4.9104
        98.258       26.041       49.274       72.697       32.617
        147.99       44.887       77.204       114.55        60.42
        198.31       63.903       105.56       156.33       88.546
        248.98        83.03       133.97       198.62          117
        299.32       102.47       162.76        242.1       145.54
        347.59       120.61       190.22       283.93       169.81
        397.66       139.91       218.31       332.83        195.8
        447.82       159.09       246.57       378.51       220.69
        497.43       178.16       275.05       408.69       248.03
        547.18       195.07       303.67       464.09       273.83
        596.7       213.67       332.82       500.29       312.02
        645.87       231.84       361.65       540.69       327.77
        694.77       253.15       390.68       591.56       356.86
        743.4       273.68       416.54       645.64       381.99
        791.45       287.17       446.76       669.08       410.92
        839.47       298.55       474.37       714.92        438.6
        ];
    mat1 = mat1 / 1000; % uW->mW
    mat1ori = mat1;
    
    % get data for the non-Ori devices
    mat0new = zeros( 20, 4 );
    mat1new = mat0new;
    for dli = 2 : 5
        L = load( sprintf( '~/Documents/graphics/dualcolor/amir_12dec19/DL_matlab_data/IP/DL%d/Res.mat', dli ) );
        mat0new( :, dli - 1 ) = table2array( L.Res.Red( :, 4 ) ); % [W]
        mat1new( :, dli - 1 ) = table2array( L.Res.Blue( :, 4 ) ); % [W]
        ilist = table2array( L.Res.Red( :, 1 ) ); % [mA]
    end
    
    % combine all
    [ ~, i1, i2 ] = intersect( I, ilist );
    mat0 = [ mat0ori( i1, : ) mat0new( i2, : ) * 1000 ];
    mat1 = [ mat1ori( i1, [ 1 3 : 5 ] ) mat1new( i2, : ) * 1000 ];
    
    %------------------------------------------------------------------------
    LW = 2;
    
    fig1( 1 ) = figure;
    
    % spectrum (1 device - ON6 or ON3..)
    subplot( 2, 2, 1 )
    ph = plot( mat( :, 1 ), mat( :, 2 ), 'k' );
    set( ph( 1 ), 'linewidth', LW )
    xlim( [ 400 700 ] )
    ylim( [ -0.005 0.6 ] )
    xlabel( 'Wavelength [nm]' ), ylabel( 'Power [AU]' )
    
    % temperature (1 device - ON6)
    subplot( 2, 2, 2 )
    ph = plot ( ( 1 : length( y ) ) * 0.5 / 60, temp, 'k' );
    set( ph( 1 ), 'linewidth', LW )
    lh = line( xlim, min( temp ) * [ 1 1 ], 'linestyle', '--', 'color', [ 0 0 0 ] );
    xlabel( 'Time [s]' ), ylabel( 'T [C]' )
    ylim( [ 23 max( ylim ) ] )
    xlim( xlim )
    
    % PI curve for red LD (4 devices - ON1,3,5,6)
    subplot( 2, 2, 3 ),
    mm = mean( mat0, 2 ); ss = std( mat0, [], 2 ) / sqrt( size( mat0, 2 ) );
    ph = plot( I, mm, '-k', I, mm + ss, '--k', I, mm - ss, '--k' ); %set( gca, 'yscale', 'log' )
    set( ph, 'color', [ 1 0 0 ] )
    set( ph( 1 ), 'linewidth', LW )
    xlabel( 'Current [mA]' ), ylabel( 'Power [mW]' )%ylabel( 'P [\muW]' )
    % add individual lines
    hold on
    ph2 = plot( I, mat0, ':' );
    set( ph2, 'color', [ 1 0 0 ] )
    line( xlim, [ 1 1 ] * 1, 'color', [ 0 0 0 ], 'linestyle', '--' )
    
    % PI curve for blue LD (4 devices - ON1,2,3,5,6)
    subplot( 2, 2, 4 ),
    mm = mean( mat1, 2 ); ss = std( mat1, [], 2 ) / sqrt( size( mat1, 2 ) );
    ph = plot( I, mm, '-k', I, mm + ss, '--k', I, mm - ss, '--b' ); %set( gca, 'yscale', 'log' )
    set( ph, 'color', [ 0 0 0.7 ] )
    set( ph( 1 ), 'linewidth', LW )
    xlabel( 'Current [mA]' ), ylabel( 'Power [mW]' )% ylabel( 'P [\muW]' )
    % add individual lines
    hold on
    ph2 = plot( I, mat1, ':' );
    set( ph2, 'color', [ 0 0 0.7 ] )
    line( xlim, [ 1 1 ] * 0.1, 'color', [ 0 0 0 ], 'linestyle', '--' )
    
    for i = 1 : 4
        subplot( 2, 2, i )
        set( gca, 'tickdir', 'out', 'box', 'off' );
        axis square
    end
    
    fprintf( 'Red: %0.3g+-%0.3g mW\n', mean( mat0( end, : ), 2 ), std( mat0( end, : ), [], 2 ) )
    fprintf( 'Blue: %0.3g+-%0.3g uW\n', mean( mat1( end, : ), 2 ) * 1000, std( mat1( end, : ), [], 2 ) * 1000 )
    %min( I( mean( mat0 >= 1, 2 ) >= 0.5 ) ), min( I( mean( mat1 >= 0.1, 2 ) >= 0.5 ) ),
    
    %-----
    fig = fig1;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG1_part%d', outdir, i );
            fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 2
% Wide band data
% ----------------------------------------------------------------------

if fignums( 2 )
    
    fig2 = [];
    
    % set plotting source and parameters
    filebase = datenum2filebase( { 'mP23', -3 } );
    par = LoadXml( filebase );
    L = load( [ filebase '.sps' ], '-mat' );
    neurochans = L.vote;
    [ chans, targets ] = get_stimchans( par );
    stimchans = chans( ismember( targets, 1 ) );
    tempchan = get_stimchans( par, [], 'temperature' );
    
    neurochans( : ) = 31; % look at this channel: 1.7; the mapping here is 1-based (as opposed to neuroscope, which is 0-based)
    tims = [ 29.17 285 250 % blue PI
        28.30 665 250 % blue: PSI
        33.18 300 1000 ]; % red PI
    savetypePlotter     = NaN;
    ilevelPlotter       = 'C';
    LPF_Plotter         = [];
    HPF_Plotter         = [];
    LPF_stim            = [];
    removeDC            = 1;
    suffix              = 'dat';
    addspikes           = 1;
    addspikewaveforms   = 1;
    poschans            = [];
    %poschans            = tempchan;
    toFlipTraces        = 0;
    spikeshanknums      = [ 1 80; 1 85; 1 86 ]; % P (PYR), S ("SOM"), I (PV)
    chColors            = repmat( [ 1 1 1 ] * 0.5, [ 6 1 ] );
    
    %--------------------------------------------
    % get the transformation matrix
    channels = stimchans;
    
    V2I                 = 0.01;
    [ pathname, filename ] = fileparts( filebase );
    homedir                     = fileparts( fileparts( pathname ) );
    
    % load transformation table
    ichans                      = [];
    tabfile                     = sprintf( '%s/prm/%s.tab', homedir, filename );
    if exist( tabfile, 'file' )
        L                       = load( tabfile, '-mat' );
    else
        fprintf( 1, '%s: No table found for %s in %s...\n', mfname, filename, tabfile )
    end
    if exist( 'L', 'var' ) && all( isfield( L, { 'tab', 'chans' } ) )
        % first, check the data in the file
        [ m, n ]                 = size( L.tab );
        if m == 0 || n == 0 || ( n - 1 ) ~= length( L.chans )
            fprintf( 1, '%s: Table in %s not constructed properly - ignored!!\n', mfname, tabfile )
        else
            % second, transform the units:
            % V command [V] -> I out [A]    divide by 100
            % P [microW] -> P [mW]          divide by 1000
            L.tab( :, 1 )           = L.tab( :, 1 ) * V2I;
            L.tab( :, 2 : n )       = L.tab( :, 2 : n ) * 0.001;
            % third, add a row of zeros at the beginning
            if sum( L.tab( 1, : ) ) ~= 0
                L.tab               = [ zeros( 1, n ); L.tab ];
                m                   = m + 1;
            end
            % fourth, populate individual matrices
            [ ichans, ix ]          = intersect( L.chans, channels );
            imats                   = zeros( m, 2, length( ix ) );
            for i                   = 1 : length( ix )
                imats( :, :, i )    = [ L.tab( :, 1 ) L.tab( :, ix( i ) + 1 ) ];
            end
        end
    end
    imats( :, 1, : ) = imats( :, 1, : ) * 1e3;
    imats( :, 2, : ) = imats( :, 2, : ) * 1e3;
    %--------------------------------------------
    
    % for Figure 1 (wide band, wide range)
    chidx = 1;
    ti = 1 : 2;
    figA = hfoAnalysisPlotter( filebase, tims( ti, : ), neurochans( chidx ), 'stimchans', stimchans( chidx ), 'exchans', poschans ...
        , 'savetype', savetypePlotter, 'figdir', '', 'ilevel', ilevelPlotter, 'LPF', LPF_Plotter, 'LPFstim', LPF_stim, 'HPF', HPF_Plotter ...
        , 'removeDC', removeDC, 'suffix', suffix, 'addspikes', addspikes, 'spikeshanknums', spikeshanknums, 'addspikewaveforms', addspikewaveforms ...
        , 'toFlipTraces', toFlipTraces, 'chColors', chColors, 'imat', imats( :, :, chidx ) );
    
    
    chidx = 2;
    ti = 3;
    figB = hfoAnalysisPlotter( filebase, tims( ti, : ), neurochans( chidx ), 'stimchans', stimchans( chidx ), 'exchans', poschans ...
        , 'savetype', savetypePlotter, 'figdir', '', 'ilevel', ilevelPlotter, 'LPF', LPF_Plotter, 'LPFstim', LPF_stim, 'HPF', HPF_Plotter ...
        , 'removeDC', removeDC, 'suffix', suffix, 'addspikes', addspikes, 'spikeshanknums', spikeshanknums, 'addspikewaveforms', addspikewaveforms ...
        , 'toFlipTraces', toFlipTraces, 'chColors', chColors, 'imat', imats( :, :, chidx ) );
    
    %fig2 = [ figA; figB ];
    
    tims = [ 29.17 385 + 25 - 350 700 % blue PI, 50 ms long
        33.18 711 - 250 700 ]; % red PI, 200 ms long
    spikeshanknums = [ 1 80; 1 86 ]; % P (PYR), S ("SOM"), I (PV)
    
    % new - include three neuronal channels, two light channels
    neurochansAll = [ 21 32 31 ]; % sites 2,5, and 7
    % [ ~, i1 ] = intersect( fliplr( par.SpkGrps( 1 ).Channels + 1 ), neurochansAll ); sort( i1 )
    
    ti = 1;
    figC = hfoAnalysisPlotter( filebase, tims( ti, : ), neurochansAll, 'stimchans', stimchans, 'exchans', poschans ...
        , 'savetype', savetypePlotter, 'figdir', '', 'ilevel', ilevelPlotter, 'LPF', LPF_Plotter, 'LPFstim', LPF_stim, 'HPF', HPF_Plotter ...
        , 'removeDC', removeDC, 'suffix', suffix, 'addspikes', addspikes, 'spikeshanknums', spikeshanknums, 'addspikewaveforms', addspikewaveforms ...
        , 'toFlipTraces', toFlipTraces, 'chColors', chColors, 'imat', imats( :, :, 1 ) );

    ti = 2;
    figD = hfoAnalysisPlotter( filebase, tims( ti, : ), neurochansAll, 'stimchans', stimchans, 'exchans', poschans ...
        , 'savetype', savetypePlotter, 'figdir', '', 'ilevel', ilevelPlotter, 'LPF', LPF_Plotter, 'LPFstim', LPF_stim, 'HPF', HPF_Plotter ...
        , 'removeDC', removeDC, 'suffix', suffix, 'addspikes', addspikes, 'spikeshanknums', spikeshanknums, 'addspikewaveforms', addspikewaveforms ...
        , 'toFlipTraces', toFlipTraces, 'chColors', chColors, 'imat', imats( :, :, 2 ) );
    
    fig2 = [ figA; figB; figC; figD ];

    %-----
    fig = fig2;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG2_part%d', outdir, i );
            figure( fig( i ) );
            pause( 0.2 )
            print( '-dpdf', figname, '-bestfit'),
            %fig_out( fig( i ), 1, figname, savetype ),
        end
    end
    
    %     % deal with the temperature traces:
    %     eegfname = [ filebase '.eeg' ];
    %     par = LoadXml( filebase );
    %     tempchan = get_stimchans( par, [], 'temperature' );
    %     a = memmapfile( eegfname, 'Format', 'int16' );
    %     n = length( a.data );
    %     nchans = par.nChannels;
    %     gainT = 2^16/10;
    %     gainS = 2^16/0.2;
    %     x1 = single( a.data( tempchan : nchans : end ) ) / gainT;
    %     y1 = single( a.data( stimchans( 1 ) : nchans : end ) ) / gainS;
    %     cc = calc_pearson( x1, y1 );
    %     t = ( 1 : length( x1 ) )' / par.lfpSampleRate;
    %     figure
    %     subplot( 2, 1, 1 ), plot( t, x1 ), ylabel( 'Voltage (temperature channel)' ), title( sprintf( 'CC=%0.3g', cc ) )
    %     subplot( 2, 1, 2 ), plot( t, y1 ), xlabel( 'Time [s]' ), ylabel( 'Current [A]' )
    % conclusion - seems that the temperature channel recorded nothing - it was probably in the air during this session
    
    %     % get the wide-band data + scale it
    %     a = memmapfile( [ filebase '.dat' ], 'Format', 'int16' );
    %     i = 1;
    %     tim = tims( i, : );
    %     t0 = min2sec( tim( 1 ) ) + tim( 2 ) / 1000;
    %     win = round( ( t0 + [ 0 tim( 3 ) / 1000 ] ) * par.SampleRate + [ 1 0 ] );
    %     chans = [ neurochans( chidx ) stimchans( 1 ) tempchan ];
    % %    chans = [ tempchan stimchans( 1 ) ]
    %     idx = ones( diff( win ) + 1, 1 ) * ( ( win( 1 ) - 1 ) * nchans + chans ) + ( 0 : nchans : ( diff( win ) * nchans ) )' * ones( 1, length( chans ) );
    %     x = double( a.data( idx ) );
    
    % to translate the D2A to temperature:
    % >> dbstop in parse1channel at 174
    % >> dbstop in parse1channel at 269
    % >> stim = parse1channel( filebase, tempchan );
    % >> chan = 40; source = { 'temp' }; voltagerange = 10;
    % >> dbcont
    % >> medx, gain
    % they turn out to be:
    % medx = -166
    % gain = 2^16/10
    
    
    
end

% ----------------------------------------------------------------------
% Figure 3
% light responses
% ----------------------------------------------------------------------

if fignums( 3 )
    
    % flags for permitting saving *pdf as object oriented file
    plot_multi              = 0;
    plot_rasts              = 1;
    use_patch               = 0;

    
    fig3 = [];
    
    
    % plot PETHs
    filebase                = datenum2filebase( { 'mP23', -3 } );
    [ ~, ~, stim ]          = LoadStims( filebase );
    par                     = LoadXml( filebase );
    
    ilev                    = 'B';
    shanknums               = 1;
    s                       = load_spikes( filebase, shanknums, ilev );
    uflag                   = 1;
    cmp0                    = 'eq';
    
    % select a unit
    shanknums = [ 1 86 ];
    %colors = [ 0 0 0.7; 1 0 0 ];
    fig = [];
    
    for stimnum = 1 : 2
        
        switch stimnum
            case 1
                % select all blue light pulses above 2.68
                %out = stim_select( stim( 1 ), 'types', 'PULSE', 'durs',  [ 0.04 0.06 ], 'vals', [ 2.68 5 ] );
                stimType                = 'PULSE';
                %stimVal = [ 2.68 5 ];
                stimVal                 = [ 0.02785 0.02945 ];
                stimDur                 = [ 0.04 0.06 ];
                stimchans               = 33;
                win                     = [ -325 375 ]; % [ms]
                gaussSD                 = 0.001;
            case 2
                stimType                = 'PULSE';
                %stimVal = [ 2.68 5 ];
                stimVal                 = [ 0.0415 0.05 ]; % three levels: [ 3.8 4.05 ]/100, [ 4.05 4.15 ]/100
                stimDur                 = [ 0.180 0.220 ];
                stimchans               = 37;
                win                     = [ -250 450 ]; % [ms]
                gaussSD                 = 0.01;
        end
        
        if plot_multi
            [ peth, bins, trigs, tims, durs, shankclu, figA ] = multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums, 'stimTypes', stimType ...
                , 'valRange', stimVal, 'durRange', stimDur, 'channels', stimchans, 'uflag', uflag, 'multi', cmp0 ...
                , 'savef', 0 );
        else
            figA = [];
        end
        
        % to plot the rasters and PSTH for each unit:
        spkFs                   = par.SampleRate;
        % (1) get the trigger times (relevant for all units)
        [ ~, tims, durs ]       = get_triggers( filebase, stimchans, uflag, cmp0, [], 'types', stimType, 'durs',  stimDur, 'vals', stimVal );
        % create a window of 325 ms before the trigger until 375 after the trigger
        wins                    = win * spkFs / 1000;
        ntrigs                  = length( tims );
        periods                 = tims * ones( 1, 2 ) + ones( ntrigs, 1 ) * wins;
        nbins                   = diff( wins ) + 1;
        mat                     = sparse( nbins, ntrigs);
        
        % (2) get the spike times (for the specific unit)
        aclu                    = s.map( ismember( s.map( :, 2 : 3 ), shanknums, 'rows' ), 1 );
        res                     = s.res( s.clu == aclu );
        
        % (3) populate the sparse array:
        % (3.1) get the spike times relative to the onset of each period
        [ ~, col, row ]         = inranges( res, periods, 1 );
        % (3.2) actually populate
        nspikes                 = length( col );
        for i                   = 1 : nspikes
            mat( row( i ), col( i ) ) = 1;
        end
        % (3.3) prepare a time vector
        timevec                 = wins( 1 ) : 1 : wins( 2 ); %[samples]
        timevec                 = timevec / spkFs * 1000; % [ms]
        
        % (4) compute PSTH w/ error bars
        % (4.1) compute IFRs
        %plot( firfilt( full( mean( mat, 2 ) ), gwin ) * spkFs )
        %firfilt( full( mean( mat, 2 ) ), gwin )
        gwin                    = makegaussfir( gaussSD, spkFs ); % gaussian kernel w/ SD = 1 ms
        smat                    = firfilt( full( mat ), gwin ) * spkFs;
        psth                    = mean( smat, 2 );
        psth_sem                = calc_sem( smat, 2 );
        
        % plot:
        figA2 = figure;
        
        if plot_rasts
            subplot( 2, 1, 1 )
            %plot_raster( mat( :, 1 : min( size( mat, 2 ), 100 ) ), timevec );
            plot_raster( mat, timevec );
            ylabel( 'Trial' )
            panels = 1 : 2;
        else 
            panels = 2;
        end
        
        subplot( 2, 1, 2 )
        acolor = colors( stimnum, : );
        acolor( acolor == 0 ) = 0.5;
        if use_patch
            patch_band( timevec, psth, psth_sem, colors( stimnum, : ), acolor );
        else
            ph = plot( timevec, psth, 'b', timevec, psth + psth_sem, 'k', timevec, psth - psth_sem, 'k' );
            set( ph( 1 ), 'color', colors( stimnum, : ) );
            set( ph( 2 : 3 ), 'color', acolor );
        end
        ylabel( 'Firing rate [spks/s]' )
        
        for i = panels
            subplot( 2, 1, i )
            ylim( [ 0 max( ylim ) ] )
            alines( [ 0 mean( durs ) * 1000 ], 'x', 'color', colors( stimnum, : ), 'linestyle', '--' );
            set( gca, 'tickdir', 'out', 'box', 'off' )
            xlabel( 'Time [ms]' )
            xlim( win )
        end
        fig = [ fig; figA; figA2 ];
        
    end
    fig3 = fig;
    
%     % out = stim_select( stim( 5 ), 'types', 'PULSE', 'durs',  [ 0.15 0.25 ], 'vals', [ ] )
%     stimType = 'PULSE';
%     stimVal = [ 4.15 5 ]; % or [ 3.8 4.05 ], [ 4.05 4.15 ]
%     stimDur = [ 0.180 0.220 ];
%     stimchans = 37;
%     
%     [ peth, bins, trigs, tims, durs, shankclu, figB ] = multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums, 'stimTypes', stimType ...
%         , 'valRange', stimVal, 'durRange', stimDur, 'channels', stimchans, 'uflag', uflag, 'multi', cmp0 ...
%         , 'savef', 0 );
%     
    %     % plot separately those that do/not undergo silencing
    %     shankclu0 = shankclu( shankclu( :, 3 ) == 0, : );
    %     shankclu1 = shankclu( shankclu( :, 3 ) == 1, : );
    %     shanknums0 = [ shankclu0; shankclu1( [ 5 6 15 ], : ) ];
    %     shanknums1 = [ shankclu0; setdiff( shankclu1, shanknums0, 'rows' ) ];
    %     shanknums0( :, 3 ) = [];
    %     shanknums1( :, 3 ) = [];
    %     [ peth, bins, trigs, tims, durs, shankclu, figC ] = multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums0, 'stimTypes', stimType ...
    %         , 'valRange', stimVal, 'durRange', stimDur, 'channels', stimchans, 'uflag', uflag, 'multi', cmp0 ...
    %         , 'savef', 0 );
    %     [ peth, bins, trigs, tims, durs, shankclu, figD ] = multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums1, 'stimTypes', stimType ...
    %         , 'valRange', stimVal, 'durRange', stimDur, 'channels', stimchans, 'uflag', uflag, 'multi', cmp0 ...
    %         , 'savef', 0 );
    
    % shanknums1 = [ 1 18 ];% apparently increases spiking
    % [ peth, bins, trigs, tims, durs, shankclu, figE ] = multipeth_make( filebase, 'ilevel', ilev, 'data', s, 'shanknums', shanknums1, 'stimTypes', stimType ...
    %         , 'valRange', stimVal, 'durRange', stimDur, 'channels', stimchans, 'uflag', uflag, 'multi', cmp0 ...
    %         , 'savef', 0 );
    
    %fig3 = [ figA; figA2; figB ];
    
    %-----
    fig = fig3;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG3_part%d', outdir, i );
            figure( fig( i ) );
            pause( 0.2 )
            print( '-dpdf', figname, '-bestfit'),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 4
% ACH and CCH of selected units
% ----------------------------------------------------------------------

if fignums( 4 )
    
    fig4 = [];
    
    filebase = datenum2filebase( { 'mP23', -3 } );
    figA = figure; plot_ss( filebase, [ 1 80 ] );
    figB = figure; plot_ss( filebase, [ 1 86 ] );
    figC = figure; plot_ss( filebase, [ 1 85 ] );
    
    figE( 1 ) = figure;
    plot_s2s( filebase, [ 1 86; 1 85 ], -1 );
    figE( 2 ) = figure;
    plot_s2s( filebase, [ 1 86; 1 85 ], -4 );
    xlim( [ -25 25 ] )
    
    figE( 3 ) = figure;
    plot_s2s( filebase, [ 1 86; 1 80 ], -1 );
    figE( 4 ) = figure;
    plot_s2s( filebase, [ 1 86; 1 80 ], -4 );
    xlim( [ -25 25 ] )
    
    figE( 5 ) = figure;
    plot_s2s( filebase, [ 1 80; 1 85 ], -1 );
    figE( 6 ) = figure;
    plot_s2s( filebase, [ 1 80; 1 85 ], -4 );
    xlim( [ -25 25 ] )
    
    fig4 = [ figA; figB; figC; figE( : ) ];
    
    %-----
    fig = fig4;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG4_part%d', outdir, i );
            figure( fig( i ) );
            pause( 0.2 )
            print( '-dpdf', figname, '-bestfit'),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 5
% Model simulation
% ----------------------------------------------------------------------
if fignums( 5 )
    
    [ d1, Pfiber ] = fiber_simulation( 0 );
    fig5 = figure;
    Red_Res_Range   = [6396.388865 11645.73241];
    Blue_Res_Range  = [298.5540018 839.4688884];
    
    subplot( 2, 2, 1 )
    plot(d1*10^3,log10( Pfiber(:,2)*1000 ),'b', d1*10^3,log10(Pfiber(:,1)*1000 ), '--b', d1*10^3,log10(Pfiber(:,3)*1000 ), '--b' )
    title( '450 nm LD' )
    ylim( log10( [0.03 120 ]) );
    xlabel( 'Optical pathway [mm]' )
    ylabel( 'Coupled light [mW]' )
    lin2log( 'y', 10 );
    set( gca, 'box', 'off', 'tickdir', 'out' )
    hold on
    plot( d1 * 1000, log10( ones( 1, length( d1 ) ) * 0.1 ),'--r');
    [ ~, y_index ] = min( abs( d1 * 10^3 - 2.9 ) );
    x_pos = d1( y_index ) * 10^3;
    plot( [x_pos x_pos], ylim, '--k');
    plot( [x_pos x_pos], log10( Blue_Res_Range / 1000 ), 'g' );
    
    subplot( 2, 2, 3 )
    plot(d1*10^3,log10(Pfiber(:,5)*1000 ),'b', d1*10^3,log10(Pfiber(:,4)*1000),'--b', d1*10^3,log10(Pfiber(:,6)*1000 ),'--b')
    title('638 nm LD')
    ylim( log10( [0.03 120 ]) );
    xlabel( 'Optical pathway [mm]' )
    ylabel( 'Coupled light [mW]' )
    lin2log( 'y', 10 );
    set( gca, 'box', 'off', 'tickdir', 'out' )
    hold on
    plot( d1 * 1000, log10( ones( 1, length( d1 ) ) * 1 ),'--r');
    
    [ ~, y_index ] = min( abs( d1 * 10^3 - 0.46 ) );
    x_pos = d1( y_index ) * 10^3;
    plot( [x_pos x_pos], ylim, '--k');
    plot( [x_pos x_pos], log10( Red_Res_Range / 1000 ), 'g' );
    
    
    %-----
    fig = fig5;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG5_part%d', outdir, i );
            figure( fig( i ) );
            pause( 0.2 )
            print( '-dpdf', figname, '-bestfit'),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 6
% Temperature and spike rate vs. pulse number
% ----------------------------------------------------------------------
if fignums( 6 )
    
    %-----
    
    filebase = filebaseLookup( 'mP23', -16 );
    filename = [ filebase '.eeg' ];
    
    chNum = 44;
    chOfIntrestTemp = 39; % temperature
    chOfIntrestCurrent = 32; % current to blue LD, shank 1
    chOfIntrestNeural = 25; % shank 2, site 2 from bottom (note - plot_ss is flipped for this session)
    
    % if y_min >100
    % max_temp_change = ((-a/b) - sqrt( (a/b)^2 - 4*(R0-y_max)/(R0*b) ))/2 -  ((-a/b) - sqrt( (a/b)^2 - 4*(R0-y_min)/(R0*b) ))/2;;
    % else
    %     max_temp_change=NaN; %10.15 in mP23_16
    % end
    % end
    
    % to analyze, we proceed as follows:
    % (1) determine all 50 ms blue light psines on channel 33 unique, define
    % these as "periods", and the mean of each period as "time" (x)
    % (2) count the number of spikes of unit 2.35 during these events, define
    % these as "spike counts" (y1)
    % (3) average the temparture of channel 40 during these events, define
    % these are "mean temperature" (y2)
    % (4) plot y1 and y2 vs. x
    
    par = LoadXml( filebase );
    spkFs = par.SampleRate;
    eegFs = par.lfpSampleRate;
    
    % step 1:
    [ ~, ~, stims ] = LoadStims( filebase );
    stm = stim_select( stims( 1 ), 'types', 'PULSE', 'durs', [ 0.09 0.11 ], 'vals', [ 0.028 0.029 ], 'index', 1 );
    stm = struct_select( stm, stm.index == 1 );
    periods = stm.times; % @ spkFs
    mp = mean( periods, 2 );
    % channels = 33;
    % uflag = 1;
    % multi = 'eq';
    % simOnly = -1;
    % params = { 'types', 'PULSE', 'durs', [ 0.09 0.11 ], 'vals', [ 0.028 0.029 ] };
    % [ ~, ~, ~, ~, stims1 ] = get_triggers( filebase, channels, uflag, multi, simOnly, params )
    % isequal( stims1, stm0 )
    
    % verify that continuous
    %plot( ( mp - mp( 1 ) ) / 20000 / 60, '.b' )
    cp = find( diff( mp ) > 2 * mean( diff( mp ) ), 1, 'last' ); % remove outliers (inter-pulse-interval larger than twice the mean)
    kidx = ( cp + 1 ): size( periods, 1 );
    periods = periods( kidx, : );
    periods0 = periods;
    
%     % sample at same dc backwards in time
%     si              = periods( 1, 1 ) - cumsum( diff( periods( :, 1 ) ) );
%     d               = flipud( diff( periods( 2 : end, : ), [], 2 ) );
%     ei              = si + d;
%     bperiods        = flipud( [ si ei ] );
%     
%     % check that the arbitrary times do not overlap with other stimuli
%     %diff( periods( 2 : end, : ), [], 2 )- diff( nperiods, [], 2 )
%     %diff( periods( 2 : end, 1 ), [], 1 ) - flipud( diff( nperiods( :, 1 ), [], 1 ) )
%     %figure, plot( stims( 1 ).times( :, 1 ), 1, '.b', periods( :, 1 ), 1, '.r', nperiods( :, 1 ), 1, 'ok' )
%     bperiods( 1 : 119, : ) = [];
%     periods        = [ bperiods; periods0 ];
    
    
%     % select "random" periods half-way between the previous ones
%     si              = periods( 1 : end - 1, 1 ) + diff( periods( :, 1 ), [], 1 ) / 2;
%     d               = diff( periods( 1 : end - 1, : ), [], 2 );
%     rperiods        = [ si si + d ];
%     periods         = rperiods;
    
    nperiods        = size( periods, 1 );
    
    % step 2:
    % load the spike train of the unit of interest
    shankclu = [ 2 35 ];
    shanknums = unique( shankclu( :, 1 ) );
    ilevel = 'B';
    spk = load_spikes( filebase, shanknums, ilevel );
    aclu = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu, 'rows' ), 1 );
    idx = spk.clu == aclu;
    %clu = spk.clu( idx );
    res = spk.res( idx ); % @ spkFs
    
    % count the number of spikes in each period:
    [ idx, pidx, out ] = inranges( res, periods );
    [ sc, pnum ] = uhist( pidx );
    counts = zeros( nperiods, 1 );
    counts( pnum ) = sc;
    counts_per_s = counts ./ ( diff( periods, [], 2 ) / spkFs ); % spk/s
    
    % count spikes in a running window
    win = ones( round( nperiods / 10 ), 1 );
    win = win / sum( win );
    scount = firfilt( counts, win ); % spk/period
    fr = scount ./ ( diff( periods, [], 2 ) / spkFs ); % spk/s
    
    % step 3:
    periods1 = resampleranges( periods, eegFs, spkFs ); % @ eegFs
    time = periods1( [ 1 end ] ) / eegFs / 60; % [min]
    [ tempar, t0 ]=LoadData_temp_amir( filename, time, chNum, chOfIntrestTemp, 'temp', eegFs );
    periods2 = periods1 - periods1( 1 ) + 1; % @ eegFs, relative to tempar onset (t0(1))
    %mdur = round( mean( diff( periods2, [], 2 ) ) ) + 1;
    mdur = min( diff( periods2, [], 2 ) ) + 1;
    %periods3 = [ periods2( :, 1 ) periods2( :, 1 ) + mdur - 1 ]; % same length for all periods
    % because mdur is short and nperiods is not large, matrix indexing is faster:
    vec = periods2( :, 1 );
    mat = vec * ones( 1, mdur ) + ones( nperiods, 1 ) * ( 0 : ( mdur - 1 ) );
    tempmat = tempar( mat ); % different periods in different rows
    mtemp = mean( tempmat, 2 );
    
    % step 4: plot
    xmode = 'time';
    dc = sum( diff( periods, [], 2 ) + 1 ) / ( periods( end, 2 ) - periods( 1, 1 ) );
    switch xmode 
        case 'pulses'
            x = ( 1 : nperiods )'; % by 
            xlbl = 'Pulse number';
        case 'time'
            x = ( periods( :, 1 ) - periods( 1, 1 ) ) / spkFs / 60;
            xlbl = 'Time [min]';
    end
    tstr = sprintf( '%d pulses, %0.2g%% dc', size( periods, 1 ), dc * 100 );
    
    figure,
    ph = plot( x, fr, 'b', x, mtemp, 'k' );
    % ph = plot( x, counts_per_s, '.k', x, fr, 'b', x, mtemp, 'k' );
    % set( ph( 1 ), 'color', [ 1 1 1 ] * 0.7 )
    ylabel( 'Firing rate [spks/s] or temperature [°C]' )
    xlabel( xlbl )
    set( gca, 'tickdir', 'out', 'box', 'off' ),
    title( tstr )
    xlim( x( [ 1 end ] ) )
    
    fig6 = figure;
    
    subplot( 2, 2, 1 )
    ph = plot( x, fr, 'b' );
    ylabel( 'Firing rate [spks/s]' )
    xlabel( xlbl )
    set( gca, 'tickdir', 'out', 'box', 'off' ),
    title( tstr )
    xlim( x( [ 1 end ] ) )
    
    subplot( 2, 2, 2 ),
    ph = plot( x, mtemp, 'k' );
    ylabel( 'Temperature [°C]' )
    xlabel( xlbl )
    set( gca, 'tickdir', 'out', 'box', 'off' ),
    title( tstr )
    xlim( x( [ 1 end ] ) )
    
    %-----
    fig = fig6;
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG6_part%d', outdir, i );
            figure( fig( i ) );
            pause( 0.2 )
            print( '-dpdf', figname, '-bestfit'),
        end
    end
    
    
end

% ----------------------------------------------------------------------
% Figure 7
% summary of population PSTHs
% ----------------------------------------------------------------------
if fignums( 7 ) || fignums( 8 )
    
    
    %-----
    if ismac
        % s = make_dualcolor_figures_psths_call; % this calls, eventually, multipeth_make
        %results_filename           = 's.mat'; % 28jan20
        %results_filename           = 's_blue.mat'; % 09feb20
        results_filename           = 'ss.mat'; % 10feb20
        load_fullname           = [ datadir '/' results_filename ];
        %load( load_fullname, 's' );
        load( load_fullname, 'ss' );

        DC_filename             = 'DC.mat'; % 10feb20
        load_fullname           = [ datadir '/' DC_filename ];
        %load( load_fullname, 's' );
        load( load_fullname, 'DC' );

%         % load the data and pad with NaNs for mB142 (last 10 cells not tested with red light)
%        load( load_fullname, 'ss' ); % 09feb20
%         s = ss;
%         clear ss
%         s.bins2 = [ s.bins2; NaN( 10, size( s.bins2, 2 ) ) ];
%         s.psth2 = [ s.psth2; NaN( 10, size( s.psth2, 2 ) ) ];
        
%         % load the list of units with significant effect?? % 09feb20
%         L1              = load( [ datadir '/blue_act.mat' ] );
%         s1( 1 )         = L1.blue_act;
%         L1              = load( [ datadir '/blue_sup.mat' ] );
%         s1( 2 )         = L1.blue_sup;
%         L1              = load( [ datadir '/red_act.mat' ] );
%         s1( 3 )         = L1.red_act;
%         L1              = load( [ datadir '/red_sup.mat' ] );
%         s1( 4 )         = L1.red_sup;
        
        % combine the two structures 
        if isequal( ss.shankclu, DC.shankclu ) && isequal( ss.filename, DC.filename )
            s = ss;
            s.pval = DC.pv; 
            s.fin = DC.fin; 
            s.fout = DC.fout;
            % pval are organized: blue_act, blue_sup, red_act, red_sup
            % fin, fout: blue, red
            eta = ( s.fin - s.fout ) ./ ( s.fin + s.fout ); % light effect index: 1 when max activation, -1 when max silencing
            s.eta = eta;
        end
        
    end
    
    %fig7 = figure;
    
    %multipeth_plot( kpeth, kbins, shankclu );

    nT          = [ -2 2 ];     % hidden in multipeth_make
    fTH         = 0.01;         % firing rate TH [spks/s], same for blue/red, but evaluated separately
    ALFA        = 0.01;
    
    for i = 1 : 2
        % get the data
        switch i
            case 1
                bins = s.bins1;
                peth = s.psth1';
            case 2
                bins = s.bins2;
                peth = s.psth2';
        end
        shankclu = s.shankclu;
        
        % choose only one stim period before and after the onset
        bins        = nanmean( bins, 1 );
        nbins       = length( bins ) - 1;
        kidx        = round( ( nbins - nbins / diff( nT ) * ( diff( nT ) - 1 ) ) + 1 : length( bins ) );
        kpeth       = peth( kidx, : );
        kbins       = bins( kidx );
        
        durs        = ones( size( peth, 2 ), 1 ) * ( bins( end ) - bins( 1 ) ) / diff( nT );
        trigs       = ones( size( durs ) );
        vals        = trigs;
        
        % plot locally
        for ct      = 0 : 1

            vidx    = s.fin( :, i ) >= fTH;
            cidx    = shankclu( :, 3 ) == ct;
            idx     = cidx & vidx;
            
            pp      = kpeth( :, idx );
            sc      = shankclu( idx, : );
            ee      = eta( idx, i );
            % ct 0, i = 1 -> 1
            % ct 1, i = 1 -> 2
            % ct 0, i = 2 -> 3
            % ct 1, i = 2 -> 4
            col     = 2 * ( i - 1 ) + ct + 1;
            pval    = s.pval( :, col );
            pval    = pval( idx, : );
            
            [ ~, sidx ] = sort( ee );
            pval    = pval( sidx, : );
            pp      = pp( :, sidx );
            ee      = ee( sidx, : );
            
            % problem - some units fire very little and thus the eta is
            % NaN; for others, the rate is sufficienly low to have non-zero
            % fout but fin is zero; this could be an effect of the light or
            % due to low rate and a small number of pulses (small sampling
            % duration T).
            % to address this fully, should decide on a minimum sensed
            % baseline lambda, e.g. 0.05 spks/s; and then compute the
            % duration required to sample to see at least some spikes, plus
            % to be able to observe a significant decrease from that;
            % handwaving estimates by Poiss(5/2)<0.01 give about 100
            % seconds of sampling; so, we can choose only units that were
            % tested with at least 100/0.05 pulses (2000) which is very
            % large.
            % in practice, we can do two other things - ignore all units
            % that do not spike at all during one or more periods, or base
            % the process on the fout (e.g. fout must be > TH for both red
            % and blue) - see vidx
            
%             pidx = pval <= ALFA;
%             fprintf( 1, 'Light=%d, ct=%d, eta=%0.3g(%0.3g), %d/%d\n' ...
%                 , i, ct, mean( ee( pidx ) ), calc_sem( ee( pidx ) ), sum( pidx ), length( ee ) )
            
            fig72( col ) = figure;
            %imagesc( kbins, 1 : size( pp, 2 ), scale( pp )' ), axis xy
            [ hh, ah ] = imagescbar( kbins, 1 : size( pp, 2 ), pp, 1 );
            subplot( ah( 1 ) )
            xlabel( 'Time [ms]' )
            ylabel( 'Unit number' )
            colormap( flipud( gray ) )
            lh = alines( [ 0 durs( 1 ) ], 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            
            subplot( ah( 2 ) )
            bh2 = barh( 1 : size( pp, 2 ), ee ); 
            set( bh2, 'FaceColor', colors( ct + 1, : ), 'EdgeColor', colors( ct + 1, : ) );
            xlabel( 'Index' )
            xlim( [ -1 1 ] )
            lh = alines( find( ee( 1 : end - 1 ) < 0 & ee( 2 : end ) > 0 ) + 0.5, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
            set( gca, 'tickdir', 'out', 'box', 'off' ),
            set( gca, 'YTickLabel', '', 'YTick', '' )
            
        end

        % sort by time of peak in PETH and plot using multipeth_plot
        
        % sort separately for PYR and INT
        vidx = s.fin( :, i ) >= fTH;
        kpeth0 = kpeth( :, vidx );
        speth = kpeth0;
        shankclu0 = shankclu( vidx, : );
        sshankclu = shankclu0;
        for ct = 0 : 1
            
%             vidx    = s.fin( :, i ) >= fTH;
            cidx    = shankclu0( :, 3 ) == ct;
%             cidx    = cidx & vidx;            

            cpeth = kpeth0( :, cidx );
            %tidx = true( size( bins ) ); % sort by any bin
            tidx = kbins >= 0 & kbins <= durs( 1 ); % sort by bins within the stim time
            if ct == 0 && i == 1 || ct == 1 && i == 2 % INT and blue light or PYR and red light
                [ ~, maxidx ] = max( cpeth( tidx, : ), [], 1 );
            else % PYR and blue light, or INT and red light
                [ ~, maxidx ] = min( cpeth( tidx, : ), [], 1 ); 
            end
            [ ~, sidx ] = sort( maxidx );
            speth( :, cidx ) = cpeth( :, sidx );
            % also sort the shankclu
            cshankclu = sshankclu( cidx, : );
            sshankclu( cidx, : ) = cshankclu( sidx, : );
        end
        
        % plot externally
        sshankclu( :, 1 ) = 1;
%         vidx = s.fin( :, i ) >= fTH;
        fig = multipeth_plot( speth, kbins, sshankclu ...
            , 'durs', durs( vidx ), 'vals', vals( vidx ), 'trigs', trigs( vidx ), 'barType', 'patch_band' );
        %colormap( myjet )
        fig7( i ) = fig;
        %multipeth_plot( s.psth1', mean( s.bins1, 1 ), s.shankclu, 'barType', 'patch_band' );
        
%         % now correspond with list of blue_act (written 09-feb-20)
%         % 1. generate a unique identifier for each filename
%         uf                  = unique( s.filename );
%         nuf                 = length( uf );
%         ufi                 = ( 1 : nuf )';
%         % 2. generate a vector of unique identifiers for each each 
%         fv0                 = zeros( size( s.filename, 1 ), 1 );
%         fv1                 = zeros( size( s1( 1 ).filename, 1 ), 1 );
%         for k               = 1 : nuf
%             idx0            = ismember( s.filename, uf( k ) );
%             fv0( idx0 )     = ufi( k );
%             idx1            = ismember( s1( 1 ).filename, uf( k ) );
%             fv1( idx1 )     = ufi( k );
%         end
%         % 3. identify the items in s that are in s1( 1 )
%         mat1                = [ fv1 s1( 1 ).shankclu s1( 1 ).ct ];
%         mat0                = [ fv0 s( 1 ).shankclu ];
%         idx                 = ismember( mat0, mat1, 'rows' );
% 
%         fig = multipeth_plot( speth( :, idx ), kbins, shankclu( idx, : ) ...
%             , 'durs', durs( idx ), 'vals', vals( idx ), 'trigs', trigs( idx ), 'barType', 'patch_band' );
%         %colormap( myjet )
        
    end
    
    % summary of statistics:
    % INT:
    for ct = 0 : 1
        fprintf( '\nSummary for ct %d (alpha=%0.2g):\n', ct, ALFA )
        % for all units:
        fprintf( 'All units:\n' )
        disp( sum( s.pval <= ALFA & ( s.shankclu( :, 3 ) == ct ) * ones( 1, 4 ) ) )
        disp( sum( s.shankclu( :, 3 ) == ct * ones( 1, 4 ) ) )
        
        
        % for valid firing rate units:
        fprintf( '\nValid units (fout>=%0.2g spks/s)\n', fTH )
        vidx = [ ( s.fin( :, 1 ) >= fTH ) * [ 1 1 ] ( s.fin( :, 2 ) >= fTH ) * [ 1 1 ] ];
        disp( sum( s.pval <= ALFA & ( s.shankclu( :, 3 ) == ct ) * ones( 1, 4 ) & vidx ) )
        disp( sum( ( s.shankclu( :, 3 ) == ct ) * ones( 1, 4 ) & vidx ) )
        % 36/80 INT increased during blue; 34/73 INT decreased during red; 18/80 PYR decreased during blue; 6/103 PYR increased during red
        
        %-----------------------------------------------------------------
        % compute mean eta for each set (01-may-20 revision):
        nmat        = ( s.shankclu( :, 3 ) == ct ) * ones( 1, 4 ) & vidx;
        vmat        = s.pval <= ALFA & ( s.shankclu( :, 3 ) == ct ) * ones( 1, 4 ) & vidx;
        etamat      = [ eta( :, 1 ) * [ 1 1 ] eta( :, 2 ) * [ 1 1 ] ];
        group       = cell( 4, 1 );
        for col     = 1 : 4
            group{ col } = etamat( vmat( :, col ), col );
            myu = mean( etamat( vmat( :, col ), col ) );
            sig = calc_sem( etamat( vmat( :, col ), col ) );
            nsig = sum( vmat( :, col ) );
            nn = sum( nmat( :, col ) );
            fprintf( 1, 'Col=%d, ct=%d, eta=%0.3g(%0.3g), %d/%d\n' ...
                , col, ct, myu, sig, nsig, nn )
        end
        %-----------------------------------------------------------------
       
        if ct == 0
            cols = [ 1 4 ];
            tstr1 = 'INT eta: act during blue vs. sil during red';
            tstr2 = 'INT with both act during blue, sil during red';
        else
            cols = [ 2 3 ];
            tstr1 = 'PYR eta: sil during blue vs. act during red';
            tstr2 = 'PYR with both sil during blue, act during red';
        end
        pvalMU      = utest( group{ cols( 1 ) }, group{ cols( 2 ) } );
        fprintf( '%s: p=%0.3g (U-test)\n', tstr1, pvalMU );
        
        mat         = s.pval <= ALFA & ( s.shankclu( :, 3 ) == ct ) * ones( 1, 4 ) & vidx;
        nDual       = sum( sum( mat( :, cols ), 2 ) == 2 );
        nDualValid  = sum( sum( vidx( :, cols ), 2 ) == 2 & s.shankclu( :, 3 ) == ct );
        fprintf( '%s: %d/%d\n', tstr2, nDual, nDualValid );
        
        %-----------------------------------------------------------------
        % 01-May-20 revision
        % observed vs. expected bi-directionally modulated units:
        idx = sum( vidx( :, cols ), 2 ) == 2 & s.shankclu( :, 3 ) == ct;
        p1 = sum( s.pval( :, cols( 1 ) ) <= ALFA & idx ) / sum( idx );
        p2 = sum( s.pval( :, cols( 2 ) ) <= ALFA & idx ) / sum( idx );
        
        expected = p1 * p2 * sum( idx );
        observed = nDual;
        ntotal = sum( idx );
        %[ p1 p2 observed expected ntotal ]
        %-----------------------------------------------------------------
        
    end
    
    %-----------------------------------------------------------------
    % compute for each unit, the actual light power used (01-May-20 revision)
    
    % get the voltages used in each session
    info_path               = '/Volumes/MacEvo1/Users/eranstark/Documents/da/dualcolor/';
    info_filename           = 'dst_dal_laser.mat';
    info_fullname           = [ info_path '/' info_filename ];
    load( info_fullname, 'dst' );
    for i = 1 : length( dst )
        [ ~, fname ] = fileparts( dst( i ).filebase );
        dst( i ).fname = fname;
    end
    list1  = { dst.fname };
    list2  = { dst.color };
    nunits = size( s.filename, 1 );
    valBlue = NaN( nunits, 2 );
    valRed = NaN( nunits, 2 );
    anum = NaN( nunits, 1 );
    for i = 1 : nunits
        idx1 = ismember( list1, s.filename( i ) );
        idx2 = ismember( list2, 'blue' );
        idx3 = ismember( list2, 'red' );
        if sum( idx1 & idx2 ) > 0
            valBlue( i, : ) = dst( idx1 & idx2 ).stimVal;
        end
        if sum( idx1 & idx3 ) > 0
            valRed( i, : ) = dst( idx1 & idx3 ).stimVal;
        end
        switch s.filename{ i }( 1 : 4 )
            case 'mB14'
                anum( i ) = 1;
            case 'mP23'
                anum( i ) = 2;
            case 'mDL5'
                anum( i ) = 3;
        end
    end
    
    % get the laser conversion tables
    V2I                         = 0.01;             % CS conversion: from V command [V] to applied I [A] from parseNchannels.m
    imatBlue = cell( 3, 1 );
    imatRed = cell( 3, 1 );
    for j = 2 : 3
        switch j
            case 2
                load_filename2           = 'mP23_03.tab';
            case 3
                load_filename2           = 'mDL5_05.tab';
        end
        load_fullname2           = [ info_path '/' load_filename2 ];
        L = load( load_fullname2, 'tab', '-mat' );
        L.tab( :, 1 )           = L.tab( :, 1 ) * V2I; % parseNchannels
        L.tab( :, 2 : end )       = L.tab( :, 2 : end ) * 0.001; % parseNchannels
        imatBlue{ j } = L.tab( :, [ 1 2 ] );
        imatRed{ j } = L.tab( :, [ 1 4 ] );
    end
    
    % convert the current to power
    powBlue = NaN( nunits, 1 );
    powRed = NaN( nunits, 1 );
    selectionMode = 'max';
    for i = 1 : nunits
        j = anum( i );
        if j == 1
            continue
        end
        switch selectionMode
            case 'mean'
                vBlue = mean( valBlue( i, : ) );
                vRed = mean( valRed( i, :  ) );
            case 'max'
                vBlue = max( valBlue( i, : ) );
                vRed = max( valRed( i, :  ) );
        end
        powBlue( i ) = interp1( imatBlue{ j }( :, 1 ), imatBlue{ j }( :, 2 ), vBlue, 'linear', 'extrap' ); % parse1channel, line 285
        powRed( i ) = interp1( imatRed{ j }( :, 1 ), imatRed{ j }( :, 2 ), vRed, 'linear', 'extrap' );
    end
%     nanmedian( powRed )
%     nanmedian( powBlue )
%     
%     figure
%     subplot( 2, 2, 1 ), pidx = s.shankclu( :, 3 ) == 0 & s.pval( :, 1 ) <= ALFA; plot( powRed( pidx ), s.eta( pidx, 1 ), '.b' ), xlabel( 'Power [mW]' ), ylabel( '{\eta}' ), title( 'PV during blue' ), lsline
%     subplot( 2, 2, 2 ), pidx = s.shankclu( :, 3 ) == 0 & s.pval( :, 4 ) <= ALFA; plot( powRed( pidx ), s.eta( pidx, 2 ), '.r' ), xlabel( 'Power [mW]' ), ylabel( '{\eta}' ), title( 'PV during red' ), lsline
%     subplot( 2, 2, 3 ), pidx = s.shankclu( :, 3 ) == 1 & s.pval( :, 2 ) <= ALFA; plot( powRed( pidx ), s.eta( pidx, 1 ), '.b' ), xlabel( 'Power [mW]' ), ylabel( '{\eta}' ), title( 'PYR during blue' ), lsline
%     subplot( 2, 2, 4 ), pidx = s.shankclu( :, 3 ) == 1 & s.pval( :, 3 ) <= ALFA; plot( powRed( pidx ), s.eta( pidx, 2 ), '.r' ), xlabel( 'Power [mW]' ), ylabel( '{\eta}' ), title( 'PYR during red' ), lsline
    
    %-----------------------------------------------------------------

    
    %-----
    fig = [ fig7 fig72 ];    
    if savef
        for i = 1 : length( fig )
            if fig( i ) == 0 || isempty( fig( i ) )
                continue
            end
            figname = sprintf( '%s/dualcolor_FIG7_part%d', outdir, i );
            figure( fig( i ) );
            pause( 0.2 )
            print( '-dpdf', figname, '-bestfit'),
        end
    end
    
end

return

% EOF


filebase            = datenum2filebase( { 'mP23', -3 } );
Fs                  = 20000;
ilev                = 'B';
shanknums           = 1;

[ ~, ~, stim ]      = LoadStims( filebase );
s                   = load_spikes( filebase, shanknums, ilev );
stimchans           = [ stim( 1 : end - 1 ).chan ];

% select all blue light pulses above 2.68
stimType            = 'PULSE';
stimVal             = [ 2.68 5 ];
stimDur             = [ 0.04 0.06 ];

stimVal             = [ 2.5 3 ];
stimDur             = [ 0.005 0.015 ];
indexVal            = 1;                    % 1 if unique, 0 if overlapping with other channels
stimchan            = 33;
chidx               = ismember( stimchans, stimchan );
out                 = stim_select( stim( chidx ), 'types', stimType, 'durs',  stimDur, 'vals', stimVal, 'index', indexVal );





clu                 = s.clu;
res                 = s.res;
tims                = out.times( :, 1 ); % [sample], onset times
trigs               = ones( size( tims ) );
halfwin             = 50; % [bins]
%binsize             = 41; % [samples]
binsize             = floor( mean( diff( out.times, 1, 2 ) + 1 ) / halfwin * 2 + 1 );
scale               = 'hz';
sdGauss             = 1; % [bins]
shankclu            = s.shankclu( ismember( s.shankclu( :, 1 ), shanknums ), : );
shankclu            = [ 1 85 1; 1 86 0 ];
clucat              = shankclu( :, [ 2 3 ] );

[ peth, bins ] = multipeth( clu, res, trigs, tims...
    , 'binsize', binsize, 'halfwin', halfwin, 'scale', scale...
    , 'sdGauss', sdGauss, 'clucat', clucat );

durs                = out.durs;
vals                = out.vals;
figtitle            = 'some figure';
figname             = '';
savetype            = NaN;
[ fig ahx ] = multipeth_plot( peth, bins, shankclu...
    , 'trigs', trigs, 'durs', durs, 'vals', vals...
    , 'str', figtitle, 'figname', figname, 'savetype', savetype );


uClu                = clucat( :, 1 );
sidx                = ismember( clu, uClu );
Clu                 = clu( sidx );
Res                 = res( sidx );
binsizeCCH          = 20;
halfwinCCH          = 50;
Fs                  = 20000;

scaleCCH            = 'count';
epochs              = out.times;
for i = 1 : 2
    if i == 1
        [ ccg bins ]        = CCG( Res, Clu, binsizeCCH, halfwinCCH, Fs, uClu, scaleCCH );
    else
        [ ccg bins ]        = CCG( Res, Clu, binsizeCCH, halfwinCCH, Fs, uClu, scaleCCH, epochs );
    end
    col = [ 1 1 1 ] * 0.5;
    figure
    subplot( 2, 2, 1 )
    bar( bins, ccg( :, 1, 2 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d x %d CCH', uClu( 1 ), uClu( 2 ) ) ),
    subplot( 2, 2, 2 )
    bar( bins, ccg( :, 2, 2 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d ACH', uClu( 2 ) ) )
    subplot( 2, 2, 3 )
    bar( bins, ccg( :, 1, 1 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d ACH', uClu( 1 ) ) )
    for j = 1 : 3
        subplot( 2, 2, j )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        xlim( [ -1 1 ] * ( halfwinCCH - 1 ) ),
    end
end

%--------
% repeat, only the CCH:
filebase            = datenum2filebase( { 'mP23', -3 } );
Fs                  = 20000;
ilev                = 'B';
shanknums           = 1;

%[ aa, bb, stim ]      = LoadStims( filebase );

padbuffer           = [ -0.01 0.01 ]; % [s]
[ Vals Trigs ]      = LoadVals( filebase );
tidx                = true( size( Trigs ) );
uvals               = dilutesegments( Vals( tidx, 1 : 2 ), 0, max( abs( padbuffer ) ) * Fs );

s                   = load_spikes( filebase, shanknums, ilev );
clu                 = s.clu;
res                 = s.res;

ridx                = inranges( res, uvals ); % remove any spike that is in any segment
res( ridx )         = [];
clu( ridx )         = [];


shankclu            = [ 1 85 1; 1 86 0 ];
clucat              = shankclu( :, [ 2 3 ] );

uClu                = clucat( :, 1 );
sidx                = ismember( clu, uClu );
Clu                 = clu( sidx );
Res                 = res( sidx );
binsizeCCH          = 20;
halfwinCCH          = 50;
Fs                  = 20000;

scaleCCH            = 'count';

[ ccg bins ]        = CCG( Res, Clu, binsizeCCH, halfwinCCH, Fs, uClu, scaleCCH );
col = [ 1 1 1 ] * 0.3;
figure
subplot( 2, 2, 1 )
bar( bins, ccg( :, 1, 2 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d x %d CCH', uClu( 1 ), uClu( 2 ) ) ),
subplot( 2, 2, 2 )
bar( bins, ccg( :, 2, 2 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d ACH', uClu( 2 ) ) )
subplot( 2, 2, 3 )
bar( bins, ccg( :, 1, 1 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d ACH', uClu( 1 ) ) )
subplot( 2, 2, 4 )
bar( bins, ccg( :, 2, 1 ), 1, 'facecolor', col, 'edgecolor', col ), title( sprintf( '%d x %d CCH', uClu( 2 ), uClu( 1 ) ) ),
for j = 1 : 4
    subplot( 2, 2, j )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlim( [ -1 1 ] * ( halfwinCCH - 1 ) ),
end

% now get the theta phase
%L = load( [ filebase '.phs' ], '-mat' );
load( [ filebase '.phs' ], '-mat', 'eegchan' );
par         = LoadXml( filebase );
nchans      = par.nChannels;
d           = memmapfile( [ filebase '.eeg' ], 'Format', 'int16' );
eeg         = d.Data( eegchan : nchans : end );
clear d

% now we have the eeg and the spikes.
% make a quick call to spikePhaseFreq.m:
[ h freqDurs phsBins freqBins b ] = spikePhaseFreq( club, resb, xb...
    , 'M', swsM, 'SD', swsSD, 'phases', phases, 'freqs', hBP, 'powTH', powTH...
    , 'verbose', 0, 'spkFs', spkFs, 'eegFs', Fs );

phases = 20;
[ h durs phsBins freqs b hRate ranges hBP ] = spikePhaseFreq( clu, res, eeg, 'phases', phases, 'freqs', { 2, 40, 20, 'linear' }, 'graphics', 1 );
cidx =  clu == 85;
[ h durs phsBins freqs b hRate ranges hBP ] = spikePhaseFreq( clu( cidx ), res( cidx ), eeg, 'phases', phases, 'freqs', { 2, 40, 20, 'linear' }, 'graphics', 1 );





%------------------------------------------------------------------------
% 14-dec-19

% weights of mice:
mat = [ 35.9 41
    31.2 35.8
    30.9 36.9 ];
mean( diff( mat, [], 2 ) ./ mat( :, 1 ) )

%------------------------------------------------------------------------