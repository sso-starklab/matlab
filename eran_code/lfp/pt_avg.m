% pt_avg                pulse-triggered averaging (spectrograms, CSD, MUA, ...)
%
% call:                 pt_avg( filebase, channels, triggers )
%
% arguments:
% 
%   filebase:       either a filebase, a par structure, or a cell array doublet { adate, fnum }
%   channels:       in source file; if negative, gets the channels in the ElecGp
%                       from the xml file
%   triggers:       time in file, in sec
%
% optional arguments (given as argument/value pairs):
%
%   window:         {[ -0.005 0.005 ]} in sec, for averaging
%   tempSD:         {0.0001} in sec, for CSD temporal smoothing (Gaussian)
%   spatBin:        {1}; side band, in channels, for CSD spatial smoothing (triangular)
%   suffix:         {'eeg'}, 'dat' (or anything else binary); 'spc' plot eeg and csd
%                           + their spectrograms; 'wlt' plot wavelet spectrograms
%   graphics:       { 1 }; if a full path, also saves the plot at that location
%   scalefactor:    { 0.12207/1000}, from A2DU to mV (assumes 1K gain, 16 bit/8V digitization)
%   specChans:      {chans}; which channels to average spectrum over 
%
%   whitenFlag:     {0} whiten signals (each channel separately); relevant for
%                        spectral analysis only
%   normalizeSWS:   {100} normalize by mean and SD during 100 random no-stim SWS epochs
%   savetype:       {'png'}
%   refchan:        {[]}; if a valid channel number, subtracts that from data
%   recomputeBL:    {0}, relevant only for 'wlt' spectrograms
%                   (1):recompute; 
%                   (0):recompute only if sample size mismatches
%                   (-2):recompute only if no file    
%
% returns:
%   avgcsd/lfp      time averages (2D matices); time x channels
%   tim             vector, [ms]
%   xcsd/lfphat     full data, 3D arrays (time x channels x instances)
%   too, foo, yoo   time/frequency vectors and spectrogram
%   yooCSD          spectrogram of the CSD
%
% calls:            ParseArgPairs, gausskernel, verb, replacetok, num3str   (general)
%                   LoadXml, WhitenSignal, mtcsg                            (blab)
%                   readbin                                                 (formats)
%                   makefir, firfilt, getWavelet, plotSpectrogram           (ssp)
%                   plotTraces, computeBaseline                             (lfp)
%                   outliers                                                (stats)
%                   imupsample, calibration, fig_out                        (graph)

% 02-jul-12 ES

% revisions 
% 15-jul-12, xcsdhat and xlfphat returned as well
% 17-jul-12 enable load merged file
% 20-sep-12 check dat file existence
% 27-sep-12 'spc' option; also spectrogram of the CSD of a subset of selected channels
% 05-oct-12 added plot of actual waveforms on the image in case of SPC
% 09-oct-12 added 'wlt' version; also normalize spectra by random no-stim SWS spectra
% 15-oct-12 added refchan option - subtract this channel from all channels
%               prior to any computation
% 17-oct-12 added -1 option for refchan, this subtracts trial-averaged mean
%               from the data prior to spectral analysis. cannot use this
%               for LFP/CSD analysis because those look at the mean, so i
%               assume semi-stationarity and removed the mean of every
%               other trial in that case
% 02-nov-12 CSD := -diff( x'' )
% 05-nov-12 allow calling with par file instead of filebase
% 27-dec-12 modified for edges-csd case
% 09-jul-13 remove outliers before calculating images (two parameters -
%               number of SDs - default 5, faction of values - default 1%)
% 10-oct-13 (1) removed the resortidx (fixed readbin bug)
%           (2) added mother input to control wavelet base
%           (3) argument handling - by pairs
% 20-oct-13 (1) expand edges for CSD
% 10-nov-13 (1) apply temporal filtering before spectral analysis 
%                  (major reduction in spike bleed-in)
% 30-dec-13 enable whitening by a FIR (to be supplied in whitenFlag). then
%               will filter the raw data with this FIR and remove the
%               filtered version from the raw data. useful for artifact
%               removal; will be applied to output as well. 
% 01-jan-14 BUG in marginal case: when CSD is computed for an array of dual
%               channels (e.g. ch1 ch2 ch2), the wavelet was computed for (ch2-ch1)/2
%               instead of for ch2-ch1.
% 20-feb-15 recomputeBL: added option to force use of the existing file
% 10-sep-18 adapatation to Stark Lab:
%               (1) organized 
%               (2) handled case of no SWS
% 18-aug-19 cleaned up

% TO DO:
% (1) fix issue of timing of spectrograms (alignment) - solved for wavelet, and
% that is in general more accurate
% (2) externalize the MUA generation (modify dat2mua to accept ranges)

function [ avgcsd, avglfp, tim, xcsdhat, xlfphat, too, foo, yoo, yooCSD ] = pt_avg( filebase, channels, trigs, varargin )

%------------------------------------------------------------------%
% initialize
%------------------------------------------------------------------%
% output
avgcsd                      = [];
avglfp                      = [];
tim                         = [];
xcsdhat                     = [];
xlfphat                     = [];
foo                         = [];
too                         = [];
yoo                         = [];
yooCSD                      = [];

% required arguments
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    fprintf( 1, '%s: filebase missing\n', mfname )
    return
end
if nargs < 2 || isempty( channels )
    fprintf( 1, '%s: channels missing\n', mfname )
    return
end
if nargs < 3 || isempty( trigs )
    fprintf( 1, '%s: No triggers\n', mfname )
    return
end

% optional arguments
[ stawinSEC, tempSD, spatBin, suffix, graphics, scalefactor, specChans...
    , whitenFlag, normalizeSWS, savetype, refchan, recomputeBL, nSD, minFract, vflag...
    , expandEdges, baselineEpochs...
    , mother, nBins, Fmin, Fmax ] = ParseArgPairs(...
    { 'stawinSEC', 'tempSD', 'spatBin', 'suffix', 'graphics', 'scalefactor', 'specChans'...
    , 'whitenFlag', 'normalizeSWS', 'savetype', 'refchan', 'recomputeBL', 'nSD', 'minFract', 'vflag'...
    , 'expandEdges', 'baselineEpochs' ...
    , 'mother', 'nBins', 'Fmin', 'Fmax' }...
    , { [ -0.005 0.005 ], 0.0001, 1, 'eeg', 1, 1 / 2^16 * 8 * 1e3 / 1e3, []...
    , 0, 1000, 'png', [], 0, 3, 0, 1 ...
    , 1, [] ...
    , 'MORLET', 64, 2, 300 }...
    , varargin{ : } );

%------------------------------------------------------------------%
% determine source file
%------------------------------------------------------------------%
if isa( filebase, 'struct' ) && isfield( filebase, 'FileName' )
    par                     = filebase;
    filebase                = filebase.FileName;
    fname                   = sprintf( '%s.%s', filebase, suffix );
else
    fname                   = sprintf( '%s.%s', filebase, suffix );
    xmlfname                = sprintf( '%s.xml', filebase );
    par                     = LoadXml( xmlfname );
end
if length( stawinSEC ) ~= 2
    fprintf( 1, '%s: stawinSEC should be a 2-element vector\n', mfname )
    return
end
suffix                      = lower( suffix );
switch suffix
    case 'mua'
        dataType            = 'MUA';
    case { 'eeg', 'spc', 'wlt' }
        dataType            = 'LFP';
    case 'dat'
        dataType            = 'RAW';
    otherwise
        fprintf( 1, '%s: Unsupported suffix\n', mfname )
        return
end     
suffix0                     = suffix;
if ~exist( fname, 'file' ) && strcmp( suffix, 'mua' )
    MAKEMUA                 = 1;
    suffix                  = 'eeg';
    fname( end - 2 : end )  = 'eeg';
else
    MAKEMUA                 = 0;
end
if ~exist( fname, 'file' ) && ( strcmp( suffix, 'spc' ) || strcmp( suffix, 'wlt' ) )
    MAKESPC                 = 1;
    suffix                  = 'eeg';
    fname( end - 2 : end ) = 'eeg';
else
    MAKESPC                 = 0;
end
if ~exist( fname, 'file' )
    verb( sprintf( '%s: missing file %s!!!!!!!!!!!!', mfname, fname ), vflag )
    return
end
[ ~, filename, extname ]    = fileparts( fname );
switch suffix
    case { 'eeg', 'mua' }
        Fs                  = par.lfpSampleRate;
    case 'dat'
        Fs                  = par.SampleRate;
end

%------------------------------------------------------------------%
% determine parameters
%------------------------------------------------------------------%

% determine channels
if any( channels < 0 ) 
    shanknum                = abs( channels( channels < 0 ) );
    chans                   = [];
    for i                   = shanknum( : ).'
        chans               = [ chans par.ElecGp{ i } + 1 ];
    end
else
    chans                   = channels;
end
nchans                      = length( chans );
totchans                    = par.nChannels;
if isempty( specChans )
    specChans               = chans;
end
if ~isempty( specChans )
    if ~sum( ismember( channels, specChans ) )
        fprintf( 1, '%s: specChans not in channel list\n', mfname )
        return
    end
end
if nchans >= 3
    doCSD                   = 1;
else
    doCSD                   = 0;
    verb( sprintf( '%s: Cannot compute CSD (%d total channels)!!!', mfname, nchans ), vflag )
end
if expandEdges
    ncsdChans               = nchans;
else
    ncsdChans               = nchans - 2;
end
if refchan == -1
    rmvMEAN                 = 1;
    refchan                 = [];
else
    rmvMEAN = 0;
end
if ~isempty( refchan ) && refchan > 0 && refchan <= totchans && refchan == round( refchan )
    doREF                   = 1;
else
    doREF                   = 0;
end
    
% spatio-temporal filtering
win1                        = gausskernel( 0, ceil( tempSD * Fs ), 1, 6 * ceil( tempSD * Fs ) + 1 );
win1                        = win1 / sum( win1( : ) );
win2                        = triang( 1 + 2 * spatBin );
win2                        = win2 / sum( win2( : ) );
% NOTE:
% Presently there is some leakage between windows
% This does not influence the spectral analyses since it is taken care of 
% However, if the same approach is adopted for any temporal filtering, 
% should load periods padded to ceil( length( win1 ) / 2 ) at each side and then clip that 

% determine periods
stawin                      = [ ceil( stawinSEC( 1 ) * Fs ) floor( stawinSEC( 2 ) * Fs ) ]; % at Fs
ntrigs                      = length( trigs );
trigs                       = round( trigs * Fs );                                                  % at Fs
periods                     = trigs * ones( size( stawin ) ) + ones( size( trigs ) ) * stawin;
[ periods, pIdx ]           = unique( periods, 'rows' );
rmvPeriods                  = periods( :, 1 ) < 0;
periods( rmvPeriods, : )    = [];
pIdx( rmvPeriods, : )       = [];
nperiods                    = size( periods, 1 );
if length( unique(  diff( periods, 1,  2) ) ) > 1
    error( 'check type casting (all should be double)' )
end

%------------------------------------------------------------------%
% MUA computation
%------------------------------------------------------------------%
% load just the periods + edges, compute mua for those, discard edges
if MAKEMUA
    % MUA-specific parameters
    MUA_PAD                 = [ -0.02 0.02 ];
    MUA_TH                  = 2;
    MUA_HPF                 = [ 300 NaN ];
    MUA_LPF                 = [ NaN 100 ];
    verb( sprintf( '%s: Computing MUA from dat...', mfname ), -vflag )
    % prepare the periods
    fname_dat               = fname;
    fname_dat( end - 2 : end ) = 'dat';
    Fs_dat                  = par.SampleRate;
    DSF                     = round( Fs_dat / Fs );
    % pad the upsampled periods
    periods_dat             = ( periods + ones( nperiods, 1 ) * MUA_PAD * Fs ) * DSF ...
        + ones( nperiods, 1 ) * [ -DSF / 2 DSF / 2 - 1 ];
    % load the data
    if ~exist( fname_dat, 'file' )
        fprintf( 1, '%s: Missing file %s!!!\n', mfname, fname_dat )
        return
    end
    x                       = readbin( fname_dat, chans, totchans, periods_dat )' * scalefactor ;
    if doREF
        xref                = readbin( fname_dat, refchan, totchans, periods_dat ) * scalefactor;   % extract periods and scale
        x                   = bsxfun( @minus, x, mean( xref, 1 ) );
    end
    [ m, n ]                = size( x );
    % comptue the MUA for each channel separately (all seg. together)
    hSpikes                 = makefir( MUA_HPF, Fs_dat, [], 'high' ); 
    hAvg                    = makefir( MUA_LPF, Fs_dat, [], 'low' ); 
    xf                      = firfilt( x, hSpikes ); 
    % //clip extreme values//
    mmx                     = [ 1 1 ]' * mean( xf ) + [ 1 -1 ]' * MUA_TH * std( xf );
    for i = 1 : n
        xf( xf( :, i ) > mmx( 1, i ) ) = mmx( 1, i );
        xf( xf( :, i ) < mmx( 2, i ) ) = mmx( 2, i );
    end
    % //clip extreme values//
    xf                      = firfilt( xf.^2, hAvg ); 
    xf( xf < 0 )            = 0;
    y                       = xf( DSF / 2 : DSF : m, : );
    y                       = y .^ 0.5;
    % remove the pads
    keep                    = diff( periods( 1, : ), [], 2 ) + 1; % assume all periods are equal
    pads                    = MUA_PAD * Fs;
    vec                     = logical( [ zeros( abs( pads( 1 ) ), 1 ); ones( keep, 1 ); zeros( abs( pads( 2 ) ), 1 ) ] );
    kidx                    = repmat( vec, [ nperiods 1 ] );
    x                       = y( kidx, : )';
end

%------------------------------------------------------------------%
% spectral analysis
%------------------------------------------------------------------%
if MAKESPC
    verb( sprintf( '%s: Computing spectrograms from %s%s (%d triggers)...'...
        , mfname, filename, extname, nperiods ), -vflag )
    specMode                = suffix0;
    % parameters:
    switch specMode
        case 'spc'
            M               = 1;
            mtNW            = 3;
            dflag           = '';
            nBins           = 4;
            % derive two functional parameters: (1) how much to pad (2) temporal bin size
            nFFT            = 2^floor( log2( Fs * M ) ); % 128 for 1250
            nWindow         = floor( ( diff( stawin ) + 1 ) / nBins );
            nWindow         = nWindow + mod( nWindow, 2 );
            verb( sprintf( 'Multi-taper temporal/spectral bin sizes: %0.3gms/%0.3gHz...'...
                , nWindow / Fs * 1000, Fs / nFFT ), vflag )
        case 'wlt'
            % wavelet parameters: 2 300 64
            nWindow         = 0;
            verb( sprintf( 'Wavelet frequency range/bins: %0.3g-%0.3g Hz, %d bins...'...
                , Fmin, Fmax, nBins + 1 ), vflag )
    end
    
    % determine the analysis window
    stawinSECraw            = stawinSEC + [ -1 1 ] ...
        * ( nWindow / Fs + diff( stawinSEC ) / 2 ); % nWindow + half the analysis window on each side (plenty!)
    stawinraw               = [ floor( stawinSECraw( 1 ) * Fs ) ceil( stawinSECraw( 2 ) * Fs ) ]; % at Fs
    periods_eeg             = trigs * ones( size( stawinraw ) ) + ones( size( trigs ) ) * stawinraw;
    periods_eeg             = unique( periods_eeg, 'rows' );
    ridx                    = periods_eeg( :, 1 ) < 0;
    periods_eeg( ridx, : )  = [];
    snipdur                 = diff( stawinraw ) + 1;
    nsnips                  = size( periods_eeg, 1 );
    
    % load, scale, whiten
    verb( sprintf( '%s: Computing LFP spectra..', mfname ), -vflag )
    x0                      = readbin( [ filebase '.eeg' ], chans, totchans, periods_eeg ) * scalefactor;
    if doREF % remove a reference channel
        xref                = readbin( [ filebase '.eeg' ], refchan, totchans, periods_eeg ) * scalefactor;
        x0                  = bsxfun( @minus, x0, mean( xref, 1 ) );
    end
    if isequal( whitenFlag, 1 ) || isequal( whitenFlag, 2 )
        % whiten each channel seprately (but all segments together); discontinuities will be ignored later
        x                   = WhitenSignal( x0', [], 0 )';
    elseif length( whitenFlag ) > 1
        % remove a FIR-filtered version of the data
        x                   = x0 - firfilt( x0', whitenFlag )';
    else
        x                   = x0;
    end
    xmat                    = reshape( x, [ length( chans ) snipdur nsnips ] );                     % arragne in a 3D array
    xmat                    = squeeze( mean( xmat( ismember( chans, specChans ), :, : ), 1 ) );     % just the spectra of the selected channels
    if nperiods == 1
        xmat                = xmat';
    end
    xmat                    = bsxfun( @minus, xmat, mean( xmat, 1 ) );                              % set each trial's mean to zero
    if rmvMEAN                                                                                      % remove a trial-averaged mean (artifacts)
        xmat                = bsxfun( @minus, xmat, mean( xmat( :, 2 : 2 : end ), 2 ) );
    end
    if length( win1 ) > 1
        xmat                = firfilt( xmat, win1 );                                                % temporal smoothing
    end

    % generate a time-frequency representation
    yoo                     = [];                                                                   % freq, time, snippet
    for j                   = 1 : size( xmat, 2 )
        xj                  = double( xmat( :, j ) );
        switch specMode
            case 'spc'
                [ yoo( :, :, j ), foo, too ] = mtcsg( xj, nFFT, Fs, nWindow, nWindow/2, mtNW, dflag );
            case 'wlt'
                [ yoo( :, :, j ), foo ] = getWavelet( xj, Fs, Fmin, Fmax, nBins, mother );
                too         = ( stawinraw( 1 ) : stawinraw( 2 ) ) / Fs;
        end
    end
    % clip to the time requested
    if strcmp( suffix0, 'spc' )
        too                 = too + stawinSECraw( 1 ) + diff( too( 1 : 2 ) ); % sec
    end
    tidx                    = too >= stawinSEC( 1 ) & too <= stawinSEC( 2 );
    yoo                     = yoo( :, tidx, : );
    too                     = too( tidx );
    
    % normalize the raw data
    if normalizeSWS
        % get the baseline decomposition
        bfname              = sprintf( '%s.%sBL.%s', filebase, suffix0, num3str( specChans( 1 ) ) );
        if exist( bfname, 'file' ) && recomputeBL ~= 1
            L               = load( bfname, '-mat' );
            if ~isfield( L, 'mm' ) || ~isfield( L, 'ss' ) || ~isfield( L, 'snipdur' ) ...
                    || ~isfield( L, 'Fs' ) || recomputeBL == 0 && L.snipdur / L.Fs < diff( stawinSEC ) ...
                    || doCSD && ~isfield( L, 'mmCSD' ) || doCSD && ~isfield( L, 'ssCSD' ) ...
                    || doCSD && all( isnan( L.mmCSD ) ) || doCSD && all( isnan( L.ssCSD ) )
                recomputeBL = 1;
            else
                mm          = L.mm;
                ss          = L.ss;
                recomputeBL = 0;
            end
        end
        if ( ~exist( bfname, 'file' ) || length( specChans ) ~= 1 || recomputeBL )
            verb( sprintf( 'Determining normalization parameters..' ), -vflag )
            mm              = [];
            switch specMode
                case 'wlt'
                    [ mm, ss, freqs, S ] = computeBaseline( filebase, specChans, specMode...
                        , 'nrand', normalizeSWS, 'windur', diff( stawinSEC )...
                        , 'epochs', baselineEpochs...
                        , 'suffix', suffix, 'scalefactor', scalefactor, 'refchan', refchan...
                        , 'fMin', Fmin, 'fMax', Fmax, 'nfBins', nBins + 1, 'Overwrite', 0 );
                case 'spc'
                    [ mm, ss, freqs, S ] = computeBaseline( filebase, specChans, specMode...
                        , 'nrand', normalizeSWS, 'windur', diff( stawinSEC )...
                        , 'epochs', baselineEpochs...
                        , 'suffix', suffix, 'scalefactor', scalefactor, 'refchan', refchan...
                        , 'M', M, 'nWindow', nWindow, 'Overwrite', 0 );
            end
            if isempty( mm )
                mm          = zeros( length( foo ), 1 );
                ss          = ones( length( foo ), 1 );
                freqs       = foo( : )';
                S.snipdur   = ceil( diff( stawinSEC ) * Fs / 2 ) * 2 + diff( stawinSEC ) * Fs;
                S.nrand     = normalizeSWS;
            end
            if length( chans ) < 3 && length( specChans ) == 1
                snipdur0    = snipdur;
                snipdur     = S.snipdur;
                nrand       = S.nrand;
                fprintf( '%s: Saving %s..\n', upper( mfilename ), bfname )
                save( bfname, 'filebase', 'specMode', 'specChans', 'Fs', 'snipdur', 'nrand'...
                    , 'mm', 'ss',  'freqs', '-v6' );
                spindur     = snipdur0;
            end
        else
            verb( sprintf( 'Loading normalization factors ' ), -vflag )
            L               = load( bfname, '-mat' );
            verb( sprintf( '(%0.3g sec) ', L.snipdur / L.Fs ), -vflag )
            mm              = L.mm;
            ss              = L.ss;
        end
        
        % actually normalize
        verb( sprintf( 'Normalizing spectra..' ), -vflag )
        for k               = 1 : length( foo )
            yoo( k, :, : )  = ( yoo( k, :, : ) - mm( k ) ) / ss( k );
        end
        verb( 'Done', vflag )
            
    end % normalizeSWS
    
    % compute spectra of CSD of those segments
    if doCSD && ( expandEdges || sum( ismember( chans( 2 : end - 1 ), specChans ) ) > 0 )
        
        verb( sprintf( '%s: Computing CSD spectra..', mfname ), -vflag )
        spcCSD              = 1;

        % prepare the CSD segments
        x                   = bsxfun( @minus, x0, mean( x0, 2 ) );
        if expandEdges
            x               = -diff( x( [ 1 1 : size( x, 1 ) size( x, 1 ) ], : ), 2 );
        else
            x               = -diff( x, 2 ); % CSD
        end
        if isequal( whitenFlag, 1 ) || isequal( whitenFlag, 2 )
            x               = WhitenSignal( x', [], 0 )';
        elseif length( whitenFlag ) > 1
            x               = x - firfilt( x', whitenFlag )';
        end
        xmat                = reshape( x, [ ncsdChans snipdur nsnips ] );
        kidx                = ismember( chans, specChans );
        if length( specChans ) == 1 && sum( ismember( chans, specChans ) ) > 1
            if chans( 1 ) == specChans
                kidx( 1 )   = 0; 
            elseif chans( end ) == specChans
                kidx( end ) = 0; 
            end
        end
        xmat                = squeeze( mean( xmat( kidx, :, : ), 1 ) );
        if nperiods == 1
            xmat            = xmat';
        end
        xmat                = bsxfun( @minus, xmat, mean( xmat, 1 ) );
        if rmvMEAN
            xmat            = bsxfun( @minus, xmat, mean( xmat( :, 2 : 2 : end ), 2 ) );
        end
        if length( win1 ) > 1
            xmat            = firfilt( xmat, win1 );
        end
        
        % compute the CSD spectra
        yooCSD              = [];
        for j               = 1 : size( xmat, 2 )
            xj              = double( xmat( :, j ) );
            switch specMode
                case 'spc'
                    yooCSD( :, :, j ) = mtcsg( xj, nFFT, Fs, nWindow, nWindow/2, mtNW, dflag );
                case 'wlt'
                    yooCSD( :, :, j ) = getWavelet( xj, Fs, Fmin, Fmax, nBins, mother );
            end
        end
        yooCSD              = yooCSD( :, tidx, : );
        
        % normalize the CSD
        if normalizeSWS
            % get the baseline decomposition
            if exist( bfname, 'file' ) && ~recomputeBL
                L           = load( bfname, '-mat' );
                if ~isfield( L, 'mmCSD' ) || ~isfield( L, 'ssCSD' ) ...
                        || all( isnan( L.mmCSD ) ) || all( isnan( L.ssCSD ) )
                    recomputeBL = 1;
                end
            end
            if ~exist( bfname, 'file' ) || length( specChans ) ~= 1 || recomputeBL
                mmCSD       = [];
                switch specMode
                    case 'wlt'
                        [ mmCSD, ssCSD, freqs, S ] = computeBaseline( filebase, specChans, specMode...
                            , 'nrand', normalizeSWS, 'windur', diff( stawinSEC )...
                            , 'epochs', baselineEpochs...
                            , 'suffix', suffix, 'scalefactor', scalefactor, 'refchan', refchan...
                            , 'doCSD', 1, 'chans', chans...
                            , 'fMin', Fmin, 'fMax', Fmax, 'nfBins', nBins + 1, 'Overwrite', 0 );
                    case 'spc'
                        [ mmCSD, ssCSD, freqs, S ] = computeBaseline( filebase, specChans, specMode...
                            , 'nrand', normalizeSWS, 'windur', diff( stawinSEC )...
                            , 'epochs', baselineEpochs...
                            , 'suffix', suffix, 'scalefactor', scalefactor, 'refchan', refchan...
                            , 'doCSD', 1, 'chans', chans...
                            , 'M', M, 'nWindow', nWindow, 'Overwrite', 0 );
                end
            if isempty( mmCSD )
                mmCSD       = zeros( length( foo ), 1 );
                ssCSD       = ones( length( foo ), 1 );
                freqs       = foo( : )';
                S.snipdur   = ceil( diff( stawinSEC ) * Fs / 2 ) * 2 + diff( stawinSEC ) * Fs;
                S.nrand     = normalizeSWS;
            end
                if length( specChans ) == 1
                    snipdur0 = snipdur;
                    snipdur = S.snipdur;
                    nrand   = S.nrand;
                    fprintf( '%s: Saving %s..\n', upper( mfilename ), bfname )
                    save( bfname, 'filebase', 'specMode', 'specChans', 'Fs', 'snipdur', 'nrand'...
                        , 'mm', 'ss', 'freqs', 'mmCSD', 'ssCSD',  '-v6' );
                    snipdur = snipdur0;
                end
            else
                verb( sprintf( 'Loading CSD normalization factors ' ), -vflag )
                L           = load( bfname, '-mat' );
                verb( sprintf( '(%0.3g sec)..', L.snipdur / L.Fs ), -vflag )
                mmCSD       = L.mmCSD;
                ssCSD       = L.ssCSD;
            end
            % acutally normalize
            verb( sprintf( 'Normalizing CSD spectra..' ), -vflag )
            for k = 1 : length( foo )
                yooCSD( k, :, : ) = ( yooCSD( k, :, : ) - mmCSD( k ) ) / ssCSD( k );
            end
            verb( 'Done', vflag )
        end % normalize
        
    else
        
        spcCSD = 0;
        
    end % CSD
    
    verb( 'done spectra.', vflag )
end

if ~MAKEMUA
    verb( sprintf( '%s: loading %s data for %s%s...'...
        , mfname, suffix0, filename, extname ), -vflag )
    if isequal( whitenFlag, 2 ) || length( whitenFlag ) > 1
        x0                  = readbin( fname, chans, totchans, periods ) * scalefactor;             % load data en-block
        if doREF % remove a reference channel
            xref            = readbin( fname, refchan, totchans, periods_eeg ) * scalefactor;       % extract periods and scale
            x0              = bsxfun( @minus, x0, mean( xref, 1 ) );
        end
        if isequal( whitenFlag, 2 )
            % whiten each channel seprately (but all segments together); discontinuities will be ignored later
            x               = WhitenSignal( x0', [], 0 )';
        elseif length( whitenFlag ) > 1
            % remove a FIR-filtered version of the data
            x               = x0 - firfilt( x0', whitenFlag )';
        end

    else
        x                   = readbin( fname, chans, totchans, periods ) * scalefactor;             % load data en-block
    end
    
end

%------------------------------------------------------------------%
% CSD computation
%------------------------------------------------------------------%
if doREF && ~MAKEMUA && ~isequal( whitenFlag, 2 )                                                   %% think about what to do if whitenFlag is a FIR
    xref                    = readbin( fname, refchan, totchans, periods ) * scalefactor ;          % load data en-block
    x                       = bsxfun( @minus, x, mean( xref, 1 ) );
end

snipdur                     = diff( stawin ) + 1;
tim                         = ( stawin( 1 ) : stawin( 2 ) )' / Fs * 1000;
nsnips                      = size( periods, 1 );

x                           = bsxfun( @minus, x, mean( x, 2 ) ); 
if length( win1 ) > 1
    x                       = firfilt( x', win1 )';
end
xlfphat                     = reshape( x, [ nchans snipdur nsnips ] );
if rmvMEAN                  % the trick - remove half-data mean!
    verb( 'Removing LFP half-mean...', -vflag )
    for i                   = 1 : nchans
        xmat                = squeeze( xlfphat( i, :, : ) );  
        xmat                = bsxfun( @minus, xmat, mean( xmat( :, 2 : 2 : end ), 2 ) ); 
        xlfphat( i, :, : )  = xmat; 
    end
end
avglfp                      =  mean( xlfphat, 3 )';

% compute the smoothed CSD
if doCSD
    verb( 'Computing CSD...', -vflag )
    if expandEdges
        xcsd                = -diff( x( [ 1 1 : size( x, 1 ) size( x, 1 ) ], : ), 2 );
    else
        xcsd                = -diff( x, 2 );
    end
    if length( win2 ) > 1 && size( xcsd, 1 ) > length( win2 )
        xcsd                = firfilt( xcsd, win2 ); % spatial smoothing
    end
    xcsdhat                 = reshape( xcsd, [ ncsdChans snipdur nsnips ] );                        % convert to 3D arrays
    if rmvMEAN
        for i               = 1 : ncsdChans
            xmat            = squeeze( xcsdhat( i, :, : ) );
            xmat            = bsxfun( @minus, xmat, mean( xmat( :, 2 : 2 : end ), 2 ) );
            xcsdhat( i, :, : ) = xmat;
        end
    end
    avgcsd                  = nanmean( xcsdhat, 3 )';
end
verb( 'done computations!', vflag )

% remove outliers (replace w/ NaNs)
if doCSD
    % detect outliers based on the CSD of the specChans
    if expandEdges
        cidx                = ismember( chans, specChans );
    else
        cidx                = ismember( chans( 2 : end - 1 ), specChans );
    end
    rmvsource               = 'csd';
end
if ~doCSD || sum( cidx ) == 0
    cidx                    = ismember( chans, specChans );
    rmvsource               = 'lfp';
end
ridx                        = false( size( periods, 1 ), 1  );
for acidx                   = find( cidx( : ).' )
    switch rmvsource
        case 'csd'
            xx              = squeeze( xcsdhat( acidx, :, : ) );
        case 'lfp'
            xx              = squeeze( xlfphat( acidx, :, : ) );
    end
    oidx                    = outliers( xx, nSD, 0, 0, 2 );
    oidx2                   = outliers( sum( abs( xx ) ), nSD, 0, 0, 2 )';
    ridx                    = ridx | oidx2 | mean( isnan( oidx ) )' > minFract;                     % more than 1% of the values are outliers
end
% replace outliers w/ NaNs
if sum( ridx )
    xlfphat( :, :, ridx )   = NaN;
    xcsdhat( :, :, ridx )   = NaN;
    if MAKESPC
        yoo( :, :, ridx )   = NaN;
        yooCSD( :, :, ridx ) = NaN;
    end
    avglfp                  = nanmean( xlfphat, 3 )';
    if doCSD
        avgcsd              = nanmean( xcsdhat, 3 )';
    end
end

%------------------------------------------------------------------%
% equate input and output sizes
%------------------------------------------------------------------%
if length( pIdx ) ~= ntrigs
    padIdx                  = setdiff( 1 : ntrigs, pIdx );
    verb( sprintf( '%s: Padding output...', mfname ), -vflag )
    sx                      = size( xlfphat );
    tmp                     = NaN * ones( [ sx( 1 : 2 ) ntrigs ] );
    tmp( :, :, pIdx )       = xlfphat;
    xlfphat                 = tmp;
    if MAKESPC
        sx                  = size( yoo );
        tmp                 = NaN * ones( [ sx( 1 : 2 ) ntrigs ] );
        tmp( :, :, pIdx )   = yoo;
        yoo                 = tmp;
    end
    if doCSD
        fprintf( '%s: Expanding edges for CSD!!\n', upper( mfilename ) )
        sx                  = size( xcsdhat );
        tmp                 = NaN * ones( [ sx( 1 : 2 ) ntrigs ] );
        tmp( :, :, pIdx )   = xcsdhat;
        xcsdhat             = tmp;
        if MAKESPC && sum( ismember( chans( 2 : end - 1 ), specChans ) ) > 0
            sx              = size( yooCSD );
            tmp             = NaN * ones( [ sx( 1 : 2 ) ntrigs ] );
            tmp( :, :, pIdx ) = yooCSD;
            yooCSD          = tmp;
        end
    end
end

%------------------------------------------------------------------%
% plot
%------------------------------------------------------------------%
% non-spectral analyses:    [ LFP traces, CSD traces; LFP image+traces, CSD image+traces ]
% spectral analyses:        [ LFP spectrogram, CSD spectrogram; LFP image+traces, CSD image+traces ]

if graphics

    verb( sprintf( '%s: Plotting...', mfname ), -vflag )

    % parameters
    yCalib                  = 0.1;                                                                  % mV
    xCalib                  = diff( stawinSEC ) * 0.1 * 1000;                                       % 10% of the time range;
    calib                   = [ xCalib yCalib ];
    CALIBCOLOR              = [ 0 0.7 0 ];                                                          % dark green
    CALIBWIDTH              = 2;
    USF                     = 10;                                                                   % for images (spectral/non-spectral)
    SF                      = 5;                                                                    % only for non-spectral traces 
    if length( graphics ) ~= 4
        fig                 = figure;
        nsp                 = 4;
        for i               = 1 : nsp
            spi( i )        = subplot( nsp / 2, 2, i );
        end
    else
        spi                 = graphics;
    end
    xlims                   = tim( [ 1 end ] );

    % mean LFP
    subplot( spi( 3 ) ),
    [ x1, y1, z1 ]          = imupsample( tim, 1 : nchans, avglfp', USF );
    imagesc( x1, y1, z1 )
    axis xy
    ylim( 0.5/USF * [ -1 1 ] + [ 1 nchans ] )
    clims                   = round( get( gca, 'clim' ) * 1000 );                                   % mV->uV
    title( sprintf( '%s:[%d,%d] uV', dataType, clims( 1 ), clims( 2 ) ) )
    xlabel( 'Time (ms)' )
    ylabel( 'Channel' )
    ytick                   = 1 : nchans;
    set( gca, 'ytick', ytick, 'yticklabel', channels( ytick ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    xlim( xlims )
    axis square
    line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
    hold on
    [ ~, ph3 ]             = plotTraces( tim, avglfp, 1, [ 1 nchans 0 nchans + 1 ], calib );
    ylim( 0.5 / USF * [ -1 1 ] + [ 1 nchans ] )
    set( ph3, 'color', [ 1 1 1 ] * 0.3 );
    ch                      = colorbar; 
    set( ch, 'tickdir', 'out', 'box', 'off' ), 
    try
        subplot( ch ), 
        title( 'mV' )
    catch
    end
    
    % mean CSD
    if doCSD
        subplot( spi( 4 ) ),
        [ x1, y1, z1 ]      = imupsample( tim, 1 : ncsdChans, avgcsd', USF );
        imagesc( x1, y1, z1 )
        axis xy
        ylim( 0.5/USF * [ -1 1 ] + [ 1 nchans ] )
        clims               = round( get( gca, 'clim' ) * 1000 );
        title( sprintf( '%s:[%d,%d] uV', 'CSD', clims( 1 ), clims( 2 ) ) )
        xlabel( 'Time (ms)' )
        ylabel( 'Channel' )
        ytick               = 1 : nchans;
        set( gca, 'ytick', ytick, 'yticklabel', channels( ytick ) )
        set( gca, 'tickdir', 'out', 'box', 'off' )
        xlim( xlims )
        axis square
        line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
        hold on
        [ ~, ph4 ]         = plotTraces( tim, avgcsd, 2, [ 1 nchans 0 nchans + 1 ], calib );       % here spacing 2 yields a scaling of 2
        ylim( 0.5 / USF * [ -1 1 ] + [ 1 nchans ] )
        set( ph4, 'color', [ 1 1 1 ] * 0.3  );
        ch                  = colorbar; 
        set( ch, 'tickdir', 'out', 'box', 'off' )
    end
    
    % time-frequency (or only time) decomposition
    if strcmp( suffix0, 'spc' ) || strcmp( suffix0, 'wlt' )

        % parameters
        calibf              = 50;                                                                   % 50 Hz
        imageType           = 'contour';                                                            % 'contour' or 'image'; necessary if wavelet
        fMaxPlot            = Fs / 2;                                                               % nyquist
        if normalizeSWS
            ustr            = 'Z';
        else
            ustr            = 'mV2';
        end
        
        % LFP spectrum 
        subplot( spi( 1 ) )
        fidx                = foo >= 0 & foo <= fMaxPlot;
        tidx                = too >= stawinSEC( 1 ) & too <= stawinSEC( 2 );
        xx                  = too( tidx ) * 1000;
        yy                  = foo( fidx );
        if normalizeSWS
            zz              = nanmean( yoo( fidx, tidx, : ), 3 );
        else
            zz              = nanmean( 20*log10( abs( yoo( fidx, tidx, : ) + eps ) ), 3 );
        end
        plotSpectrogram( xx, yy, zz', USF, imageType, 'linear' );
        xlabel( 'Time (ms)' )
        ylabel( 'Freq (Hz)' )
        clims               = get( gca, 'clim' );
        title( sprintf( '%s%s, %d/%dT, C%s; [%0.2g,%0.2g]%s'...
            , replacetok( filename, '\_', '_' ), extname, nsnips - sum( ridx )...
            , nsnips, num2str( specChans( : ).' ), clims( 1 ), clims( 2 ), ustr ) )
        xlim( xlims )
        axis square
        line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );       
        calibration( [ calib( 1 ) calibf ], 1i* 0.05 * [ 1 1 ], gca, 0, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH );
        ch                  = colorbar;
        set( ch, 'tickdir', 'out', 'box', 'off' )
        
        % CSD spectrum 
        if spcCSD
            if normalizeSWS
                zz          = nanmean( yooCSD( fidx, tidx, : ), 3 );
            else
                zz          = nanmean( 20*log10( abs( yooCSD( fidx, tidx, : ) + eps ) ), 3 );
            end
            if nansum( zz( : ) )
                subplot( spi( 2 ) )
                plotSpectrogram( xx, yy, zz', USF, imageType, 'linear' );
                xlabel( 'Time (ms)' )
                ylabel( 'Freq (Hz)' )
                clims       = get( gca, 'clim' );
                title( sprintf( 'CSD: [%0.2g,%0.2g]%s'...
                    , clims( 1 ), clims( 2 ), ustr ) )
                xlim( xlims )
                axis square
                line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
                calibration( [ calib( 1 ) calibf ], 1i* 0.05 * [ 1 1 ], gca, 0, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH );
                ch          = colorbar;
                set( ch, 'tickdir', 'out', 'box', 'off' )
            end
        else
            axis( spi( 2 ), 'off' )
        end
                
    else
        
        % raw traces
        subplot( spi( 1 ) ),
        plotTraces( tim, avglfp, [], SF );
        ylims               = ylim;
        ylims               = [ ylims( 1 ) - 0.1 * abs( ylims( 1 ) ) ylims( 2 ) * 1.1 ];
        ylim( ylims )
        title( sprintf( '%s%s, %d snips', replacetok( filename, '\_', '_' ), extname, nsnips ) )
        axis square
        line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
        
        if doCSD
            subplot( spi( 2 ) ),
            csd_pad         = NaN * ones( size( tim ) );
            plotTraces( tim, [ csd_pad avgcsd csd_pad ], [], SF );
            ylim( ylims )
            axis square
            line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
        end
        
    end

    % finalize and save
    colormap( jet )
    figpath                 = graphics;
    if isa( figpath, 'char' ) && exist( figpath, 'dir' )
        figname             = sprintf( '%s/%s%s.%s.pta%d.%s', figpath, filename, extname, num3str( specChans, 3 ), nperiods, suffix0 );
        fig_out( fig, 1, [ figname '.' savetype ], savetype ),
        fprintf( 1, 'saving %s', figname )
    end
    verb( 'Done!', vflag );
    
end

return

% EOF

