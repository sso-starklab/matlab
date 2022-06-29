% parse1channel                 parse stimulus events from continous data (single channel)
%
% CALL                          stim = parse1channel( filebase, channel )
%
% ARGUMENTS
% filebase          full path or par structure
% channel           any channel
%
% optional:
% suffix            {'eeg'} or 'dat'
% Overwrite:        1 to recompute and overwrite
%                   0 to compute (with writing but not overwriting)
%                   -1 to just load/compute if a file exists/doesn't (no writing)
%                   {-2} load if exists, compute and save if doesn't
% minAmpRelative    {0.001} of peak-to-peak/2 voltage range (i.e. for 20V
%                       peak-to-peak, this corresponds to 100 mV threshold)
% minDurationSEC    {0.0005}, i.e. 0.5 ms (note hard coded 2-sample min duration)
% minDutyCycle      {0.5}; relevant only for sines/chirps (sindetect.m argument)
% sdGaussSEC        {0}; filter with a Gaussian fir
%
% CALLS
%       ParseArgPairs, LoadXml
%       stim_make, get_stimchans, num3str, makegaussfir, makeblocks, calc_com
%       fft_upsample, resampleSig, firfilt, tempmatch, sindetect, parseSchmitt, mtcsg
%       lin_downsample, parse, resampleranges, stim_get
% 
% CALLED BY
%       parseNchannels
%
% NOTE
% the transformation from A2D units to physical values is done for the driving current

% 01-aug-19 ES based on parseOneChannel

% revisions
% 04-aug-19 (1) added support for transformation matrices in formatted file (see details in interp_mat_example.m)
%           (2) removed historical DEBUG, UI
% 17-aug-19 (1) annotated, cleaned up
% 11-dec-19 (1) WN is thresholded by minAmp( 1 )
% 05-oct-20 (1) changed mveeg from [2 2] to [0 0] to support small value of uLED

% future extensions
% (2) allow template matching for user-defined templates (not only the canonical WN)

function stim = parse1channel( filebase, channel, varargin )

% control parameters
vflag                       = 1;

%-----------------------------------------------------------------%
% constants
%-----------------------------------------------------------------%
% general processing
blocksize                   = 2^20;                 % [samples] block-wise loading of wide-band / eeg data
maxDur_DEFAULT              = 20;                   % [s], determines block overlap ONLY
spkFs                       = 20000;              	% [Hz], for time stamps of output structure

% adaptive median and threshold
medBuf                      = 10;                   % scalar, for determining minAmplitude multiples for THing median
BLweight                    = [ 2 1 ];              % hi will be the weighted mean of the BL and the minAmplitude
hilo                        = [ 1 1 / sqrt( 2 ) ];  % lo will be a fixed fraction of hi

% sine wave/zap/white noise detection
sineCC_eeg                  = [ 0.85 0.75 ];        % hi-frequency chirps (0-100) with eeg filter overshoot and are distorted
sineCC_dat                  = [ 0.85 0.85 ];        % reasonable fit for piecewise sine (sometimes clipped), monotonous zap
wnCC_DEFAULT                = 0.8;                  % template matching

% classification of idealized waveforms of invidividual events:
minDurationAbsolute             = 2;                % general minimum duration
mveeg                       = [ 0 0 ];              % eeg data: samples to remove from beginning/end of segment for pulses/psines...
mvdat                       = [ 0 0 ];              % dat data: "
% si := |(mean(x1)-mean(x2))/mean(x)|
% ai := mean(x)/max(x)
% ramp:             si = 0.5;   ai = 0.5;
% pulse:            si = 0;     ai = 1;
% triangle:         si = 0;     ai = 0.5;
% sine:             si = 0;     ai = 0.5;
% psine             si = 0;     ai between 0.5 and 1 (i.e. between sine and pulse)
% the differentiation between sine, psine, and triangle is done based on
% the cv of the 1st derivative: very low for a triangle, intermediate for a
% sine, high for a pulse (0.1, 0.4, and 2.1 for noiseless signals)
sErr                        = 0.25;                 % permitted error for symmetry index
aErr                        = 0.25;                 % permitted error amplitude index
pErr                        = 0.05;                 % deviation from max for psine plateau
nDownsample                 = -13;
wf_separatrix               = [ 0.6 0.25 ];

% frequency characterization of ZAP/SINES
udTH                        = [ 0.75 0.25 ];        % frequency determination of sine waves
mtNW                        = 3;                    % time-freq parameter
dflag                       = '';                   % no detrending - remove mean for each segment separately
twin                        = [ 0.25 0.75 ];        % robust regression of chirp slope

% supported stimulus types; the order here is critical
stimTypes                   = { 'WN', 'ZAP', 'SINES', 'MISC' };  % misc :=  psines, pulses, ramps, triangulars

%-----------------------------------------------------------------%
% initialization
%-----------------------------------------------------------------%
stim                        = stim_make;
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 2 || isempty( filebase ) || isempty( channel )
    return
end
chan                        = channel( 1 );
[ suffix, Overwrite...
    , minAmpRelative, minDurationSEC, minDutyCycle, sdGaussSEC, tbaseline ...
    , imat, wnCC, sineCC, maxDur ] = ParseArgPairs(...
    { 'suffix', 'Overwrite'...
    , 'minAmpRelative', 'minDurationSEC', 'minDutyCycle', 'sdGaussSEC', 'tbaseline' ...
    , 'imat', 'wnCC', 'sineCC', 'maxDur' }...
    , { 'eeg', -2 ...
    , 0.001, 0.0005, 0.5, 0, [] ...
    , [], wnCC_DEFAULT, [], maxDur_DEFAULT }...
    , varargin{ : } );

% *.prm.xml file
prmfile                     = [ filebase '.prm.xml' ];
if ~exist( prmfile, 'file' )
    fprintf( 1, '%s: Missing parameter file %s\n', mfname, prmfile )
    return
end
par                         = LoadXml( prmfile );

% files and paths
% eeg file
eegExists                   = 0;
eegfile                     = [ filebase '.eeg' ];
if exist( eegfile, 'file' )
    eegExists               = 1;
end
if ~eegExists
    eegfile                 = [ filebase '.lfp' ];
    if exist( eegfile, 'file' )
        eegExists           = 1;
    end
end
if ~eegExists
    fprintf( 1, '%s: Missing %s!\n', mfname, eegfile )
    return
end
% source file
sourcefile                  = [ filebase '.' suffix ];
if ~exist( sourcefile, 'file' )
    fprintf( 1, '%s: Missing %s!\n', mfname, sourcefile )
    return
end
[ ~, filename, extname ]    = fileparts( filebase );

% par-dependent parameters
switch suffix
    case 'eeg'
        Fs                  = par.lfpSampleRate;
        mvwin               = mveeg;
        if isempty( sineCC )
            sineCC          = sineCC_eeg;
        end
    case 'dat'
        Fs                  = par.SampleRate;
        mvwin               = mvdat;
        if isempty( sineCC )
            sineCC          = sineCC_dat;
        end
end
nchans                      = par.nChannels;
nFFT                        = 2 ^ floor( log2( Fs ) ); % resolve 1 Hz
fROI                        = [ 0 Fs / 2 ];

% chan specific parameters
[ chan, target, voltagerange, source ] = get_stimchans( par, chan );
if isempty( chan ) || isempty( source )
    [ chan, target, voltagerange, source ] = get_stimchans( par, chan, 'trig' );
end
if isempty( chan ) || isempty( source )
    return
end
source                      = source{ 1 };
if target > 0
    shank                   = target;
else
    shank                   = [];
end
gain                        = 2 .^ par.nBits / voltagerange;                 % a2du / gain -> V

% basic detection parameters
if minAmpRelative < 0
    minAmplitude            = abs( minAmpRelative );
else
    minAmplitude            = abs( voltagerange ) / 2 * minAmpRelative;    % 0.1% of the maxAmp
end
minDuration                 = minDurationSEC * Fs;                          % [samples]
minDuration                 = ceil( max( minDuration, minDurationAbsolute ) );
if tbaseline > 0
    nwin                    = round( tbaseline * Fs );
    tbwin                   = ones( nwin, 1 ) / nwin;                       % MA
else
    tbwin                   = [];
end

% check whether to run or just load
savename        = sprintf( '%s.stm.%s', filebase, num3str( channel ) );
if Overwrite < 0 && exist( savename, 'file' )
    fprintf( 1, '%s: loading %s...\n', mfname, savename )
    load( savename, '-mat' )
    if stim_check( stim )
        [ nums, types ]     = uhist( stim.types );
        str = '';
        for i               = 1 : length( nums )
            str             = sprintf( '%s%d %s; ', str, nums( i ), types{ i } );
        end
        if i > 0
            fprintf( 1, '%ss.\n', str( 1 : end - 2 ) )
        else
            fprintf( 1, 'No stimuli detected.\n' )
        end
        return
    else
        fprintf( 1, 'format mismatch!\n' )
        Overwrite           = abs( Overwrite );
        stim                = stim_make;
    end
end

% determine smoothing window
if sdGaussSEC == 0
    g                       = [];
    f                       = 1;
else
    g                       = makegaussfir( sdGaussSEC, Fs );               % fir
    f                       = 1 / ( sqrt( sdGaussSEC * Fs ) * 1.96 );       % variance reduction factor
end

%------------------------------------------------------------------------%
% determine median by 1st pass through the data
%------------------------------------------------------------------------%
fprintf( 1, '%s: parsing stimulus data for %s channel %d, shank %d (%s): \n'...
    , mfname, [ filename extname ], channel, shank, source )

% compute median
fprintf( 1, '\t\tComputing median... ' ) 
e                           = memmapfile( sourcefile, 'Format', 'int16' );
n                           = length( e.Data ) / nchans;
d                           = zeros( 2 ^ par.nBits, 1 );                % allocate memory for the histogram of the data
edges                       = ( ( -2 ^( par.nBits - 1 ) - 0.5 ) : ( 2 ^( par.nBits - 1 ) - 0.5 ) )';
b                           = edges( 1 : end - 1 ) + 0.5;               % bin centers
blocks                      = makeblocks( n, blocksize, 0 );
nblocks                     = size( blocks, 1 );
fprintf( 1, '%d blocks: ', nblocks )
for bidx                    = 1 : nblocks
    fprintf( '%d ', bidx )
    idx                     = ( blocks( bidx, 1 ) * nchans - nchans + chan ) : nchans : ( blocks( bidx, 2 ) * nchans - nchans + chan );
    xb                      = single( e.data( idx ) );
    h                       = histc( xb, edges );
    h( end )                = [];
    d                       = d + h;
end
clear e
idx                         = b < ( medBuf * minAmplitude * gain );     % compute noise based only on low-amplitude (medBuf * minAmplitude) data
bnoise                      = b( idx );
dnoise                      = d( idx );
cnoise                      = cumsum( dnoise );
if cnoise( end ) == 0
    medx                    = 0;
else
    medx                    = bnoise( find( cnoise >= ( cnoise( end ) / 2 ), 1, 'first' ) );
end

% determine minAmp
if minAmpRelative < 0
    % force using some external value
    minAmp( : )             = abs( minAmpRelative ) * hilo';
else
    % adapt detection threshold
    b                       = ( b - medx ) / gain;           	% shift bin centers to physical measures (Iout)
    idx                     = b < minAmplitude;
    [ mNoise, sNoise ]      = calc_com( b( idx ), d( idx ) );   % COM of histogram yields mean and SD of raw data
    BL                      = mNoise + 3 * f * sNoise;          % includes correction for variance reduction if Gaussian filter will be applied to data
    minAmp                  = calc_com( [ BL minAmplitude ], BLweight ) * hilo;         % minAmp is a dual-TH for Schimdt triggering
end

% convert threshold if required (make sure not two equal, no zeros)
if ~isempty( imat )
    minAmp1                 = imat( find( imat( :, 1 ) >= minAmp( 1 ), 1, 'first' ) + 1, 2 );
    if isempty( minAmp1 )
        minAmp1             = interp1( imat( :, 1 ), imat( :, 2 ), minAmp( 1 ), 'linear', 'extrap' );
    end
    minAmp2                 = minAmp1 / minAmp( 1 ) * minAmp( 2 );
    minAmp                  = [ minAmp1 minAmp2 ];
end
fprintf( 1, ': %0.3g A2DU\n\t\tminAmp: %0.3g/%0.3g V\n'...
    , medx, minAmp( 1), minAmp( 2 ) )

%------------------------------------------------------------------------%
% parse the relevant stimulation channel (2nd pass through the same data)
%------------------------------------------------------------------------%
a                           = memmapfile( sourcefile, 'Format', 'int16' );
n                           = length( a.Data ) / nchans;
noverlap                    = 2 ^ ceil( log2( maxDur * Fs ) );
blocks                      = makeblocks( n, blocksize, noverlap );
nblocks                     = size( blocks, 1 );

% start filling output structure
stim.filebase               = filebase;
stim.suffix                 = suffix;
stim.chan                   = chan;
stim.voltagerange           = voltagerange;         % note: this refers to range/DAQ gain
stim.source                 = source;               % LED/LD/DPSS..
stim.duration               = n / Fs;               % [s]
stim.median                 = medx;                 % [a2du]
stim.generator              = { computer, datestr( now, 'ddmmmyy' ) };

fprintf( 1, '\t\tDetecting events: %d blocks...', nblocks )
ntypes                      = length( stimTypes );
for i                       = 1 : ntypes
    
    stimType                = upper( stimTypes{ i } );
    
    % detect WN by template matching
    WN                      = [];
    if isempty( WN )
        WN                  = load( 'WN10k' );      % frozen 1 s WN @ 24414.1625
        WN.template         = WN.WN;                % mean 0, SD 1
    end
    USF                     = Fs / WN.Fs;
    if USF == 1
        template            = WN.template;
    else
        if USF              == round( USF )
            template        = fft_upsample( WN.template, USF, 1 );
        else
            template        = resampleSig( WN.template, WN.Fs, Fs );
            if USF == 0.125
                USF     = 1;                    % flag; USF will be used for fine-tuning the detection
            end
        end
    end
    
    % computation
    switch stimType
            
        case 'WN'
            mat1                = [];       % sines
            mat2                = [];       % zap
            mat3                = [];       % misc
            onset               = [];       % wn
            wntau               = [];       % for/rev
            
            for bidx = 1 : nblocks
                fprintf( 1, '%d ', bidx )
                if vflag && bidx == nblocks, fprintf( 1, '\n' ); end
                idx             = ( blocks( bidx, 1 ) * nchans - nchans + chan ) : nchans : ( blocks( bidx, 2 ) * nchans - nchans + chan );
                xb              = ( single( a.data( idx ) ) - medx ) / gain;                        % load the data segment and transform A2DU to current [mA]
                if ~isempty( imat )
                    xb          = interp1( imat( :, 1 ), imat( :, 2 ), xb, 'linear', 'extrap' );    % apply reverse non-linear transform to P [mW]
                end
                if ~isempty( tbwin )
                    xbase       = firfilt( xb, tbwin );     % HPF: correct for slow baseline drift
                    xb          = xb - xbase;
                end
                if ~isempty( g )
                    xb          = firfilt( xb, g );         % LPF: correct for very fast noise
                end
                % detect WN
                c               = tempmatch( template, xb );
                ob              = round( mean( parse( find( c > wnCC ) ), 2 ) );
                % threshold WN according to minAmp
                tlen            = length( template );
                midx            = ob * ones( 1, tlen ) + ones( size( ob, 1 ), 1 ) * ( ( 1 : tlen ) - 1 );
                mx              = mean( xb( midx ), 2 );
                ridx            = mx < minAmp( 1 );
                if ~(length (ob)< length(ridx) )
                ob( ridx, : )   = [];
                else
                    ob = [];
                end
                % detect sines/zap:
                [ pidx, aidx, trains, meandurs, si, zi, clu ] = sindetect( xb, minAmp, sineCC, minDutyCycle, minDuration );
                tb1             = trains( clu == 1, : ); % sines
                mat1            = [ mat1; blocks( bidx, 1 ) + tb1 - 1 ];
                tb2             = trains( clu == 2, : ); % zap
                mat2            = [ mat2; blocks( bidx, 1 ) + tb2 - 1 ];
                % detect misc
                tb              = parseSchmitt( xb, minAmp( 1 ), minAmp( 2 ) );
                mat3            = [ mat3; blocks( bidx, 1 ) + tb - 1 ];
                % fine-tune WN detection ( ----not clear why required----- )
                if ~isempty( ob ) && USF ~= 1
                    cand        = [ ob ob + length( template ) - 1 ];
                    ob          = tb( isoverlap( tb, cand ), 1 );
                end
                onset           = [ onset; blocks( bidx, 1 ) + ob - 1 ];
                wntau           = [ wntau; ones( length( ob ), 1 ) ];
            end
            % summarize WN
            onset               = unique( onset );
            trainidx            = [ onset onset + length( template ) - 1 ];
            trainidx            = uniteranges( trainidx );
            % mid-summarize MISC
            mat3( ( diff( mat3, [], 2 ) + 1 ) < minDuration, : ) = [];
            if isempty( mat3 ) && isempty( trainidx )
                fprintf( '\t\tNo threshold-crossing events !!!\n' )
                break 
            end
            % dilute non-WN
            trainidx1           = uniteranges( mat1 );
            trainidx2           = uniteranges( mat2 );
            trainidx3           = uniteranges( mat3 );
        case 'ZAP'
            trainidx2( isoverlap( trainidx2, stim.times ), : ) = [];
            trainidx            = trainidx2;
        case 'SINES'
            trainidx1( isoverlap( trainidx1, stim.times ), : ) = [];
            trainidx            = trainidx1;
        case 'MISC'
            trainidx3( isoverlap( trainidx3, stim.times ), : ) = [];
            trainidx            = trainidx3;
    end % stimType
    
    % determine the parameters for each stimulus type
    nt          = size( trainidx, 1 );
    verb( sprintf( '\t\t%s: \t\t%d events detected. ', stimType, nt ), -vflag )
    if nt       == 0
        verb( '', vflag );
        continue
    end
    verb( sprintf( 'Computing parameters... ' ), -vflag )
    durs                = diff( trainidx, [], 2 ) + 1; % samples
    types               = cell( nt, 1 );
    slopes              = zeros( nt, 2 );
    maxvals             = zeros( nt, 1 );
    means               = zeros( nt, 1 );
    SDs                 = zeros( nt, 1 );
    franges             = zeros( nt, 2 );
    plateaus            = zeros( nt, 1 );
    for j               = 1 : nt
        period          = trainidx( j, : ) + [ -1 1 ];
        period          = [ max( period( 1 ), 1 ) min( period( 2 ), n ) ]; % enable derivative
        idx             = ( period( 1 ) * nchans - nchans + chan ) : nchans : ( period( 2 ) * nchans - nchans + chan );
        xj              = ( double( a.data( idx ) ) - medx ) / gain;
        slopes( j, : )  = [ diff( xj( [ 1 2 ] ) ) diff( xj( end - 1 : end ) ) ];
        if ( diff( period ) - 1 ) == durs( j )
            xj( [ 1 end ] ) = [];
        end
        % pulses,psines @ eeg - take central part only (due to filtering-induced distortion)
        if ~ismember( stimType, { 'ZAP', 'SINES', 'WN' } ) && length( xj ) > sum( mvwin )
            sidx        = [ 1  length( xj ) ] + mvwin .* [ 1 -1 ];
            xj          = xj( sidx( 1 ) : sidx( 2 ) );
        end
        means( j )      = mean( xj );
        SDs( j )        = std( xj );
        maxvals( j )    = max( xj );
        
        switch stimType
            case 'ZAP'
                if vflag && ~mod( j, 10 )
                    verb( '.', -vflag )
                end
                xj                  = xj - mean( xj );
                if floor( length( xj ) / nFFT ) > 1
                    nfft            = nFFT;
                else % shorten nFFT to have at least two sampling points
                    nfft            = 2 ^( floor( log2( length( xj ) ) ) - 1 );
                end
                [ yoo, foo, too ]   = mtcsg( xj, nfft, Fs, nfft, nfft / 2, mtNW, dflag );
                fidx                = foo >= fROI( 1 ) & foo <= fROI( 2 );
                com                 = calc_com( foo( fidx ), yoo( fidx, : ) )';
                tidx                = round( length( too ) * twin( 1 ) ) : round( length( too ) * twin( 2 ) );
                b                   = polyfit( too( tidx ), com( tidx ), 1 );
                b                   = [ b( 2 ) b( 1 ) ];
                f                   = b( 1 ) + b( 2 ) * too( [ 1 end ] )';
                [ minval, minidx ]  = min( f );
                [ maxval, maxidx ]  = max( f );
                f( minidx )         = max( floor( minval ), 0 );
                f( maxidx )         = ceil( maxval );
                franges( j, : )     = f;
                types{ j }          = stimType;
            case 'SINES'
                if vflag && ~mod( j, 10 )
                    verb( '.', -vflag )
                end
                types{ j }          = stimType;
                cros                = SchmittTrigger( xj, maxvals( j ) * udTH( 1 ), maxvals( j ) * udTH( 2 ) );
                f                   = round( Fs / mean( diff( cros ) ) );
                franges( j, : )     = f;
            case 'WN'
                if vflag && ~mod( j, 10 )
                    verb( '.', -vflag )
                end
                types{ j }          = stimType;
                franges( j, : )     = wntau( j );
            otherwise
                % classify by splitting the signal into two and using the descision rules noted above
                if vflag && ~mod( j, 10 )
                    verb( '.', -vflag )
                end
                nj                  = length( xj );
                if nj               <= 3
                    types{ j }      = 'PULSE';
                    continue
                end
                if ( nj / 2 )       ~= round( nj / 2 )
                    nj              = nj - 1;
                end
                xj                  = firfilt( xj, ones( 3, 1 ) / 3 );
                si                  = abs( ( mean( xj( 1 : nj / 2 ) ) - mean( xj( nj / 2 + 1 : nj ) ) ) / mean( xj ) );
                ai                  = ( mean( xj ) - minAmp( 1 ) ) / ( max( xj ) - minAmp( 1 ) );
                cv                  = NaN;
                
                if si               > ( 0.5 - sErr ) && ai >= ( 0.5 - aErr ) && ai <= ( 0.5 + aErr )
                    type = 'RAMP';                                                  % si = 0.5, ai = 0.5
                elseif si           < ( 0 + sErr )
                    if ai           >= ( 1 - aErr )
                        type        = 'PULSE';                                             % si = 0, ai = 1
                    else
                        newLen      = abs( nDownsample );
                        oldLen      = length( xj );
                        if oldLen   > ( 2 * newLen + 1 )
                            % if too long, downsample to the cannonical reference
                            x1      = diff( lin_downsample( xj, nDownsample ), 1 );
                        else % if too short, expand
                            xj0     = interp1( 1 : oldLen, xj, 1 : ( ( oldLen - 1 ) / ( newLen - 1 ) ) : oldLen );
                            x1      = reshape( xj0, [], 1 );
                        end
                        cv          = std( abs( x1 ) ) / mean( abs( x1 ) );
                        if cv       > wf_separatrix( 1 )
                            type    = 'PSINE';                                         % cv = 2.1
                        elseif cv   > wf_separatrix( 2 )
                            type    = 'SINE';                                          % cv = 0.4
                        else
                            type    = 'TRIANG';                                        % cv = 0.1
                        end
                    end
                else
                    type            = 'UNKNOWN';
                end
                types{ j }          = type;
                if strcmp( types{ j }, 'PSINE' )
                    vec             = parse( find( xj >= maxvals( j ) * ( 1 - pErr ) ) );
                    if ~isempty( vec )
                        plateaus( j ) = ( diff( vec( 1, : ) ) + 1 ) / Fs;
                    end
                end
        end % switch stimType
        
    end % for j = 1 : nt
    
    % accumulate
    stim.types                  = [ stim.types; types ];               %
    stim.times                  = [ stim.times; trainidx ];            % [ Fs samples ]
    stim.slopes                 = [ stim.slopes; slopes ];             % [ dV/dS ]
    stim.vals                   = [ stim.vals; maxvals ];              % [ V ]
    stim.stats                  = [ stim.stats; [ means SDs ] ];       % [ V, V ]
    stim.franges                = [ stim.franges; franges ];           % [ Hz Hz ]
    stim.plateaus               = [ stim.plateaus; plateaus ];         % [ s ]
    
    verb( sprintf( 'Done. ' ), vflag )
    
end % for i = 1 : ntypes

stim.types                          = upper( stim.types );
stim.times                          = resampleranges( stim.times, spkFs, Fs );   % [ spkFs samples ]
stim.durs                           = ( diff( stim.times, [], 2 ) + 1 ) / spkFs;  % [ s ], precise
clear a

%------------------------------------------------------------------------%
% finalize - sort and save
%------------------------------------------------------------------------%
% sort the events according to onset time
if ~isempty( stim.times )
    [ ~, sidx ]                     = sort( stim.times( :, 1 ) );
    stim                            = stim_get( stim, sidx );
end
% organize data for evt/val
mat                                 = [ stim.times stim.slopes stim.vals ];
% make fields for multi-channel information
nevents                             = size( mat, 1 );
if nevents
    stim.index                      = false( nevents, 1 );
    stim.category                   = zeros( nevents, 3 );
end

% save stim file
if channel                          < 100
    suf                             = sprintf( 't%s', num3str( channel, 2 ) );
else
    suf                             = num2str( channel );
end
evtfname                            = sprintf( '%s.evt.%s', filebase, suf );
valfname                            = sprintf( '%s.val.%s', filebase, suf );
if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( savename, 'file' ) )
    verb( sprintf( '\t\tsaving %s...', savename ), vflag )
    save( savename, 'stim', '-v6' )
    sstm                            = 1;
else
    sstm                            = 0;
end

% save evt file (ndm format - see Hazan et al., 2006, JNM)
if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( evtfname, 'file' ) )
    if nevents > 0
        evtime                      = mat( :, 1 : 2 ) / spkFs * 1000;               % convert to ms
        evlabel                     = ones( size( evtime, 1 ), 1 ) * [ 1 2 ];       % assign labels: 1-onset; 2-offset
        fid                         = fopen( evtfname, 'wt' );
        fprintf( fid, '%10.2f  %2.0f\n', [ evtime( : )'; evlabel( : )' ] );         % full resolution <=100h of recording...
        fclose( fid );
        verb( sprintf( '\t\tWritten evt file for channel %d (%d events)', channel, size( mat, 1 ) ), vflag )
    else        
        if exist( evtfname, 'file' )
            verb( sprintf( '\t\tRemoving existing evt file for channel %d!!!', channel ), vflag )
            if ispc
                system( sprintf( 'del /f %s', evtfname ) );
            else
                system( sprintf( 'rm -f %s', evtfname ) );
            end
        end
    end
elseif sstm
    verb( sprintf( '%s: NOTE potential mismatch between evt and stm files!', mfname ), vflag )
end

% save val file 
if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( valfname, 'file' ) )
    if nevents > 0
        fid                         = fopen( valfname, 'wt' );
        fprintf( fid, '%10.0f %10.0f %10.5f %10.5f %10.5f\n', mat' );                       % full resolution for up to 18 bits
        fclose( fid );
        verb( sprintf( '\t\tWritten val file for channel %d (%d events)', channel, size( mat, 1 ) ), vflag )
    else        
        if exist( valfname, 'file' )
            verb( sprintf( '\t\tRemoving val file for channel %d!!!', channel ), vflag )
            if ispc
                system( sprintf( 'del /f %s', valfname ) );
            else
                system( sprintf( 'rm -f %s', valfname ) );
            end
        end
    end
elseif sstm
    verb( sprintf( '%s: NOTE potential mismatch between val and stm files!', mfname ), vflag )
end

return

% EOF

