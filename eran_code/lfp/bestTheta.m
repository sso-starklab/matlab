% bestTheta             determine channel with the "best" theta signal
%
% [ eegchan, phs ] = bestTheta( filebase, par, neurochans, periods, thetaBP, nonthetaBP, Overwrite )
%
% filebase      full path or date/num cell array
% par           structure from LoadXml
% neurochans    channels to consider, defaults to par.SpkGrps
% periods       {[]}; periods to consider, eeg Fs (see readbin.m)
% thetaBP       {[5 11]}, Hz
% nonthetaBP    {[2 4]}, Hz
% Overwrite      1: compute and overwrite
%                0: compute, do not overwrite (but do write if does not exist)
%               -1: load/compute if a file exists/doesn't (no writing at all)
%               -2: load if existing, compute and save if not (i.e. make sure exists)
% graphics      {0}; graphical summary fig
% specMode      {'welch'}; can use 'mt' or 'timedomain' with linearly more time (not necessary)
%
% does:         determine the channel with the "best" theta signal from all recorded
%               channels. an optional argument periods allows to focus on specific
%               epochs, e.g. high speed/head accleration, otherwise all data are
%               considered; an optional argument chans focuses on a subset of channels
%
% calls:        LoadXml, get_egroup, verb, makeblocks
%               my_spectrum, mtcsd, outliers
%
% files:
% requires:     filebase.eeg
% optional:     filebase.xml
% output:       filebase.phs
%
% eegchan       channel with highest theta/non-theta ratio
% phs           theta phase (by 2-pole Butter)

% 18-nov-12 ES

% algorithm:
% (1) according to the par file, go over all neuro channels, compute the welch
% spectrum (fastest; alternatively could compute the multi-taper or wavelet
% spectrum, but not needed for this implementation)
% (2) determine for each channel the theta/non-theta ratio, and
% define "best" by the theta-nontheta ratio
% (3) for this channel, compute the theta phase using band-pass filtering
% (or wavelet, but that's again the same)
% (4) save the channel number and it's theta phase in a filebase.phs file
% (see phAnalysis for the format). Use bandpass filtering with an iir and
% hilbert

% revisions
% 05-dec-12 partitioned into blocks (just for whole file case, not for periods)
% 08-jan-13 indexing bug corrected
% 03-apr-13 saved pxx0 as well (the raw spectra of all neuro-channels)
%           full support for periods
% 07-apr-13 (1) block-wise phase computation
%           (2) apply spatial smoothing
% 30-aug-18 (1) modified legend call to allow compatibility with R2018a
%           (2) modified neurochans dimension to allow compatibility with R2018a
% 17-aug-19 cleaned up

% to do:
% (1) compute the phase block-wise /DONE/
% (2) check whether a theta-band oscillation actually exists (e.g. by MP)

function [ eegchan, phs, pxx0, f ] = bestTheta( filebase, par, neurochans, periods, thetaBP, nonthetaBP, Overwrite, graphics, specMode )

%------------------------------------------------------------------------
% constants
verbose                     = 1;
npoles                      = 2;             % relevant for the phase estimation/timedomain
NW                          = 3;                 % relevant only for mt

% to speed things up
BLOCKSIZE                   = 2^20;       % bytes/block; use 1 MB blocks (of the file)
nbytes                      = 2;             % for an int16

% defaults
THETA                       = [ 5 11 ];       % theta
NONTHETA                    = [ 2 4 ];     % delta

%------------------------------------------------------------------------
% initialize output
eegchan                     = [];
phs                         = [];

%------------------------------------------------------------------------
% arguments
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
eegfname                    = sprintf( '%s.eeg', filebase );
[ ~, filename ]             = fileparts( eegfname );
if ~exist( eegfname, 'file' )
    fprintf( 1, '%s: cannot proceed without an eeg file\n', mfname )
    return
end
if nargs < 2 || isempty( par )
    xmlfname                = sprintf( '%s.xml', filebase );
    if ~exist( xmlfname, 'file' )
        fprintf( 1, '%s: cannot proceed without an xml file\n', mfname )
        return
    end
    par                     = LoadXml( xmlfname);
end
if nargs < 3 || isempty( neurochans )
    neurochans              = [];
end
if nargs < 4 || isempty( periods )
    periods                 = [];
end
if nargs < 5 || isempty( thetaBP ) || length( thetaBP( : ) ) ~= 2
    thetaBP                 = THETA;
end
if nargs < 6 || isempty( nonthetaBP ) || length( nonthetaBP( : ) ) ~= 2
    nonthetaBP              = NONTHETA;
end
if nargs < 7 || isempty( Overwrite )
    Overwrite               = 0;
end
if nargs < 8 || isempty( graphics )
    graphics                = 0;
end
if nargs < 9 || isempty( specMode )
    specMode                = 'welch';
end
if NW == 0
    specMode = 'welch'; 
end
specMode                    = lower( specMode );
if ~ismember( specMode, { 'welch', 'mt', 'timedomain' } )
    fprintf( 1, '%s: unsupported mode\n', mfname )
    return
end

% check if existing
phsfname                    = sprintf( '%s.phs', filebase );
if Overwrite < 0 && exist( phsfname, 'file' )
    verb( sprintf( '%s: Loading theta phase data, %s...', mfname, phsfname ), -verbose )
    load( phsfname, '-mat' )
    if nargout > 2 && ~exist( 'pxx0', 'var' )
        verb( sprintf( 'Missing field pxx0; recomputing...' ), verbose );
        Overwrite           = 1;
    else
        verb( sprintf( 'C%d; Done!', eegchan ), verbose );
        return
    end
end

% parameters
nchans                      = par.nChannels;
Fs                          = par.lfpSampleRate;
nFFT                        = 2^floor( log2( Fs ) );
if isempty( neurochans )
    [ ~, neurochans ]       = get_egroup( par );
end
neurochans                  = intersect( neurochans, 1 : nchans );
[ btheta, atheta ]          = butter( npoles, thetaBP / Fs * 2, 'bandpass' );
if strcmp( specMode, 'timedomain' )
    [ bnontheta, anontheta ] = butter( npoles, nonthetaBP / Fs * 2, 'bandpass' );
end

% load the data (one channel at a time) and compute the spectrum
verb( sprintf( '%s: Determining spectrum for %s: ', mfname, filename ), -verbose )
pxx0                        = zeros( ceil( ( nFFT + 1 ) / 2 ), length( neurochans ) );
power0                      = zeros( 1, length( neurochans ) );
power1                      = power0;
a                           = memmapfile( eegfname, 'Format', 'int16' );

% prepare blocks and indices
totsamples                  = length( a.data ) / nchans;
sampblk                     = floor( BLOCKSIZE / nbytes );
if isempty( periods )
    blocks                  = makeblocks( totsamples, sampblk, 0 );
else
    blocks                  = periods;
    blocks( blocks < nFFT, : ) = [];
end
nsampsblk                   = diff( blocks, [], 2 ) + 1;
nsamples                    = sum( nsampsblk );
nblocks                     = size( blocks, 1 );
switch specMode
    case { 'welch', 'mt' }
        pxxBlock            = zeros( ceil( ( nFFT + 1 ) / 2 ), nblocks );
    case 'timedomain'
        power1Block         = zeros( 1, nblocks );
        power0Block         = zeros( 1, nblocks );
end
if sum( nsampsblk > sampblk ) > 0
    fprintf( '%s: long blocks; should sub-partition!!!', upper( mfilename ) );
end

% actually compute the spectra
for i                       = 1 : length( neurochans )
    chan                    = neurochans( i );
    verb( sprintf( '%d ', chan ), -verbose );
    
    % block-wise:
    t0                      = clock;
    for bidx                = 1 : nblocks
        idx                 = ( blocks( bidx, 1 ) * nchans - nchans + chan ) : nchans : ( blocks( bidx, 2 ) * nchans - nchans + chan );
        eeg                 = single( a.data( idx ) );
        eeg                 = eeg - mean( eeg );
        switch specMode
            case 'welch'
                [ pxxBlock( :, bidx ), f ]  = my_spectrum( eeg, nFFT, Fs, nFFT, nFFT / 2, 'none' );
            case 'mt' % ~nTapers longer (linear)
                [ pxxBlock( :, bidx ), f ]  = mtcsd( eeg, nFFT, Fs,  nFFT, nFFT / 2, NW, '' );
            case 'timedomain' % ~3 times longer
                theta                       = filtfilt( btheta, atheta, double( eeg ) );
                nontheta                    = filtfilt( bnontheta, anontheta, double( eeg ) );
                power1Block( bidx )         = sum( abs( hilbert( theta ) ) );
                power0Block( bidx )         = sum( abs( hilbert( nontheta ) ) );
        end
    end
    
    % summarize
    switch specMode
        case { 'welch', 'mt' }
            pxx0( :, i )    = sum( pxxBlock .* ( ones( size( pxxBlock, 1 ), 1 ) * nsampsblk' ), 2 ) / sum( nsampsblk );
        case 'timedomain'
            power1( i )     = sum( power1Block );
            power0( i )     = sum( power0Block );
    end
    et( i )                 = etime( clock, t0 );
    
end
verb( '', verbose );

% determine the channel with the best theta/nontheta ratio
verb( sprintf( '%s: Best theta/non-theta ratio: ', mfname ), -verbose )
if ~strcmp( specMode, 'timedomain' )
    fidx1                   = f >= thetaBP( 1 ) & f <= thetaBP( 2 );
    fidx0                   = f >= nonthetaBP( 1 ) & f <= nonthetaBP( 2 );
    power1                  = sum( pxx0( fidx1, : ), 1 );
    power0                  = sum( pxx0( fidx0, : ), 1 );
end
tdratio                     = power1 ./ power0;
rmv                         = outliers( power0, 3 ) | outliers( power1, 10 ); % remove outliers
win                         = triang( 5 ); 
win                         = win / sum( win );
tdratio( ~rmv )             = firfilt( tdratio( ~rmv ), win ); % smooth spatially
[ sRatio, sidx ]            = sort( tdratio( ~rmv ) );
validchans                  = neurochans( ~rmv );
eegchanidx                  = sidx( end );
eegchan                     = validchans( eegchanidx );
if strcmp( specMode, 'timedomain' )
    pxx                     = [];
    pxx0                    = [];
else
    pxx                     = pxx0( :, eegchanidx );
end
verb( sprintf( 'C%d; global ratio: %0.3g', eegchan, sRatio( end ) ), verbose )

% compute the phase for that channel
verb( sprintf( '%s: Computing theta phase data...', mfname ), -verbose )
if nsamples <= ( 10 * BLOCKSIZE )
    eeg                     = double( a.data( eegchan : nchans : end ) );
    theta                   = filtfilt( btheta, atheta, eeg );
    phs                     = single( angle( hilbert( theta ) ) );
else
    noverlap                = ceil( Fs / min( thetaBP ) * pi );
    blocks2                 = blocks;
    blocks2( 2 : nblocks, 1 )       = blocks( 2 : nblocks, 1 ) - noverlap;
    blocks2( 1 : nblocks - 1, 2 )   = blocks( 1 : nblocks - 1, 2 ) + noverlap;
    chan                    = eegchan;
    phs                     = zeros( nsamples, 1, 'single' );
    for bidx                = 1 : nblocks
        idx                 = ( blocks2( bidx, 1 ) * nchans - nchans + chan ) : nchans : ( blocks2( bidx, 2 ) * nchans - nchans + chan );
        eeg                 = double( a.data( idx ) );
        thetaBlock          = filtfilt( btheta, atheta, eeg );
        phsBlock            = angle( hilbert( thetaBlock ) );
        prePad              = ( blocks( bidx, 1 ) - blocks2( bidx, 1 ) );
        ridx = [ 1 : prePad ( ( 1 + prePad + diff( blocks( bidx, : ) ) ) : diff( blocks2( bidx, : ) ) ) + 1 ];
        phsBlock( ridx )    = [];
        phs( blocks( bidx, 1 ) : blocks( bidx, 2 ) ) = phsBlock;
    end
    theta = thetaBlock;
end
verb( sprintf( 'Done!' ), verbose )

% clean up
if isempty( periods )
    clear a
end

% graphical summary
if graphics
    
    figure,
    
    subplot( 2, 2, 1 )
    [ h0, w ]               = freqz( btheta, atheta, 2^nextpow2( Fs ), Fs );
    if strcmp( specMode, 'timedomain' )
        [ h1, w ]           = freqz( bnontheta, anontheta, 2^nextpow2( Fs ), Fs );
        ph                  = plot( w, abs( h0 ), 'b', w, abs( h1 ), 'r' );
    else
        ph                  = plot( w, abs( h0 ), 'b' );
    end
    set( ph, 'linewidth', 2 )
    xlim([ 0 max( [ thetaBP( : ); nonthetaBP( : ) ] ) * 2 ] )
    grid on
    title( 'Frequency response of the phase filter' );
    
    subplot( 2, 2, 2 )
    bh                      = bar( neurochans, tdratio );
    set( bh, 'edgecolor', [ 0 0 1 ], 'facecolor', [ 0 0 1 ] );
    hold on
    bh                      = bar( eegchan, tdratio( neurochans == eegchan ) );
    set( bh, 'edgecolor', [ 1 0 0 ], 'facecolor', [ 1 0 0 ] );
    neurochansRange         = neurochans( [ 1 end ] );
    xlim( neurochansRange( : )' + [ -1 1 ] )
    xlabel( 'Channel' ), ylabel( 'T/D ratio' )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    title( replacetok( filename, '\_', '_' ) )
    
    subplot( 2, 2, 3 )
    if strcmp( specMode, 'timedomain' )
        bh                  = bar( [ mean( nonthetaBP ) mean( thetaBP ) ], [ power0( neurochans == eegchan ) power1( neurochans == eegchan ) ] );
        set( bh, 'edgecolor', [ 0 0.7 0 ], 'facecolor', [ 0 0.7 0 ] );
        xlim( xlim + [ -1 1 ] )
        ylabel( 'Power' )
    else
        line( f, log10( pxx + eps ), 'linewidth', 2 )
        xlim( [ 0 100  ] )
        ylabel( 'Power (log)' )
    end
    xlabel( 'Frequency (Hz)' )
    title( sprintf( 'Best channel: %d', eegchan ) )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    
    subplot( 2, 2, 4 )
    tidx                    = 1 : min( length( eeg ), 1 * Fs );
    plot( tidx / Fs, [ eeg( tidx ) theta( tidx ) phs( tidx ) * std( theta ) ] );
    xlim( tidx( [ 1 end ] ) / Fs )
    legend( 'wideband', 'theta', 'phase' );%, 0 );
    xlabel( 'Time (sec)' )
    set( gca, 'box', 'off', 'tickdir', 'out' )
    
    
end

% save phs file
if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( phsfname, 'file' ) )
    verb( sprintf( '%s: saving %s...', mfname, phsfname ), verbose )
    save( phsfname, 'phs', 'eegchan', 'pxx', 'f', 'pxx0', 'nsamples', '-v6' );
end

return

% EOF
