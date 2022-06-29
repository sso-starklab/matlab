% thetaRatio        theta/non-theta ratio
%
% [ tdRatio theta nontheta ] = thetaRatio( eegfname, eegchan, nchans, thetaBP, nonthetaBP, Fs, filtMode, npoles, tdWindow, tdMode )
%
% -low level funciton intended to be called from within segmentBehavior, no argument checking, all arguments required
% -can be used to compute ratio of power in any two bands within the Nyquist
% -fir filtering is done with a linear phase fir with zero phase
% distortion; iir filtering is done with a Butterworth using filtfilt so
% the actual magnitude response is doubled
%
% INPUT:
% eegfname          full path + name + suffix; can be *.dat; or data
% eegchan           one channel
% nchans            number of channels in the source file
% thetaBP           [ hipass lowpass ], [Hz], e.g. [ 5 11 ]
% nonthetaBP        [ hipass lowpass ], [Hz], e.g. [ 2 4 ]
% Fs                [Hz]; of the source file, e.g. 1250 for 'eeg' files
% filtMode          'iir' or 'fir'
% npoles            relevant only for 'iir', use a low number (e.g. 2)
% tdWindow          [sec], use a number sufficiently large to accomodate
%                       several cycles of the lowest frequency, e.g. 2 sec
% tdMode            'rms' or 'hil'; use 'rms' (faster; identical results)
%
% OUTPUT:
% tdRatio           the ratio between the two bands in the running window
% theta             the bandpass filtered data
% nontheta          the bandpass filtered data
%
% calls             verb, makefir, makeblocks, firfilt, ma_rms
% 
% see also          segmentBehavior

% 07-apr-13 ES

% revisions
% 16-apr-19 edge condition in blocks handled
% 04-aug-19 edge condition in blocks handled

function [ tdRatio, theta, nontheta ] = thetaRatio( eegfname, eegchan, nchans, thetaBP, nonthetaBP, Fs, filtMode, npoles, tdWindow, tdMode )

verbose         = 1;
blocksize       = 2^20;

verb( sprintf( '%s: computing theta/delta ratio...', upper( mfilename ) ), -verbose )

% prepare filters/windows
switch lower( filtMode( 1 : 3 ) )
    case 'fir'
        hTheta      = makefir( thetaBP, Fs, [], 'bandpass' );
        hnonTheta   = makefir( nonthetaBP, Fs, [], 'bandpass' );
    case 'iir'
        [ btheta, atheta ]          = butter( npoles, thetaBP / Fs * 2, 'bandpass' );
        [ bnontheta, anontheta ]    = butter( npoles, nonthetaBP / Fs * 2, 'bandpass' );
end
win             = ones( tdWindow * Fs, 1 ) / Fs / tdWindow;

% determine block structure
if isa( eegfname, 'char' )
    a           = memmapfile( eegfname, 'Format', 'int16' );
    neeg        = length( a.Data ) / nchans;
elseif isa( eegfname, 'numeric' )
    eeg         = eegfname( : );
    neeg        = length( eeg );
end
if neeg         <= ( 10 * blocksize )
    blocks      = [ 1 neeg ];
    blocks0     = blocks;
    nblocks     = 1;
else
    noverlap                        = ceil( 3 * Fs * tdWindow );
    blocks0                         = makeblocks( neeg, blocksize, 0 );
    blocks                          = blocks0;
    nblocks                         = size( blocks, 1 );
    blocks( 2 : nblocks, 1 )        = blocks0( 2 : nblocks, 1 ) - noverlap;
    blocks( 1 : nblocks - 1, 2 )    = blocks0( 1 : nblocks - 1, 2 ) + noverlap;
    m                               = find( any( blocks > neeg, 2 ), 1 ); 
    if ~isempty( m )
        blocks( m, 2 )                  = blocks( nblocks, 2 ); 
        blocks( m + 1 : nblocks, : )    = []; 
        nblocks                         = m;
    end
end

verb( sprintf( 'Computing ratio for eeg channel %d; processing in %d block(s): '...
    , eegchan, nblocks ), -verbose )

% initialize output
tdRatio         = zeros( neeg, 1 );
if nargout      > 1
    theta       = zeros( neeg, 1 );
    nontheta    = zeros( neeg, 1 );
end

% compute
for bidx = 1 : nblocks
    
    verb( sprintf( '#%d ', bidx ), -verbose )
    
    % load eeg
    if exist( 'a', 'var' )
        idx     = ( ( blocks( bidx, 1 ) - 1 ) * nchans + eegchan ) : nchans : ( ( blocks( bidx, 2 ) - 1 ) * nchans + eegchan );
        eeg     = double( a.data( idx ) );
    end
    
    % filter the theta, non-theta
    if nblocks  == 1
        verb( 'filtering...', -verbose )
    end
    switch lower( filtMode( 1 : 3 ) )
        case 'fir'
            thetaBlock      = firfilt( eeg, hTheta );
            nonthetaBlock   = firfilt( eeg, hnonTheta );
        case 'iir'
            thetaBlock      = filtfilt( btheta, atheta, eeg );
            nonthetaBlock   = filtfilt( bnontheta, anontheta, eeg );
    end
    
    % compute the running ratio
    if nblocks == 1
        verb( 'computing ratio...', -verbose )
    end
    switch lower( tdMode( 1 : 3 ) ) % the ratios are identical
        case 'hil'
            % exactly twice the MS obtained with the same window
            thetaPower0     = firfilt( abs( hilbert( thetaBlock ) ).^2, win );
            nonthetaPower0  = firfilt( abs( hilbert( nonthetaBlock ) ).^2, win );
        case 'rms'
            thetaPower0     = ma_rms( thetaBlock, Fs * tdWindow, 0 );
            nonthetaPower0  = ma_rms( nonthetaBlock, Fs * tdWindow, 0 );
    end
    tdRatio0                = thetaPower0 ./ nonthetaPower0;
    tdRatioBlock            = firfilt( tdRatio0, win );
    
    % accumulate the results
    prePad                  = ( blocks0( bidx, 1 ) - blocks( bidx, 1 ) );
    ridx                    = [ 1 : prePad ( ( 1 + prePad + diff( blocks0( bidx, : ) ) ) : diff( blocks( bidx, : ) ) ) + 1 ];
    tdRatioBlock( ridx )    = [];
    tdRatio( blocks0( bidx, 1 ) : blocks0( bidx, 2 ) ) = tdRatioBlock;
    if nargout > 1
        thetaBlock( ridx )  = [];
        nonthetaBlock( ridx ) = [];
        theta( blocks0( bidx, 1 ) : blocks0( bidx, 2 ) ) = thetaBlock;
        nontheta( blocks0( bidx, 1 ) : blocks0( bidx, 2 ) ) = nonthetaBlock;
    end
    
end

% clean up
clear a
verb( 'Done!', verbose )

return

% EOF
