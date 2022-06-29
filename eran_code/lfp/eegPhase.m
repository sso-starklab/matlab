% eegPhase          estimate the 1-frequency phase of 1 channel from a binary file
%
% call              [ phs, t ] = eegPhase( infile, chan, periods, fir )
%                   [ ..., uphs, cycs, talign ] = ( ..., outfile, Fin, Fout, nchans, method )
%                       ( ..., expand, trigs, graphics )
%
% gets              infile              filebase (taken as 'eeg') or full name w/ path
%                   channel             any channel
%                   periods             [s]; default: the entire file
%                   fir                 determines filtering. either 
%                                           band ([Hz]; 2-element vector),
%                                           vector (a fir), or
%                                           cell array of vectors (sequentially applied firs)
%                                       defaults by method: extrema:    [ 1 60 14 ] (bandpass 1-60, then peak detection 14 Hz)
%                                                           others:     [ 50 NaN ] (highpass 50 Hz)
%
% optional arguments:
%                   outfile             optional, then output will be saved to disk
%                   Fin                 [Hz]; if not supplied, taken from the xml file
%                                           according to infile suffix
%                   Fout                [Hz]; if not supplied, same as Fin
%                   nchans              if not supplied, taken from the xml file
%                   method              {'hilbert'}; string passed to calcPhase ('hilbert', 'wavelet')
%                                           or to calcPhaseExtrema ('extrema')
%                   expand              [s]; default: 0.001 s
%                   trigs               trigger points
%                   graphics            {0}; interval (blocks) to plot
%
% returns           phs                 phase (0:2pi; 0-peak)
%                   t                   sample number (at Fout)
%                   uphs                unwrapped phase (unwrapPhase)
%                   cycs                cycle number (monotonic)
%                   talign              first peak below cycle0 trough (unwrapPhase)
%
% GENERAL
%   this routine is optimized for both periods (multiple ripples, for
%   instance) and huge data sets (e.g. theta from a long merged file). It
%   supports any source/output sampling frequency, e.g. eeg or dat files,
%   with output at any equal or higher frequency than the input. Thus, it
%   does not downsample; to do that, simply:
%   >> mod( downsample( unwrap( phs ), DSF ), 2 * pi )
% 
% DOES
%   wrapper for calcPhase/calcPhaseExtrema, which compute a phase from 
%   a vector/matrix with filtering but without any disk/memory handling 
%   or edge padding
%
% DOES NOT
%   compute phases at multiple frequency bands. Thus to obtain several phases
%   this  routine needs to be called multiple times (and the data has to be 
%   loadded multiple times), potentially wasteful if the objective is to 
%   obtain phases of multiple different bands; use eegPhaseBank for this.
% 
% SAVED FORMAT
%   the format implemented saves in a binary format only the phases
%   themselves. this corresponds to the assumption that the entire file is
%   analyzed. 
%   thus, when 'periods' are used, the binary saved does not correspond
%   to the entire file. however, indices (t) can be returned. then, it is 
%   the responsibility of the calling routine to handle indices.
%
% calls             LoadXml                         (blab)
%                   makeblocks, monotonic, verb     (general)
%                   resampleranges                  (sets)
%                   calcPhase, calcPhaseExtrema, 
%                           makefir, monotonic_phase, 
%                           unwrapPhase             (ssp) 
%
% see also          calcPhase
%                   calcPhaseExtrema

% 10-feb-13 ES

% revisions
% 11-feb-13 upsampling errors corrected
% 26-jul-13 (1) integrated the upsampling properly
%           (2) integrated unwrapping the phase (according to an external
%           trigger)
% 28-jul-13 moved phase unwrapping to an external routine
% 10-dec-14 modified to support phase computation of a long file
% 21-mar-21 cleaned up
% 10-apr-21 (1) bug fixes
%           (2) implemented extrema method
%           (3) graphics argument added
% 12-apr-21 (1) bug fix (block spillover)
%           (2) error checking for file size
% 29-jun-21 (1) added smoothing resolution of waveform phase
%           (2) added enforcement of monotonic phase
% 01-jul-21 (1) added default "fir" for extrema

function [ phs, t, uphs, cycs, talign ] = eegPhase( infile, chan, periods, fir, outfile, Fin, Fout, nchans, method, expand, trigs, graphics )

% constants
vflag                           = 1;
blocksize                       = 2 ^ 16;                                   % [samples]
iirDP                           = 2 * pi / 2;                               % [rad]

% output arguments
phs                             = [];
t                               = [];
uphs                            = [];
cycs                            = [];
talign                          = [];

% arguments
nargs                           = nargin;
nout                            = nargout;
if nargs < 2
    return
end
if nargs < 3 || isempty( periods )
    periods                     = [];
end
if nargs < 4
    fir                         = [];
end
if nargs < 5 || isempty( outfile )
    outfile                     = '';
end
if nargs < 6 || isempty( Fin )
    Fin                         = [];
end
if nargs < 7 || isempty( Fout )
    Fout                        = [];
end
if nargs < 8 || isempty( nchans )
    nchans                      = [];
end
if nargs < 9 || isempty( method )
    method                      = 'hilbert';
end
method                          = lower( method );
if ~ismember( method, { 'hilbert', 'peaks', 'troughs', 'wavelet', 'extrema' } )
    error( 'unsupported method!' )
end
if isempty( fir )
    switch method
        case 'extrema'
            fir                 = [ 1 60 14 ];
        otherwise
            fir                 = [ 50 NaN ];
    end
end
if nargs < 10 || isempty( expand )
    expand                      = 0.001;
end
if nargs < 11 || isempty( trigs ) 
    trigs                       = [];
end
if nargs < 12 || isempty( graphics ) 
    graphics                  	= 0;
end

mfname                          = mfilename;

%------------------------------------------------------------------%
% determine parameters
%------------------------------------------------------------------%
pathname                        = fileparts( infile );
if ~exist( infile, 'file' )
    infile                      = [ infile '.eeg' ];
end
if ~exist( infile, 'file' )
    verb( sprintf( '%s: missing input file %s', mfname, infile ), vflag )
    return
end
if isempty( outfile ) && nout == 0
    verb( sprintf( '%s: no output arguments and no output file', mfname ), vflag )
    return
end
if isempty( Fin ) || isempty( nchans )
    xml                         = dir( [ pathname '/*.xml' ] );
    if isempty( xml )
        verb( sprintf( '%s: cannot determine Fin / nchans', mfname ), vflag )
        return
    end
    xmlfname                    = [ pathname '/' xml( 1 ).name ];
    par                         = LoadXml( xmlfname );
end
if isempty( Fin )
    suffix                      = infile( end - 2 : end );
    switch suffix
        case 'eeg'
            Fin                 = par.lfpSampleRate;
        case 'dat'
            Fin                 = par.SampleRate;
        otherwise
            verb( sprintf( '%s: cannot determine Fin', mfname ), vflag )
            return
    end
end
if isempty( nchans )
    nchans                      = par.nChannels;
end
if isempty( Fout )
    Fout                        = Fin;
end
expand                          = ceil( expand * Fout );
method                          = lower( method );
usf                             = ceil( Fout / Fin );

%------------------------------------------------------------------%
% determine filter
%------------------------------------------------------------------%
mindur                          = [];
if isequal( method, 'extrema' )
    if isa( fir, 'cell' )
        verb( sprintf( '%s: for method ''extrema'', an IIR filter is used (''fir'' should be a 3 element vector)', mfname ), vflag )
        return
    end
    if length( fir ) ~= 3
        verb( sprintf( '%s: for method ''extrema'', an IIR filter is used (''fir'' should be a 3 element vector)', mfname ), vflag )
        return
    end
    iirBP                       = fir( 1 : 2 );
    iirDT                       = Fout / fir( 3 );
else
    if isa( fir, 'cell' )
        fdurs                   = zeros( length( fir ), 1 );
        for i                   = 1 : length( fir )
            fdurs( i )          = numel( fir{ i } );
        end
        mindur                  = max( fdurs );
    elseif isvector( fir )
        if strcmp( method, 'wavelet' )
            if length( fir ) <= 2
                mindur          = Fin / min( fir );
            end
            fir                 = [ fir( 1 : 2 ) Fout ];
        else
            fir                 = fir( ~isnan( fir ) );
            if length( fir ) == 1
                fir             = makefir( [ fir NaN ], Fout, [], 'hipass' );
            elseif length( fir ) == 2
                fir             = makefir( fir, Fout, [], 'bandpass' );
            end
            mindur              = length( fir );
        end
    end
    if isempty( mindur )
        verb( sprintf( '%s: input format mismatch', mfname ), vflag )
        return
    end
    iirBP                       = [];
    iirDT                       = [];
end
mindur                          = 3 * mindur;

%------------------------------------------------------------------%
% determine block structure
%------------------------------------------------------------------%
a                               = memmapfile( infile, 'Format', 'int16' );
n                               = length( a.Data ) / nchans;
if n ~= round( n )
    verb( sprintf( '%s: file %s contains %d int16, which is not an integer multiple of %d' ...
        , mfname, infile, length( a.Data ), nchans ) ...
        , vflag )
    return
end

if isempty( periods )
    % divide into blocks for manageable memory handling
    ewin                        = 1;                                        % 1 s - a worst case scenario.. should be 1/min(freq)
    noverlap                    = ceil( Fin * ewin * pi );
    blocks0                     = makeblocks( n, blocksize, 0 );
    blocks                      = blocks0;
    nblocks                     = size( blocks, 1 );
    blocks( 2 : nblocks, 1 )    = blocks0( 2 : nblocks, 1 ) - noverlap;
    blocks( 1 : nblocks - 1, 2 ) = blocks0( 1 : nblocks - 1, 2 ) + noverlap;
    blocks( blocks > n )        = n;
    nsamples                    = n;
    mode                        = 'int';
else
    % externally-specified periods:
    periods                     = resampleranges( periods * Fin, Fout, Fin );
    nperiods                    = size( periods, 1 );
    durs                        = diff( periods, [], 2 ) + 1;               % [samples]
    fidx                        = durs < mindur;
    pad                         = zeros( nperiods, 1 );
    if sum( fidx )
        pad( fidx, : )          = ceil( ( mindur - durs( fidx ) ) / 2 );
    end
    if sum( ~fidx )
        pad( ~fidx, : )         = expand;
    end
    periods                     = periods + pad * [ -1 1 ];
    cidx                        = pad * [ 1 1 ] + [ ones( nperiods, 1 ) durs ];
    blocks                      = resampleranges( periods, 1, usf );
    blocks( blocks < 1 )        = 1; 
    blocks( blocks > n )        = n;
    nblocks                     = size( blocks, 1 );
    nsamples                    = sum( diff( periods, 1, 2 ) + 1 );
    mode                        = 'ext';
end

%------------------------------------------------------------------%
% analyze data (blockwise) and save to output
%------------------------------------------------------------------%
if graphics
    fig1                        = figure;
    fig2                        = figure;
    set( fig1, 'position', [ 20   715   840   630 ] )
    set( fig2, 'position', [ 861   715   840   630 ] )
end

if ~isempty( outfile )
    fprintf( '%s: extracting %s phases from %d blocks (from %s->%s)\n'...
        , upper( mfilename ), method, size( blocks, 1 ), infile, outfile )
end
if ~isempty( outfile )
    fp                          = fopen( outfile, 'w' );
end
blocksUps                       = resampleranges( blocks, Fout, Fin );      % this is a bit more data than the periods
if isequal( mode, 'int' )
    if Fin ~= Fout
        error( 'presently not supported' )
    end
end
if isempty( trigs )
    trigsUps                    = blocksUps( :, 1 );
else
    trigsUps                    = round( trigs * Fout );
end

% allocate memory
talign                          = NaN * ones( nblocks, 1 );
if nout > 0 && isequal( mode, 'int' )
    phs                         = zeros( nsamples, 1, 'single' );
    if nout > 1
        t                       = phs;
    end
    if nout > 2
        uphs                    = phs;
    end
    if nout > 3
        cycs                    = phs;
    end
end

% go over blocks
for bidx                        = 1 : nblocks
    
    % grab
    idx                         = round( ( blocks( bidx, 1 ) * nchans - nchans + chan ) ) : nchans : round( ( blocks( bidx, 2 ) * nchans - nchans + chan ) );
    xb                          = double( a.Data( idx ) );
    
    % compute
    if isequal( method, 'extrema' )
        y                       = calcPhaseExtrema( xb, iirBP, Fout, iirDT, iirDP, 0 );
    else
        y                       = calcPhase( xb, fir, method, usf, 0 );
    end
    
    % enforce monotonic phase
    y                           = monotonic_phase( y );
    
    % clip
    if isequal( mode, 'int' )
        prePad                  = ( blocks0( bidx, 1 ) - blocks( bidx, 1 ) );
        ridx                    = [ 1 : prePad ( ( 1 + prePad + diff( blocks0( bidx, : ) ) ) : diff( blocks( bidx, : ) ) ) + 1 ];
        y( ridx )               = [];
        z                       = single( y );
    else
        kidx                    = ( cidx( bidx, 1 ) : cidx( bidx, 2 ) ) +  periods( bidx, 1 ) - blocksUps( bidx, 1 );
        z                       = single( y( kidx ) );
    end
    
    if graphics && ~mod( bidx, round( graphics ) )
        if local_plotter( xb, blocks, bidx, Fout, nchans, chan, a, iirBP, iirDT, iirDP, fir, method, usf, fig1, fig2 )
            return
        end
    end
    
    % write out to file
    if ~isempty( outfile )
        fwrite( fp, z, 'float32' );
    end
    
    % accumulate results
    if nout > 0
        if isequal( mode, 'int' )
            idx                 = blocks0( bidx, 1 ) : blocks0( bidx, 2 );
            phs( idx )          = z;
        else
            phs                 = [ phs; z ];
        end
    end
    if nout > 1 
        tidx                    = ( blocksUps( bidx, 1 ) : blocksUps( bidx, 2 ) )';
        if isequal( mode, 'int' )
            tidx( ridx )        = [];
            t( idx )            = tidx;
        else
            tidx                = tidx( kidx );
            t                   = [ t; tidx ];
        end
    end
    if nout > 2
        if isempty( trigs )
            gtrf                = 1;
        else
            gtrf                = find( tidx == trigsUps( bidx ) );         % the global trough by external determination (e.g. eeg trough closest to peak power)
        end
        [ uz, talign( bidx ) ]  = unwrapPhase( z, gtrf );
        if isequal( mode, 'int' )
            if bidx == 1
                euphs           = 0;
            else
                euphs           = uphs( blocks0( bidx - 1, 2 ) );
            end
            uphs( idx )         = ceil( euphs / ( 2 * pi ) ) * 2 * pi + uz;
        else
            uphs                = [ uphs; uz ];
        end
    end
    if nout > 3
        cycnum                  = floor( uz / ( 2 * pi ) );
        nnans                   = ~isnan( cycnum );
        cycnumt                 = [ -flipud( monotonic( flipud( -cycnum( cycnum < 0 ) ) ) ); cycnum( cycnum == 0 ); monotonic( cycnum( cycnum > 0 ) ) ];
        cycnum( nnans )         = cycnumt;
        if isequal( mode, 'int' )
            cycs( idx )         = cycnum;
        else
            cycs                = [ cycs; cycnum ];
        end
    end
    
end

clear a
if ~isempty( outfile )
    fclose(fp);
end

return


%------------------------------------------------------------------%
% local_plotter - internal routine
%------------------------------------------------------------------%

function rc = local_plotter( xb, blocks, bidx, Fout, nchans, chan, a, iirBP, iirDT, iirDP, fir, method, usf, fig1, fig2 )

xidxB                           = length( xb ) - 1 * Fout + ( 1 : Fout );   % last second
xidxA                           = xidxB + blocks( bidx, 1 ) - 1;

idxAll                          = round( ( blocks( 1, 1 ) * nchans - nchans + chan ) ) : nchans : round( ( blocks( bidx, 2 ) * nchans - nchans + chan ) );
xbAll                           = double( a.Data( idxAll ) );

figure( fig1 )
subplot( 2, 1, 1 )
if isequal( method, 'extrema' )
    [ phsB, xfB ]               = calcPhaseExtrema( xb, iirBP, Fout, iirDT, iirDP, 1 );
    xlim( xidxB( [ 1 end ] ) / Fout )
else
    [ phsB, xfB ]               = calcPhase( xb, fir, method, usf, 1 );
    xlim( xidxB( [ 1 end ] ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end
title( sprintf( 'Block#%d: %d to %d; samples %d to %d', bidx ...
    , blocks( bidx, 1 ), blocks( bidx, 2 ) ...
    , blocks( bidx, 1 ) + xidxB( 1 ) - 1 ...
    , blocks( bidx, 1 ) + xidxB( end ) - 1 ) )

subplot( 2, 1, 2 )
if isequal( method, 'extrema' )
    [ phsA, xfA ]               = calcPhaseExtrema( xbAll, iirBP, Fout, iirDT, iirDP, 1 );
    xlim( xidxA( [ 1 end ] ) / Fout )
else
    [ phsA, xfA ]               = calcPhase( xbAll, fir, method, usf, 1 );
    xlim( xidxA( [ 1 end ] ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
end

title( sprintf( 'Entire dataset: %d to %d', blocks( 1, 1 ), blocks( bidx, 2 ) ) )

figure( fig2 )
subplot( 3, 1, 1 )
plot( xidxB, xb( xidxB ), 'b', xidxB, xbAll( xidxA ), '--r' )
xlim( xidxB( [ 1 end ] ) )
set( gca, 'tickdir', 'out', 'box', 'off' )
title( sprintf( 'Block#%d: %d to %d; samples %d to %d', bidx ...
    , blocks( bidx, 1 ), blocks( bidx, 2 ) ...
    , blocks( bidx, 1 ) + xidxB( 1 ) - 1 ...
    , blocks( bidx, 1 ) + xidxB( end ) - 1 ) )

subplot( 3, 1, 2 )
plot( xidxB, xfB( xidxB ), 'b', xidxB, xfA( xidxA ), '--r' )
xlim( xidxB( [ 1 end ] ) )
set( gca, 'tickdir', 'out', 'box', 'off' )
title( 'filtered' )

subplot( 3, 1, 3 )
ph                              = plot( xidxB, phsB( xidxB ), 'b', xidxB, phsA( xidxA ), '--r' );
legend( ph, { 'block', 'full' } )
xlim( xidxB( [ 1 end ] ) )
set( gca, 'tickdir', 'out', 'box', 'off' )
title( 'Phase [rad]' )
xlabel( 'Within-block [sample]' )

str = input( sprintf( 'Block $%d; continue?', bidx ), 's' );
rc                              = 0;
switch str
    case 'k'
        keyboard
    case 'q'
        rc                      = 1;
    otherwise
end

return

% EOF
