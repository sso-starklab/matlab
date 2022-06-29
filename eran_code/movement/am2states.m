% AM2STATES             segment time into movement/immobility by AM recording
%
% [ mat0, mat1, am ] = am2states( filebase, channels, nchans, Fs, TH, graphics, Overwrite )
%
% filebase      full path to files
% channels      to consider as AM. defaults to 3 last channels
% nchans        in eeg file. defaults to xml file parameter nChannels
% Fs            of eeg file. defaults to xml file parameter lfpSampleRate
% TH            to use for partitioning data into states
%                   defaults to 0.03/0.05G
%                   negative to recompute
% Overwrite     {0}/1/-1/-2; applies only to am file
%
% mat0          immobility segments
% mat1          mobility segments
% am            RMS of vector sum of AM (see am2rms)
%
% does:
% (1) get RMS of resultant vector (compute and save, or load)
%   -load the raw AM data for the file
%   -scale to G
%   -compute RMS (see AM2RMS.M) + saves into *.am file
% (2) determine classification TH (optional)
%   -determine threshold by the dip of the bimodal distribution (should be
%   around 0.01-0.2) or by an external TH (default: 0.05 m/sec^2; this is conservative)
% (3) segment
%   -segment the data into immobility/mobility (95% of the values in a 1sec moving window below/above TH)
%
% files
% input:        filebase.eeg
% optional      filebase.xml
% output        filebase.am
%
% calls         ParseArgPairs, LoadXml, get_stimchans
%               am2rms
%               local_max
%               firfilt, parse, dilutesegments, setdiffranges
%               replacetok, getdatainranges, plotranges

% 16-jul-12 ES

% revisions
% 17-jul-12 saving *.am file added
% 08-nov-12 (1) check for eeg file before going on
%           (2) faster loading methods using memmapfile
% 19-nov-12 (1) external arguments, including external TH
%           (2) always convert to G
%           (3) if TH is two-element vector, use different threshold for
%               immobility/mobility
% 02-sep-13 default amchans from *prm.xml file
% 31-aug-18 (1) modified to support externally-given bias and V2G
%           (2) get the P2P and nBits from par
%           (3) save binary locally
%           (4) TH determined using log-scale binning - gives nice trough
%           (5) graphics modified, saved
%           (6) TH output argument
% 10-sep-18 (1) modified amTH algorithm
%           (2) modified argument call
% 14-oct-19 (1) handled case of continuous mobility 

% to do:
% blockwise.

function [ mat0, mat1, am, TH ] = am2states( filebase, varargin )

%--------------------------------------------------------------------%
% preps
%--------------------------------------------------------------------%
% scaling constants (irrelevant if adaptive TH)
DEFAULT_biasV               = 1.5;                                                                  % 3V supply
DEFAULT_V2G                 = 1 / 0.3;                                                              % scaling of ~300mV/g
DEFAULT_TH                  = [ 0.02 -1 ];

% arguments
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    error( 'filebase required' )
end
[ channels, nchans, Fs, TH, graphics...
    , Overwrite, WINDUR, FRACTVALS, mindurSEC, minisiSEC...
    , amZ, amMA, amLP, amT...
    , biasV, V2G ] = ParseArgPairs(...
    { 'channels', 'nchans', 'Fs', 'TH', 'graphics' ...
    , 'Overwrite', 'WINDUR', 'FRACTVALS', 'mindurSEC', 'minisiSEC' ...
    , 'amZ', 'amMA', 'amLP', 'amT'...
    , 'biasV', 'V2G' }...
    , { [], [], [], DEFAULT_TH, [ 1 0 ]...
    , 0, 1, 0.95, 0.1, 0.5...                                                                       % logical, s, fractions, s, s
    , 0, 2, 10, 0.05 ...                                                                            % logical s, Hz, s
    , DEFAULT_biasV, DEFAULT_V2G } ...
    , varargin{ : } );

% parameters
if length( graphics ) == 1
    graphics                = [ graphics 0 ];
end
minTH                       = [];
if TH( 2 ) <= 0
    minTH                   = TH( 1 );
    TH                      = [];
end
if length( TH ) > 2
    TH                      = TH( 1 : 2 );
end
xmlfname                    = [ filebase '.xml' ];
eegfname                    = [ filebase '.eeg' ];
amfile                      = [ filebase '.am' ];
[ ~, fname ]                = fileparts( amfile );
cpathi                      = strfind( filebase, 'dat' );
cpath                       = filebase( 1 : ( cpathi - 1 ) );
fpath                       = sprintf( '%sfigs/', cpath );
if ~exist( fpath, 'dir' )
    eval( sprintf( '!mkdir %s', fpath ) )
end
figfname                    = sprintf( '%s%s.am', fpath, fname );
if ( isempty( nchans ) && ~exist( xmlfname, 'file' ) ) || ~exist( eegfname, 'file' )
    error( 'missing files' )
end
if isempty( nchans ) || isempty( Fs )
    par                     = LoadXml( xmlfname );                                                  % default to xml file parameters
    nchans                  = par.nChannels;
    Fs                      = par.lfpSampleRate;
end
if isempty( channels )
    if ~exist( 'par', 'var' )
        par                 = LoadXml( xmlfname );
    end
    channels                = get_stimchans( par, [], 'am' );
    if isempty( channels )
        mat0                = [];
        mat1                = [];
        am                  = NaN;
        return
    end
    fprintf( 1, '%s: %s: Using channels %s\n', mfname, filebase, num2str( channels( : ).' ) );
end
if ~exist( 'par', 'var' )
    par                     = LoadXml( xmlfname );
end
nBits                       = par.nBits;
try
    [ ~, ~, voltageRanges ] = get_stimchans( par, [], 'am' );
    P2P                     = mean( voltageRanges );
catch
    P2P                     = par.VoltageRange;
end

%--------------------------------------------------------------------%
% load am data and compute rms vector
%--------------------------------------------------------------------%
if Overwrite == 1 || Overwrite == -1 || ~exist( amfile, 'file' )
    % load data
    fprintf( 1, '%s: loading eeg data...', mfname );
    a                       = memmapfile( eegfname, 'Format', 'int16' );
    n                       = length( a.data );
    idx                     = ( channels( : ) * ones( 1, n / nchans ) + ones( length( channels ), 1 ) * [ 0 : nchans : ( n - nchans ) ] )';
    x                       = single( a.data( idx ) );
    clear a
    % compute rms
    x                       = ( ( x / 2.^nBits * P2P ) - biasV ) * V2G;                             % convert bits->g (about +-5)
    fprintf( 1, 'computing RMS...' );
    am                      = am2rms( x, Fs, graphics( 2 ), amZ, amMA, amLP, amT );
    clear x
    % save to file
    if Overwrite == 1 || ( ~exist( amfile, 'file' ) && Overwrite ~= -1 )
        fprintf( 1, 'saving BINARY %s\n', amfile );
        fp                  = fopen( amfile, 'w' );
        fwrite( fp, single( am ), 'float32' );
        fclose( fp );
    end
else
    fprintf( 1, '%s: loading %s\n', mfname, amfile );
    a                       = memmapfile( [ filebase '.am' ], 'Format', 'single' );
    am                      = single( a.data );
    clear a
end

%--------------------------------------------------------------------%
% determine THs if not given
%--------------------------------------------------------------------%

% find the first real trough in the smoothed (bimodal) distribution
if isempty( TH ) || graphics( 1 )
    fprintf( 1, '%s: determining threshold..\n', mfname );
    nbins                   = 200; 
    bins                    = logspace( log10( min( am ) ), log10( max( am ) ), nbins );
    h                       = hist( am, bins );
    win                     = triang( 11 );
    win                     = win / sum( win );
    h                       = firfilt( h, win );
    val                     = local_max( h, 'min' );
end
if isempty( TH )
    val( bins( val ) < minTH ) = [];
    if isempty( val )
        val                 = find( bins >= minTH, 1, 'first' );
    end
    if length( val ) == 1
        TH                  = bins( val );
    else
        if val( 1 ) == 1                                                                            % take the second minima
            TH              = bins( val( 2 ) );
        else                                                                                        % take the first minima
            TH              = bins( val( 1 ) );
        end
    end
    TH                      = TH * [ 0.95 1.05 ];
end

%--------------------------------------------------------------------%
% partition to states
%--------------------------------------------------------------------%

% look for time segments that the animal was mainly im/mobile
fprintf( 1, '%s: post-processing..', mfname );
win                         = ones( round( WINDUR * Fs ), 1 ) / Fs;
if length( TH ) == 1                                                                                % single TH
    im                      = firfilt( am <= TH, win ) / WINDUR;
    mat0                    = parse( find( im >= FRACTVALS ) );                                     % immobile
    mat1                    = parse( find( im <= ( 1 - FRACTVALS ) ) );                             % mobile
elseif length( TH ) == 2                                                                            % two THs
    mat0                    = parse( find( am <= TH( 1 ) ) );
    mat1                    = parse( find( am >= TH( 2 ) ) );
end

% apply post-hoc logic (combine/reject short events):
mat0                        = dilutesegments( mat0, mindurSEC * Fs, minisiSEC * Fs );
mat1                        = dilutesegments( mat1, mindurSEC * Fs, minisiSEC * Fs );

% make sure mutually exclusive
mat1tmp                     = setdiffranges( mat1, mat0 );
mat0                        = setdiffranges( mat0, mat1 );
mat1                        = mat1tmp;

% validate immobility, >95% of values/segment are <TH( end )
if ~isempty( mat0 )
    amSeg                   = zeros( size( mat0 ) );
    for i                   = 1 : size( mat0, 1 )
        amSeg( i, : )       = [ mean( am( mat0( i, 1 ) : mat0( i, 2 ) ) ) std( am( mat0( i, 1 ) : mat0( i, 2 ) ) ) ];
    end
    nstd                    = abs( norminv( min( 1 - FRACTVALS, FRACTVALS ) / 2, 0, 1 ) );
    hi                      = ( amSeg( :, 1 ) + nstd * amSeg( :, 2 ) ) > TH( end );
    mat0( hi, : )           = [];
end
fprintf( 1, 'done!\n' );

%--------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%
if graphics( 1 )
    
    tidx                    = 1 : size( am, 1 );
    tim                     = tidx / Fs;
    [ ~, filename ]         = fileparts( filebase );
    fig                     = figure;
    
    subplot( 3, 1, 1 )
    plot( tim( 1 : 10 : end ), am( tidx( 1 : 10 : end ) ) ),
    xlim( tim( [ 1 end ] ) )
    title( sprintf( '%s: RMS( AM )', replacetok( filename, '\_', '_' ) ) )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    subplot( 3, 1, 2 ),
    hold on
    plot( bins, h ),
    line( TH( 1 ) * [ 1 1 ], ylim, 'color', [ 1 0 0 ] );
    if length( TH ) == 2
        line( TH( 2 ) * [ 1 1 ], ylim, 'color', [ 0 0 1 ] );
    end
    title( sprintf( 'Distribution of RMS(AM); TH=%0.3g/%0.3g m/s^2 %0.3g%%/%0.3g%%' ...
        , TH( 1 ), TH( end ), 100 * sum( am <= TH( 1 ) ) / length( am ), 100 * sum( am >= TH( end ) ) / length( am ) ) );
    xlim( bins( find( ( cumsum( h ) / sum( h ) ) > 0.99, 1, 'first' ) ) * [ -0.01 1 ] )
    h0                      = hist( getdatainranges( am, mat0 ), bins );
    h1                      = hist( getdatainranges( am, mat1 ), bins );
    b1                      = bar( bins( : ), h1( : ) ); set( b1, 'edgecolor', [ 0 0 1 ], 'facecolor', [ 0 0 1 ] );
    b0                      = bar( bins( : ), h0( : ) ); set( b0, 'edgecolor', [ 1 0 0 ], 'facecolor', [ 1 0 0 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    subplot( 3, 1, 3 ),
    hold on
    if exist( 'im', 'var' )
        plot( tim( 1 : 10 : end ), im( tidx( 1 : 10 : end ) ), '-r' )
    end
    if exist( 'mo', 'var' )
        plot( tim( 1 : 10 : end ), mo( tidx( 1 : 10 : end ) ), '-b' )
    end
    if ~isempty( mat0 )
        ph                  = plotranges( ceil( mat0 / Fs ), 1.1 );
        set( ph, 'linewidth', 2, 'color', [ 1 0 0 ] )
    end
    if ~isempty( mat1 )
        ph                  = plotranges( ceil( mat1 / Fs ), 1.2 );
        set( ph, 'linewidth', 2, 'color', [ 0 0 1 ] )
    end
    hold off
    title( 'im/mobility index' )
    xlabel( 'Time [s]' )
    xlim( tim( [ 1 end ] ) )
    ylim( [ 0 1.3 ] );
    set( gca, 'tickdir', 'out', 'box', 'off' )
    
    % save the figure
    fprintf( '%s: Saving %s\n', mfname, figfname )
    print( fig, '-dpng', [ figfname '.png' ] )
    
end


return

% EOF
