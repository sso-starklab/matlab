% plotProbeHFOs     summary plot of LFP and CSD in time-space and time-freq domains
% 
% call              [ x, trigs ] = plotProbeHFOs( filebase, toflip, ignoredChannels, varargin  )
%
% gets              filebase            full path
%                   toflip              {0}, flip channels (to account for convention mismatchs)
%                   ignoredChannels     {[]}
%                   
% optional arguments
%                   plotSpec            {1}
%                   savetype            {'png'}
%                   ripbase             {'spw'}
%                   xPostInj            1
%                   staWin              {[ -0.05 0.05 ]}
%                   padBuffer           {[ -0.01 0.01 ]}
%                   toSmooth            1
%
% calls
%
%   preparations:   ParseArgPairs, LoadXml, get_probe, LoadStims, LoadVals, resampleranges, get_egroup, isoverlap
%   computations:   pt_avg
%   plots:          tilefig, makegaussfir, firfilt, imupsample, plotTraces, plotSpectrogram, calibration
%   saving:         replacetok, fig_out
%
% returns           x                   lfp time domain matrices (nsamples x nchannels x ngroups)
%                   trigs               event time vectors (cell array: ngroups x 1)

% 11-oct-13 ES

% revisions
% 11-nov-13 padBuffer added
% 12-dec-14 optional smoothing added
% 21-jul-15 (1) case of nsites == 1 handled...
%           (2) case of par.SpkGrps missing handled...\
% 31-aug-18 (1) cleaned
% 08-jul-20 (1) modified so that probe and egroups will be based on
%               Anatomical groups (and not spike groups)
% 19-jul-20 (1) changed call to get_egroup (that was not changed in 08-jul)
%               so that probe and egroups will be based on Anatomical groups

function [ x, trigs ] = plotProbeHFOs( filebase, toflip, ignoredChannels, varargin  )

%----------------------------------------------------------------------%
% preps
%----------------------------------------------------------------------%

% constants
tempSD                      = 0;
spatBin                     = 1;

% arguments
nargs                       = nargin;
if nargs < 2 || isempty( toflip )
    toflip                  = 0;
end
if nargs < 3 || isempty( ignoredChannels )
    ignoredChannels         = [];
end

[ plotSpec, savetype, ripbase, xPostInj, staWin, padBuffer, toSmooth ] = ParseArgPairs(...
    { 'plotSpec', 'savetype', 'ripbase', 'xPostInj', 'staWin', 'padBuffer', 'toSmooth' }...
    , { 1, 'png', 'spw', 1, [ -0.05 0.05 ], [ -0.01 0.01 ], 1 }...
    , varargin{ : } );


if plotSpec == 1
    avgMode                 = 'wlt';
else
    avgMode                 = 'eeg';
end

% paths
delim                       = strfind( filebase, '/dat/' );
if isempty( delim )
    fprintf( '%s: Cannot save fig\n', upper( mfilename ) )
end
figdir                      = [ filebase( 1 : delim ) 'figs/hfo' ];
if ~exist( fileparts( figdir ), 'dir' )
    mkdir( fileparts( fileparts( figdir ) ), 'figs' )
end
if ~exist( figdir, 'dir' )
    mkdir( fileparts( figdir ), 'hfo' )
end
[ pathname, filename, extname ] = fileparts( filebase );
filename = [ filename extname ];
if strcmp( ripbase, 'spw' )
    figname                 = [ figdir '/' filename '.hfo_spontaneous' ];
else
    figname                 = [ figdir '/' filename '.' ripbase '_spontaneous' ];
end

% par
par                         = LoadXml( filebase );
try
    ngroups                 = length( par.SpkGrps );
catch
    ngroups                 = 1;
end
probe                       = get_probe( par, 0 );
if isempty( probe )
    probe                   = 1;
end
nsites                      = size( probe, 1 );
if ngroups > size( probe, 2 )
    error( 'mismatch' )
end
scalefactor                 = par.VoltageRange / 2^par.nBits / par.Amplification * 1e3; % A2DU to mV (not uV)
if toflip
    probe                   = flipud( probe );
end

% stim
[ Vals, Trigs ]             = LoadStims( filebase );
if isempty( Vals )
    [ Vals Trigs ]          = LoadVals( filebase );
end
if ~isempty( Vals )
    vals                    = resampleranges( Vals( :, 1 : 2 ), par.lfpSampleRate, par.SampleRate );
    pad                     = [ floor( padBuffer( 1 ) * par.lfpSampleRate ) ceil( padBuffer( 2 ) * par.lfpSampleRate ) ];
    vals                    = [ vals( :, 1 ) + pad( 1 ) vals( :, 2 ) + pad( 2 ) ];
else
    vals                    = [];
end

% exclude post-injection data
if xPostInj
    ichan                   = get_stimchans( par, [], 'nanoject' );
    if ~isempty( ichan )
        evtfname            = [ filebase '.evt.t' num2str( ichan ) ];
    else
        evtfname            = [];
    end
    if exist( evtfname, 'file' )
        iTimes              = load( evtfname );
        a                   = memmapfile( [ filebase '.eeg' ], 'Format', 'int16' );
        nsamples            = length( a.Data ) / par.nChannels;
        clear a
        xperiods            = [ round( min( iTimes( : ) ) / 1000 * par.lfpSampleRate ) nsamples ];
        vals                = sortranges( [ vals; xperiods ] );
    end
end

% ripple event times
load( [ filebase '.sps' ], '-mat', 'vote' );
trigs                       = cell( ngroups, 1 );
chans                       = cell( ngroups, 1 );
for i                       = 1 : ngroups
    chan                    = vote( i );
    if isnan( chan )
        continue
    end
    [ egroup, ~, chans{ i } ]  = get_egroup( par, chan, 0 );
    if isempty( egroup )
        error( 'misamtch' )
    end
    chans{ i }( ismember( chans{ i }, ignoredChannels ) ) = [];
    spwfname                = sprintf( '%s.%s.%s', filebase, ripbase, num3str( chan ) );
    if ~exist( spwfname, 'file' )
        continue
    end
    fprintf( 'C%d: Loading HFO times...\n', chan );
    load( spwfname, '-mat', 'rips' );
    if ~isequal( chan, rips.chans )
        error( 'mismatch' )
    end
    if ~isempty( vals )
        ridx                = isoverlap( rips.edges, vals );
        rips                = rips_select( rips, ~ridx );
    end
    trigs{ i }              = rips.trigs/rips.Fs;
end


%----------------------------------------------------------------------%
% computations
%----------------------------------------------------------------------%

% get the mean LFP/wavelets:
stawins                     = [ ceil( staWin( 1 ) * par.lfpSampleRate ) floor( staWin( 2 ) * par.lfpSampleRate ) ];
nsamples                    = diff( stawins ) + 1;
x                           = NaN * ones( nsamples, nsites, size( probe, 2 ) );
y                           = x;
yoo                         = cell( 1, ngroups );
yooCSD                      = yoo;
for i                       = 1 : ngroups
    chan                    = vote( i );
    if isnan( chan )
        continue
    end
    if isempty( trigs{ i } ) || isempty( chans{ i } )
        continue
    end
    [ avgcsd, avglfp, tim, ~, ~, too, foo, yoo0, yoo1 ] = pt_avg( par...
        , chans{ i }, trigs{ i }, 'specChans', vote( i )...
        , 'stawinSEC', staWin, 'tempSD', tempSD, 'spatBin', spatBin...
        , 'suffix', avgMode, 'graphics', 0 ...
        , 'scalefactor', scalefactor );
    if length( probe( : ) ) == 1
        ii                  = 1;
        uj                  = 1;
    else
        [ ii, jj ]          = find( ismember( probe, chans{ i } ) );
        uj                  = unique( jj );
    end
    x( :, ii, uj )          = avglfp;
    if length( chans{ i } ) >= 3
        y( :, ii, uj )      = avgcsd;
    end
    if plotSpec
        yoo{ i }            = nanmean( yoo0 + eps, 3 );
        if ~isempty( yoo1 )
            yooCSD{ i }     = nanmean( yoo1 + eps, 3 );
        end
    end
end

if toflip
    x                       = x( :, nsites : -1 : 1, : );
    y                       = y( :, nsites : -1 : 1, : );
end

%----------------------------------------------------------------------%
% plot
%----------------------------------------------------------------------%

if plotSpec
    [ ah, fig ]             = tilefig( 2, ngroups, -2, 0.85, 'right' );
else
    [ ah, fig ]             = tilefig( 2, ngroups, 1, 0.85, 'right' );
end

if toSmooth
    gwin                    = makegaussfir( toSmooth, 1 );
end

% LFP, CSD. each divided into time-freq, time-space
USF                         = 10;
iNaN                        = 1; % how to treat missing channels in time-space plots

% parameters
yCalib                      = 0.1;                              % mV
xCalib                      = diff( staWin ) * 0.1 * 1000;      % 10% of the time range;
calib                       = [ xCalib yCalib ];
calibf                      = 50;
CALIBCOLOR                  = [ 0 0.7 0 ];
CALIBWIDTH                  = 2;

clims1                      = NaN * ones( ngroups, 2, 2 );
clims2                      = clims1;

% go over groups/shanks
for i                       = 1 : ngroups
    
    % decide whether to flip
    if toflip
        fchans              = fliplr( chans{ i }( : ).' );
    else
        fchans              = chans{ i };
    end
    
    for j                   = 1 : 2 % LFP then CSD
        
        % time-space plots (traces)
        subplot( ah( i + ( j - 1 ) * ngroups, 1 ) )
        nsites              = length( fchans );
        if j == 1
            xplot           = x;
            str             = 'LFP';
            ytick           = 1 : nsites;
            if nsites == 1 && size( xplot, 2 ) > 1
                xplot       = xplot( :, chans{ i } );
            end
            scaleTraces     = 1;
        else
            xplot           = y;
            str             = 'CSD';
            ytick           = 1 : nsites;
            scaleTraces     = 2;
        end
        if isempty( ytick ) || sum( sum( isnan( x( :, :, i ) ) ) ) == numel( xplot( :, :, i ) ) ...
                || j == 2 && nsites < 2 
            axis off
            subplot( ah( i + ( j - 1 ) * ngroups, 2 ) )
            axis off
            continue
        end
        if toSmooth > 0
            xplotI          = firfilt( xplot( :, :, i ), gwin );
        else
            xplotI          = xplot( :, :, i );
        end
        if nsites == 1
            plot( tim, xplotI );
        else
            [ x1, y1, z1 ]  = imupsample( tim, ytick, xplotI( :, ytick )', USF, iNaN );
            imagesc( x1, y1, z1 ), axis xy
            clims1( i, :, j ) = get( gca, 'clim' );
            if i == 1
                if ~plotSpec
                    xlabel( 'Time [ms]' )
                end
                ylabel( [ str ' channel' ] )
            end
            set( gca, 'xticklabel', '' )
            line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
            hold on
            [ ~, ph3 ]      = plotTraces( tim, xplotI, scaleTraces, [ 1 nsites 0 nsites + 1 ], calib );
            ylim( [ 1 nsites ] + 0.5 / USF * [ -1 1 ]  )
            set( ph3, 'color', [ 1 1 1 ] * 0.3 );
            tomark  = fchans == vote( i );
            set( ph3( tomark ), 'color', [ 1 0 0 ], 'linewidth', 2 );
            set( gca, 'yticklabel', probe( :, i ) )
            set( gca, 'ytick', ytick, 'yticklabel', fchans( ytick ) )
        end
        axis square
        set( gca, 'tickdir', 'out', 'box', 'off' )
        
        % time-freq (wavelets)
        if plotSpec
            
            subplot( ah( i + ( j - 1 ) * ngroups, 2 ) )
            if j == 1
                xplot       = yoo{ i };
            else
                xplot = yooCSD{ i };
                if isempty( xplot )
                    axis off
                    continue
                end
            end
            
            plotSpectrogram( too * 1000, foo, xplot', USF, 'contour', 'linear' );
            alines( [ 80 160 240 ], 'y', 'color', [ 1 1 1 ], 'linestyle', '--' );
            if i == 1 && j == 2
                ylabel( 'Freq (Hz)' )
                xlabel( '' )
                set( gca, 'xticklabel', '' )
            else
                xlabel( '' )
                ylabel( '' )
                set( gca, 'xticklabel', '', 'yticklabel', '' )
            end
            title( '' )
            clims2( i, :, j ) = get( gca, 'clim' );
            line( [ 0 0 ], ylim, 'color', [ 1 1 1 ] * 1, 'linewidth', 1, 'linestyle', '--' );
            calibration( [ calib( 1 ) calibf ], 1i* 0.05 * [ 1 1 ], gca, 0, 'color', CALIBCOLOR, 'linewidth', CALIBWIDTH );
            axis square
            ah0 = gca;
            set( ah0, 'tickdir', 'out', 'box', 'off' )
            
            % add a colorbar
            pos0 = get( ah0, 'position' );
            h = colorbar( 'horiz' );
            set( h, 'tickdir', 'out', 'box', 'off' )
            set( ah0, 'position', pos0 )
            posh = get( h, 'position' );
            newposh = [ posh( 1 ) max( 0, pos0( 2 ) - 1.5 * posh( 4 ) / 2 ) posh( 3 ) posh( 4 ) / 2 ];
            set( h, 'position', newposh )
            
        end % plotSpec
        
    end % LFP/CSD
    
end % groups

% scale
for i                       = 1 : ngroups
    for j                   = 1 : 2
        % time-space (not all same scale)
        subplot( ah( i + ( j - 1 ) * ngroups, 1 ) )
        zstr                = sprintf( '[%d %d]'...
                            , round( 1000 * clims1( i, 1, j ) )...
                            , round( 1000 * clims1( i, 2, j ) ) );
        th                  = text( min( xlim ) + diff( xlim ) * 0.9, min( ylim ) + diff( ylim ) * 0.9, zstr );
        set( th, 'horizontalAlignment', 'right', 'color', [ 1 0 1 ], 'fontsize', 9 )
        % LFP
        if j == 1 
            title( sprintf( 'C%d; %d'...
                , vote( i ), length( trigs{ i } ) ) )
        end
    end
end
tstr                        = sprintf( '%s', replacetok( filename, '\_', '_' ) );
textf( 0.5, 0.975, tstr );

%----------------------------------------------------------------------%
% save
%----------------------------------------------------------------------%
if ~isempty( savetype ) && ( isa( savetype, 'cell' ) || ~all( isnan( savetype ) ) ) && ~isempty( figname )
    if ~isa( savetype, 'cell' )
        savetype            = { savetype };
    end
    for j                   = 1 : length( savetype )
        for i               = 1 : length( fig )
            if fig( i ) == 0
                continue
            end
            fig_out( fig( i ), 1, [ figname '.' savetype{ j } ], savetype{ j } );
        end
    end
end

return

% EOF