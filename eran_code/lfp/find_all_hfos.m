% find_all_hfos       detect channels with max amp spontaneous HFOs and compute stats + figs for those
%
% [ stats ] = find_all_hfos( filebase, Overwrite )
% 
% filebase      full file base or 2-element cell array
% Overwrite:    {1} to recompute and overwrite
%               0 to just compute (with writing but not overwriting)
%               -1 to just load/compute if a file exists/doesn't (no writing)
%               -2 load if exists, compute and save if doesn't
%
% stats         1 row per electrode group (shank): 
%               [ number best_channel nrips occurrence_rate details ], 
%               where details are mean: [ f pow sd ncycles durs ]
%
% notes:
%   (1) it is not ideal to run this (or any state-dependent analyses) on
%       individual files, because then SWS detection (and TH determination) is noisy
%   (2) although the detection is based on the SWS statistics, this routine
%       considers spontaneous events ONLY (no overlap with any stimulus) during any brain state
%       if the total duration of the stimulus is short this assumption is reasonable
%
% post-hoc verification:
%   the channels are chosen by the high-pass filtered mean LFP, but a
%   well-localized CSD source (in time-space) or spectral event (in
%   time-freq) is not required. Thus, this info can be used for indpendent
%   verification (look at the png)
%
% calls: 
%   ParseArgPairs, LoadXml, LoadStims, LoadVals, get_egroup                     (formats)
%   find_hfos, rips_select, segmentBehavior                                     (lfp)
%   verb, uhist                                                                 (general)
%   dilutesegments, intersectranges, isoverlap, resampleranges, setdiffranges   (sets) 
%   firfilt, makegausslpfir, pt_avg                                             (ssp)
%
% files:
%   input: *eeg, *xml, (*whl)
%   intermediate:   segmentBehavior: *phs, *am/*mov.mat, *sts*
%                   detect_hfos: *.spw.NNN, *.evt.rNN (optional)
%   final: *wltBL.NNN, *.spw.NNN, *.evt.rNN - for the selected channels
%          *sps (statistics)

% 07-nov-12 ES

% revisions
% 11-dec-12 argument handling, behavioral segmentation, val file merging..
% 07-feb-13 (1) added buffer around stim detections
%           (2) use LoadStims instead of LoadVals
% 30-mar-13 constants moved to varargin
% 03-apr-13 added prior. This is needed when two oscillators are
%               included in the same electrode group (e.g. linear probe spanning CA3 and
%               CA1)
% 07-apr-13 (1) overwrite individual spw and wltBL files if Overwrite is 1
%           (2) modified selection algorithm 
% 18-sep-13 (1) median filtering added
%           (2) bug in second iteration corrected (did not compute baseline
%           based on SWS in final computations)
% 01-oct-13 (1) changes in algorithm (see detect_hfos) and parameters (5
%               SD)
%           (2) ignoredChannels added
% 23-oct-13 (1) support of external bperiods
%           (2) imporved algorithms for channel selection
%           (3) manual exceptions added..
% 13-aug-19 cleaned a bit and renamed as find_hfos.m
% 18-aug-19 cleaned up
% 15-may-20 added optional argument bstate (defaults to 'sws')
% 14-jul-20 changed automated channel selection to rely on Anatomical
%           groups and not spike groups (see get_egroup)
%           this allows supporting *xml files in which a spike group is
%           only a subset of an anatomical groups

% DETAILED ALGORITHM:
% 1. Make sure data are segmented into states (segmentBehavior)
%       for this, a theta signal is selected, theta/delta is computed,
%       acceleration/movment is computed, and data are partitioned into states
% 12 Detect HFOs independently on each and every channel. 
%       calls: detect_hfos.
%       results (optional): *.spw.NNN, *.evt.rNN
% 2. For each electrode group, determine the channel for which the detected events are of maximal amplitude
% 3. Compute the statistics and plot CSD/spectrograms for the selected channels
% 4. Compute stats

function [ stats ] = find_all_hfos( filebase, varargin )

stats                       = [];

%------------------------------------------------------------------%
% globals
%------------------------------------------------------------------%
global BatchMode
if isempty( BatchMode )
    BatchMode               = 0;
end

%------------------------------------------------------------------------%
% constants
%------------------------------------------------------------------------%
% segmentBehavior:
behaviorDurSEC              = 2;             % [s]
OverwriteStates             = -2;           % do not overwrite if already existing
minSWS                      = 10;                    % [s] - minimum SWS duration. otherwise ALL data are used as a baseline!!

% general
vflag                       = 1;

%------------------------------------------------------------------------%
% arguments
%------------------------------------------------------------------------%
mfname                      = upper( mfilename );
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ Overwrite, normalizeSWS, ignoredChannels, iperiods...
    , bperiods, bstate ...
    , padBuffer, ripBP, ripTH, mindurSECrip, minisiSECrip, diffOrd...
    , filtMode, powerMode, clipBase...
    , OverwriteSPWsearch, bpROI, OverwriteSPW, prior, selectionMode...
    , staWin, tempSD, spatBin, whitenFlag, savetype, graphics ] = ParseArgPairs(...
    { 'Overwrite', 'normalizeSWS', 'ignoredChannels', 'iperiods'...
    , 'bperiods', 'bstate' ...
    , 'padBuffer', 'ripBP', 'ripTH', 'mindurSECrip', 'minisiSECrip', 'diffOrd'...
    , 'filtMode', 'powerMode', 'clipBase'...
    , 'OverwriteSPWsearch', 'bpROI', 'OverwriteSPW', 'prior', 'selectionMode'...
    , 'staWin', 'tempSD', 'spatBin', 'whitenFlag', 'savetype', 'graphics' }...
    , { -2, 1e3, [], [] ...
    , [], 'sws' ...
    , [ -0.01 0.01 ], [ 80 250 ], [ 2 5 15 ], 0.015, 0.015, 0 ...
    , 'dog', 'LP', 1 ...
    , 0, [], 1, [], 'lfpPower' ...
    , [ -0.05 0.05 ], 0.0001, 1, 0, 'png', 1 }...
    , varargin{ : } );

% get the filebase and the par file
par                         = LoadXml( filebase );
[ pathname, filename, extname ] = fileparts( filebase );
filename                    = [ filename extname ];
if ~exist( pathname, 'dir' )
    fprintf( 1, '%s: missing directory %s\n', mfname, pathname )
    return
end
Fs                          = par.lfpSampleRate;
if isempty( bpROI )
    bpROI                   = [ -1 1 ] * ceil( 1/mean( ripBP ) * Fs * pi / 2 ); % pi cycles (was +-20 ms, which is way too much)
end

% directory for figures
homedir                     = strfind( pathname, 'dat' );
if isempty( homedir )
    homedir                 = [ fileparts( pathname ) '/' ];
else
    homedir                 = pathname( 1 : homedir - 1 );
end
figdir                      = sprintf( '%sfigs', homedir );
if ~exist( figdir, 'dir' )
    mkdir( homedir, 'figs' )
end

% determine whether to run
spsfname                    = sprintf( '%s.sps', filebase );
if Overwrite < 0 && exist( spsfname, 'file' )
    verb( sprintf( '%s: loading %s...', mfname, spsfname ), vflag )
    L = load( spsfname, '-mat' );
    if Overwrite <= -1 || graphics <= 0
        stats               = L.stats;
        vote                = L.vote;
        if graphics <= 0
            return
        end
    end
end
if Overwrite == -1
    figdir                  = 1; 
end

% get the val/evt files
Vals                        = LoadStims( filebase );
if isempty( Vals )
    Vals                    = LoadVals( filebase );
end

% make sure segmented into behavioral states
segmentBehavior( filebase, 'windur', behaviorDurSEC, 'Overwrite', OverwriteStates );
verb( sprintf( '%s: loading sts.sws file...', mfname ), vflag )
if isempty( bperiods )
    bperiods                = load( [ filebase '.sts.' bstate ] );
end
totalSWS                    = sum( diff( bperiods, 1, 2 ) + 1 ) / par.lfpSampleRate;
a                           = memmapfile( [ filebase '.eeg' ], 'Format', 'int16' );
neeg                        = length( a.data ) / par.nChannels;
clear a
if totalSWS < minSWS
    fprintf( '\n\n%s: NOTE: total time in SWS is %0.3g sec; may use entire dataset for baseline!!\n\n\n'...
        , upper( mfilename ), totalSWS )
    disp( [ 1 neeg ] )
end
bperiods                    = intersectranges( bperiods, [ 1 neeg ] );
if ~isempty( iperiods )
    bperiods                = intersectranges( bperiods, iperiods );
end

% stimulus times
if isempty( Vals )
    verb( sprintf( '%s: no stimuli for %s', mfname, filebase ), vflag )
else
    verb( sprintf( '%s: acquiring stimulus times...', mfname ), vflag )
    vals                    = resampleranges( Vals( :, 1 : 2 ), par.lfpSampleRate, par.SampleRate );
    pad                     = [ floor( padBuffer( 1 ) * par.lfpSampleRate ) ceil( padBuffer( 2 ) * par.lfpSampleRate ) ];
    vals                    = [ vals( :, 1 ) + pad( 1 ) vals( :, 2 ) + pad( 2 ) ];
    bperiods                = setdiffranges( bperiods, vals );
    bperiods                = dilutesegments( bperiods, behaviorDurSEC * par.lfpSampleRate, 0 );
end

% parameters for site/shank selection:
gwin                        = makegausslpfir( mean( ripBP ) / pi, Fs, 4 );
winsamp                     = [ ceil( staWin( 1 ) * Fs ) floor( staWin( 2 ) * Fs ) ];
cidx                        = abs( winsamp( 1 ) ) + 1  + round( bpROI * Fs / 1000 );
cidx                        = cidx( 1 ) : cidx( 2 ); % ideally, should use a tapered window here..

%------------------------------------------------------------------------%
% search stage
%------------------------------------------------------------------------%

% initialize
ngroups                     = length( par.SpkGrps );
%ngroups                     = length( par.AnatGrps );
% get_stimchans( par, [], 'neuronal' )
if exist( 'vote', 'var' )
    search                  = 0;
else
    search                  = 1;
    vote                    = NaN * ones( 1, ngroups );
end
%[ ~, allchans ]             = get_egroup( par );
[ ~, ~, allchans ]           = get_egroup( par );
allchans( ismember( allchans, ignoredChannels ) ) = [];
mat                         = zeros( length( allchans ), 7 );
mat( :, 1 )                 = allchans( : );
k                           = 0;

% prior
if isempty( prior )
    prior                   = [ ( 1 : par.nChannels )' ones( par.nChannels, 1 ) ];
end

% go over groups
for i                       = 1 : ngroups
    if ~search
        break
    end
    %chans                   = par.SpkGrps( i ).Channels + 1;
    chans                   = par.AnatGrps( i ).Channels + 1;
    chans( ismember( chans, ignoredChannels ) ) = [];
    if isempty( chans ) 
        continue
    end
    if sum( ismember( prior( :, 1 ), chans ) ) == length( chans )
        pr                  = prior( ismember( prior( :, 1 ), chans ), 2 );
    else
        pr                  = ones( length( chans ), 1 );
    end
    nchans                  = length( chans );
    for j                   = 1 : nchans
        k                   = k + 1;
        mat( k, 2 )         = i;
        chan                = chans( j );        
        spwfname            = sprintf( '%s.spw.%s', filebase, num3str( chan ) );
        verb( sprintf( '%s: C%d: ', mfname, chan ), -1 );
        if OverwriteSPW < 0 && exist( spwfname, 'file' )
            verb( sprintf( 'loading HFO times...' ), -1 );
            load( spwfname, '-mat' );
            reDo = 0;
            if reDo == 0
                % difference between span of troughs and span of entire ripple, max two (2.5..) cycles, i.e.
                dt          = zeros( length( rips.t ), 2 );
                for h       = 1 : length( rips.t )
                    dt( h, : ) = minmax( rips.peaks( ismember( rips.peaks( :, 1 ), h ), 2 ) - rips.t( h ) );
                end
                if max( diff( rips.edges, 1, 2 ) - diff( dt, 1, 2  ) ) > ceil( rips.Fs / min( rips.bp ) * 2.5 )
                    reDo    = 1;
                end
            end
        else
            reDo            = 1;
        end
        if reDo
            verb( sprintf( 'determining HFO times...' ), -1 );
            rips            = find_hfos( filebase, chan, 'bperiods', bperiods...
                , 'Overwrite', OverwriteSPWsearch, 'ripBP', ripBP, 'ripTH', ripTH...
                , 'mindurSECrip', mindurSECrip, 'minisiSECrip', minisiSECrip...
                , 'diffOrd', diffOrd, 'behaviorDurSEC', behaviorDurSEC...
                , 'padBuffer', padBuffer, 'filtMode', filtMode, 'powerMode', powerMode...
                , 'clipBase', clipBase );
        end
        if par.lfpSampleRate ~= rips.Fs
            vals            = round( Vals( :, 1 : 2 ) / par.SampleRate * rips.Fs );
        end
        if ~isempty( Vals )
            ridx            = isoverlap( rips.edges, vals, 0 );
            rips            = rips_select( rips, ~ridx );
        end
        if isempty( rips.trigs ) 
            continue
        end
        
        % statistics for each channel during the detected events:
        fprintf( '%d events detected. P2P estimation..', length( rips.trigs ) )

        % selection is based on (1) individual events (2) Gaussian filtering 
        trigtimes           = rips.trigs/rips.Fs;
        [ ~, ~, ~, xcsdhat, xlfphat ] = pt_avg( par, chans, trigtimes...
            , 'specChans', rips.chans( 1 )...
            , 'stawinSEC', staWin, 'tempSD', tempSD, 'spatBin', spatBin...
            , 'suffix', 'eeg'...
            , 'graphics', 0  );
        
        % for each triggering channel, select the channel that yielded
        % the highest amp events, based on the individual events:
        xi                  = permute( xlfphat, [ 2 3 1 ] ); 
        xi                  = reshape( xi, [ diff( winsamp ) + 1 length( rips.trigs ) * nchans ] );
        xhp                 = xi - firfilt( xi, gwin ); 
        mp2p                = nanmean( reshape( range( xhp( cidx, : ) ), [ length( rips.trigs ) nchans ] ), 1 );
        p2p{ i }( j, : )    = mp2p;
        [ ~, maxidx ]       = max( pr( : ).' .* p2p{ i }( j, : ) );

        % on the csd: here we are looking for a source of max amp, not of max amp of the bandpass
        if length( chans ) > 1 && ~isempty( xcsdhat ) % group is 1 channel only
            xi              = permute( xcsdhat, [ 2 3 1 ] );
            xi              = reshape( xi, [ diff( winsamp ) + 1 length( rips.trigs ) * nchans ] );
            switch selectionMode
                case 'csdPower'
                    xi      = reshape( xi, [ diff( winsamp ) + 1 length( rips.trigs ) * nchans ] );
                    xhp     = xi - firfilt( xi, gwin );
                    meanCSD = nanmean( reshape( range( xhp( cidx, : ) ), [ length( rips.trigs ) nchans ] ), 1 );
                otherwise % csdSource, lfpPower
                    meanCSD = nanmean( reshape( mean( xi( cidx, : ) ), [ length( rips.trigs ) nchans ] ) );
            end
            [ ~, maxidxCSD ] = max( pr( : ).' .* meanCSD );
            csdSelection    = [ chans( maxidxCSD ) meanCSD( maxidxCSD ) ];
        else
            csdSelection    = [ NaN NaN ];
        end
        mat( mat( :, 2 ) == i & mat( :, 1 ) == chan, 3 : 7 ) = [ chans( maxidx ) length( rips.trigs ) mp2p( maxidx ) csdSelection ];

    end
    
    % the simplest way is to select the channel with the highest LFP (or CSD) amplitude, 
    % regardless of the number of detections
    switch selectionMode
        case 'csdSource'
            cmat            = mat( mat( :, 2 ) == i, [ 6 4 7 ] ); % CSD
            cmat            = [ cmat( :, 1 ) 0.1 * cmat( :, 2 ) .* cmat( :, 3 ) / max( cmat( :, 2 ) ) + 0.9 * cmat( :, 3 ) ]; % weighted slightly by number of events
        case 'csdPower'
            cmat            = mat( mat( :, 2 ) == i, [ 6 4 7 ] ); % CSD
            cmat            = [ cmat( :, 1 ) 0.5 * cmat( :, 2 ) .* cmat( :, 3 ) / max( cmat( :, 2 ) ) + 0.5 * cmat( :, 3 ) ]; % weighted slightly by number of events
        case 'lfpPower'
            cmat            = mat( mat( :, 2 ) == i, [ 3 4 5 ] ); % LFP
            cmat            = [ cmat( :, 1 ) 0.1 * cmat( :, 2 ) .* cmat( :, 3 ) / max( cmat( :, 2 ) ) + 0.9 * cmat( :, 3 ) ]; % weighted slightly by number of events
    end
    [ ~, maxidx ]           = max( cmat( :, 2 ) );
    vote( i )               = cmat( maxidx, 1 );

    if vote( i ) == 0
        vote( i )           = NaN;
    end
    verb( sprintf( '%s: Group %d: selected channel: [ %d ]...', mfname, i, vote( i ) ), 1 );
    
end
verb( sprintf( '%s: Selected channels: [ %s ]...', mfname, num2str( vote( : ).' ) ), 1 );

%------------------------------------------------------------------------%
% second stage: recursion
%------------------------------------------------------------------------%
% now we have a set of voted channels, one per shank, that supposedly
% maximize the local ripple power. however that is not verified, this is
% done here:

% first get all group members (ignore no-detection channels)
nodetectChannels            = mat( mat( :, 4 ) <= 0, 1 );
chans                       = cell( ngroups, 1 );
for i                       = 1 : ngroups
    chan                    = vote( i );
    if isnan( chan )
        continue
    end
    % the other channels on the same shank:
    [ ~, ~, chans{ i } ]    = get_egroup( par, chan, 0 ); % use anatomical groups for detecting chan
    chans{ i }( ismember( chans{ i }, ignoredChannels ) ) = [];
    chans{ i }( ismember( chans{ i }, nodetectChannels ) ) = [];
end

% now go over groups:
vote0                       = vote;
trigs                       = cell( ngroups, 1 );
for i                       = 1 : ngroups
    chan                    = vote( i );
    if isnan( chan )
        continue
    end
    % the event times:
    optchans                = NaN * ones( size( chans{ i } ) );
    oldchan                 = NaN;
    j                       = 0;
    while ~isequal( oldchan, chan ) && j <= length( chans{ i } )
        j                   = j + 1;
        
        spwfname            = sprintf( '%s.spw.%s', filebase, num3str( chan ) );
        load( spwfname, '-mat' );
        if ~isempty( Vals )
            ridx            = isoverlap( rips.edges, vals );
            rips            = rips_select( rips, ~ridx );
        end
        trigs{ i }          = rips.trigs/rips.Fs;
        
        if isempty( trigs{ i } ) || length( chans{ i } ) < 3
            continue
        end
        oldchan = chan;
        % compute the max amp channel based on the mean LFP (not CSD!)
        [ avgcsd, avglfp ]  = pt_avg( par, chans{ i }, trigs{ i }...
            , 'whitenFlag', whitenFlag...
            , 'stawinSEC', staWin, 'tempSD', tempSD, 'spatBin', spatBin...
            , 'suffix', 'eeg'...
            , 'graphics', 0  );
        
        switch selectionMode
            case 'csdSource'            % max CSD source
                pp          = nanmean( avgcsd( cidx, : ) );
            case 'csdPower'             % max ripple-band CSD power
                xhp         = avgcsd - firfilt( avgcsd, gwin );
                pp          = range( xhp( cidx, : ) );
            case 'lfpPower'             % max ripple-band power
                xhp         = avglfp - firfilt( avglfp, gwin );
                pp          = range( xhp( cidx, : ) );
        end
        [ ~, maxidx ]       = max( pp );
        chan                = chans{ i }( maxidx );
        optchans( j )       = chan;
        % tried a few things here. the ripple-band power in the CSD is often 
        % much higher for artifacts/non-CA1pyr channels, so it makes more
        % sense to use the CSD mean instead (i.e. look for a source)
        % however, a physiological result that I see consistently is that 
        % the max CSD source is slightly above the max ripple power (in
        % both linear and dense probes), so since the objective here is to
        % detect HFOs (and not the CSD source itself), the
        % 'csdPower'/'lfpPower' mode is preferred
    end
    if j > length( chans{ i } )
        [ cnt, chn ]        = uhist( optchans );
        candidx             = cnt == max( cnt );
        cands               = chn( candidx );
        if sum( candidx ) == 1
            chan            = cands;
        elseif ismember( cands, vote( i ) )
            chan            = cands( ismember( cands, vote( i ) ) );
        elseif all( isnan( optchans ) )
            chan            = NaN; % no ripples
        else
            chan            = cands( 1 ); % arbitrary..
        end
    end
    if ~isequal( vote( i ), chan )
        fprintf( 'MODIFYING optimal channel: %d->%d (%d''th->%d''th)!\n'...
            , vote( i ), chan, find( chans{ i } == vote( i ) )...
            , find( chans{ i } == chan ) )
        vote( i ) = chan;
    end
end

verb( sprintf( '%s: Updated selection: [ %s ]...', mfname, num2str( vote( : ).' ) ), 1 );

%------------------------------------------------------------------------%
% compute stats + plot for the selected channels (CSD, spectra)
%------------------------------------------------------------------------%
paramnames                  = { 'egroup', 'best_channel', 'n', 'HFO/sec' ...
    , 'Freq (Hz)', 'Power (uV)', 'Power (SD)', 'Cycles/event', 'Duration (ms)' };
stats                       = zeros( ngroups, length( paramnames ) );
trigs                       = cell( ngroups, 1 );
chans                       = cell( ngroups, 1 );
for i                       = 1 : ngroups
    chan                    = vote( i );
    if isnan( chan )
        continue
    end
    spwfname                = sprintf( '%s.spw.%s', filebase, num3str( chan ) );
    if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( spwfname, 'file' ) )
        verb( sprintf( '%s: C%d: Computing HFO times...', mfname, chan ), 1 );
        rips                = find_hfos( filebase, chan, 'bperiods', bperiods...
            , 'Overwrite', Overwrite, 'ripBP', ripBP, 'ripTH', ripTH...
            , 'mindurSECrip', mindurSECrip, 'minisiSECrip', minisiSECrip...
            , 'diffOrd', diffOrd, 'behaviorDurSEC', behaviorDurSEC...
            , 'padBuffer', padBuffer, 'filtMode', filtMode, 'powerMode', powerMode...
            , 'clipBase', clipBase );
    else
        verb( sprintf( '%s: C%d: Loading HFO times...', mfname, chan ), 1 );
    end
    load( spwfname, '-mat' );
    if ~isempty( Vals )
        ridx                = isoverlap( rips.edges, vals );
        rips                = rips_select( rips, ~ridx );
    end
    rips.durs               = diff( rips.edges, [], 2 ) / rips.Fs * 1000; % ms
    if isempty( rips.durs )
        rips.ncycles        = [];
    else
        rips.ncycles        = rips.durs .* rips.f / 1000; % ncycles
    end
    orate                   = length( rips.t ) / rips.nsamps * rips.Fs; % occurrence rate
    stats( i, : )           = [ i chan length( rips.t ) orate mean( rips.f )...
        mean( rips.pow ) mean( rips.sd ) mean( rips.ncycles ) mean( rips.durs ) ];
    trigs{ i }              = rips.trigs/rips.Fs;
    [ ~, ~, chans{ i } ]    = get_egroup( par, rips.chans, 0 );
    chans{ i }( ismember( chans{ i }, ignoredChannels ) ) = [];
end

%------------------------------------------------------------------------%
% save
%------------------------------------------------------------------------%
% sps file: votes + n-events, mean frequency, mean amp (SD/power)
if Overwrite >= 1 || ( Overwrite ~= -1 && ~exist( spsfname, 'file' ) )
    verb( sprintf( '%s: saving %s...', mfname, spsfname ), vflag )
    save( spsfname, 'filebase', 'vote', 'stats', 'vote0', 'mat', '-v6' )
end

%------------------------------------------------------------------------%
% plot
%------------------------------------------------------------------------%
if graphics
    if Overwrite == 1
        recomputeBL         = 1;
        baselineEpochs      = bperiods;
    else
        recomputeBL         = 0;
        baselineEpochs      = [];
    end
    for i                   = 1 : ngroups
        if isnan( vote( i ) )
            continue
        end
        pt_avg( par, chans{ i }, trigs{ i }, 'specChans', vote( i )...
            , 'stawinSEC', staWin, 'tempSD', tempSD, 'spatBin', spatBin...
            , 'suffix', 'wlt', 'baselineEpochs', baselineEpochs ...
            , 'whitenFlag', whitenFlag, 'normalizeSWS', normalizeSWS, 'recomputeBL', recomputeBL...
            , 'graphics', figdir, 'savetype', savetype );
        if gcf >= 10 || BatchMode
            close all
        end
    end
end

return

% EOF
