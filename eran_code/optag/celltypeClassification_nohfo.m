% celltypeClassification        by CC, optics, and waveform/temporal stats
%
% call              s = celltypeClassification( filebase )
%
% optional arguments (name/value pairs):
% 
%                   ilevel      {'B'}; isolation level
%                   suffix      {'s2s.nohfo'}; spikes2spikes output on which to run check_mono 
%                   supFlag     { 0 }; by default, tag optically units with rate increase
%                                   if supFlag is 1, optical tagging according to rate decrease
%                                   this argument is also transferred to selectIntensity
%                   celltype    { 1 }; by default, selectIntensity is
%                                   called with celltype 1 (PYR); 0 corresponds to INT
%                   opsinType   {'ChR2'}; only for text purposes
%
%                   durRange    {[0.04 0.08]}; [s], duration of stimuli to consider in DCanalysis
%                   wavRange    {[400 500]}; [nm], wavelengths to consider 
%                   stimType    { 'PULSE', 'PSINE' }
%                   sourceType  { 'LED' }
%
%                   uflag       { 1 }; arguments for get_triggers
%                   multi       { 'eq' }; arguments for get_triggers
%                   simOnly     { 0 }; arguments for get_triggers
%
%                   slevel      { [] }; [A], if not empty, forces a range of current level
%                   jitflag     { 0 }; argument for classify_waveform_plot (set to 1 if many units)
%                   nbins       { 100 }; argument for opticalTagging (PSTH)
%                   pTH         { 0.01 }; local and argument for selectIntensity
%                   minSpikes   { 10 }; local and argument for selectIntensity
%
%                   graphics    { 1 }
%                   Overwrite   { -2 }
%                   savetype    { 'png' }
%
% returns           s           structure with waveform, connectivity, and
%                                   optical information for each unit. if
%                                   multiple stimulation channels were
%                                   used, some of s fields will be 3D
%
% does
%                   (1) gets waveform features (spikes_stats, *sst)
%                   (2) gets connectivity data (check_mono, *s2s)
%                   (3) gets optical stimulution data (opticalTagging)
%                   (4) organizes in two structures, saves, and plots
%
% calls
%           LoadXml                                             (blab)
%           get_stimchans, get_triggers, LoadStims              (formats)
%           ParseArgPairs                                       (general)
%           fig_out, replacetok, textf                          (graph)
%           DCanalysis, opticalTagging, selectIntensity         (optag)
%           check_mono, classify_waveform, classify_waveform_plot, 
%               determine_units, load_spikes, 
%               spikes_stats_depth, wfeatures                   (spikes)
%           struct_select                                       (structs)
% 
% notes about using this routine
%
% (1) this is a high-level routine, intended to be used at the level of 
% the recording session, isolation level, and stimulation channel. 
%
% For enhanced user control, there are also options to control stimulation
% parameters such as durRange, uflag/multi/simOnly, and slevel. However,
% changing these parameters will NOT modify the name of the output file or
% figure, but will overwrite them. In contrast, the detailed information is
% generated by DCanalysis and saved in *mat files with unique names.
%
% (2) this routine does not actually classify
% To actually classify, the routine classify_waveform.m is called at the
% end of this routine, generating the figure *celltypeClassification
% 
% If you want to create a new classify_waveform, then:
% 	1. use this routine (celltypeClassification) to process the data from all relevant sessions
%   2. then make a classifier from the tagged data (train data)
%   3. then classify each individual session (test data)
%
% see also          classify_waveform (actually classifies)
%                   make_gm_classifier (makes a classifier)
%
% assumptions       check_mono pruned the data

% 16-may-13 ES

% revisions
% 14-dec-19 cleaned up
% 15-dec-19 added celltype argument
% 16-dec-19 (1) changed call to opticalTagging
%           (2) changed figure name
% 17-dec-19 (1) distinct channels of the same type (e.g. four blue LEDs)
%               were previously pooled into the same set of stimuli,
%               however this may dilute effects. changed to analyze all
%               clusters for each channel separately, and Bonferroni
%               correct the results. thus some of the fields in s may be 3D
%           (2) nbins and pTH added as input arguments
% 19-dec-19 (1) if no stimuli in range, fill the relevant parts of s with empty fields
%           (2) optical tagging done for all units, for each stim channel separately
% 23-dec-19 (1) opticalTagging called with Overwrite, since saving is not ambivalent any more

% what this routine does is accumulate the statistics for each unit, namely
% 1. basic waveform stats (trough-peak (d); width (t); half-width (hw); asymetry index (asy))
% 2. short-time cross-correlation stats (excitatory, inhibitory, excited, inhibited)
% 3. optical properties (p-value for being activated/suppressed by light; optical tagging flag (act))
%       also keeps the 201-bin PETH for each cell plus some additional
%       stats (pinfo: shank/clu/opt:1act,-1sup,NaN/ntrials/mean_duration[samples])
%             (pstats: n_out/n_in/r_out/r_in/p_act/p_sup/latency[ms],activated only)
% 4. elementary spiking statistics (rate, ACH center-of-mass[ms] over the first 50 ms)
%
% 
% to do:
% (3) go over all *DC.PULSE.*u1* files in the path . then take the decision for
% each local cell. If not all cells were examined under uniuqe conditions,
% (i.e. there is no *u1* files, or the existing do not apply to all cells),
% check if there is a *u-1* file; if there is, tag the unexamined cells
% based on that

function [ s, stats ] = celltypeClassification( filebase, varargin )

%------------------------------------------------------------------
% i/o
%------------------------------------------------------------------

% initialize
stats                       = [];

% constants
fTrials                     = 0.1;              % should be a spike on average on at least 10% of the trials

% arguments
nargs                       = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ ilevel, suffix, opsinType, celltype, supFlag ...
    , durRange, wavRange, stimType, sourceType ...
    , uflag, multi, simOnly ...
    , slevel, jitflag, nbins, pTH, minSpikes ...
    , graphics, Overwrite, savetype ] = ParseArgPairs(...
    { 'ilevel', 'suffix', 'opsinType', 'celltype', 'supFlag' ...
    , 'durRange', 'wavRange', 'stimType', 'sourceType' ...
    , 'uflag', 'multi', 'simOnly' ...
    , 'slevel', 'jitflag', 'nbins', 'pTH', 'minSpikes' ...
    , 'graphics' ,'Overwrite', 'savetype' }...
    , { 'B', 's2s.nohfo', [], 1, 0 ...
    , [ 0.04 0.08 ], [ 400 500 ], { 'PULSE', 'PSINE' }, 'LED' ...
    , 1 'eq', 0 ...
    , [], 0, 100, 0.01, 10 ...
    , 1, -2, 'png' }...
    , varargin{ : } );

%------------------------------------------------------------------
% set up - general
%------------------------------------------------------------------

% paths
delim                       = strfind( filebase, '/dat/' );
if isempty( delim )
    fprintf( '%s: Cannot save fig and/or data\n', upper( mfilename ) )
end
if isa( graphics, 'char' ) && exist( graphics, 'dir' )
    figdir                  = graphics;
    graphics                = 1;
else
    figdir                  = [ filebase( 1 : delim ) 'figs/dc' ];
    if ~exist( fileparts( figdir ), 'dir' )
        mkdir( fileparts( fileparts( figdir ) ), 'figs' )
    end
    if ~exist( figdir, 'dir' )
        mkdir( fileparts( figdir ), 'dc' )
    end
end
matdir                      = [ filebase( 1 : delim ) 'mat/dc' ];
if ~exist( fileparts( matdir ), 'dir' )
    mkdir( fileparts( fileparts( matdir ) ), 'mat' )
end
if ~exist( matdir, 'dir' )
    mkdir( fileparts( matdir ), 'dc' )
end
[ ~, filename, extname ]    = fileparts( filebase );
filename                    = [ filename extname ];

% determine the stimulation channels
[ stimchans, ~, ~, stimsources, stimwavelengths ] = get_stimchans( filebase );
kidx                        = inrange( stimwavelengths, wavRange );
kidx                        = kidx & ismember( stimsources, sourceType );
stimchans( ~kidx )          = [];
stimsources( ~kidx )        = [];
if isempty( stimchans )
    fprintf( '%s: No relevant stim chans in %s, aborting.\n', filebase, upper( mfilename ) )
    return
end
tstr                        = ( [ repmat( 'T', size( stimchans( : ) ) ) num2str( stimchans( : ) ) ]' );
tstr                        = tstr( : ).';

% define the corename
corename                    = sprintf( '%s.%s_c%s.celltypeClassification', filename, tstr, ilevel );
figname                     = [ figdir '/' corename ];
savename                    = [ matdir '/' corename ];

if Overwrite < 0 && exist( savename, 'file' )
    fprintf( '%s: loading %s...\n', upper( mfilename ), savename )
    load( savename, 's', 'stats', '-mat' );
    return
end

% necessary files
par                         = LoadXml( filebase );
spkFs                       = par.SampleRate;
opsinType                   = lower( opsinType );

%--------------------------------------------------------------------%
% determine which units
%--------------------------------------------------------------------%
shankclu                    = determine_units( filebase, [], ilevel );
nclu                        = size( shankclu, 1 );
pidx                        = shankclu( :, 3 );


%--------------------------------------------------------------------%
% waveform features
%--------------------------------------------------------------------%
fprintf( '%s: Getting waveform features (requires *sst file)...\n', upper( mfilename ) )
if ~exist( [ filebase '.sst' ], 'file' )
    msg                     = sprintf( '%s: Missing %s.sst file!!', upper( mfilename ), filebase );
    throw( MException( 'hfoAnalysis:MissingSSTfile', msg ) );
end
sstfname                    = [ filebase '.sst' ];
try % added 11-dec-19
    sst                     = spikes_stats_depth( filebase );
catch
    load( sstfname, 'sst', '-mat' )
    sst.depth               = NaN * ones( size( sst.shankclu, 1 ), 1 );
end
shankclu0                   = sst.shankclu;
gidx                        = ismember( shankclu0( :, 1 : 2 ), shankclu( :, 1 : 2 ), 'rows' );
if sum( gidx ) ~= nclu
    error( 'mismatch!!' )
end
[ sst, rc, fieldnames ]     = struct_select( sst, gidx ); % dilute to get same units
if sum( rc == 2 ) > 0
    fieldname               = fieldnames( rc == 2 );
    if isequal( fieldname{ 1 }, 'frateb' )
        sst.frateb          = sst.frateb( :, gidx );
    end
    if isequal( fieldname{ 1 }, 'max' )
        sst.max( :, ~gidx ) = [];
    end
    fieldname{ 1 }          = [];
    if ~isempty( fieldname{ 1 } ) && ( length( fieldname ) > 1 )%|| ~isequal( fieldname{ 1 }, 'max' ) )
        error( 'mismatch' )
    end
end
t                           = 1 ./ sst.fmax * 1000;         % the total spike width
acom                        = sst.ach_com;                  % [ms]
[ hw1, asy1, tp1 ]          = wfeatures( sst.max ) ;
if ~isequal( length( hw1 ), length( asy1 ), length( tp1 ), length( t ) )
    fprintf( '%s: CHECK wfeatures output!!!\n', upper( mfilename ) )
    keyboard
end
d                           = tp1( : ) / spkFs * 1000;      % [ms]
hw                          = hw1( : ) / spkFs * 1000;      % [ms]
asy                         = asy1( : );

%--------------------------------------------------------------------%
% the functional connectivity data:
%--------------------------------------------------------------------%
fprintf( '%s: Checking mono-synaptic connectivity...(*%s file)...\n', upper( mfilename ), suffix )
mono                        = check_mono_nohfo( filebase, gidx, 'suffix', 's2s.nohfo' );   % use only those units
if ~isempty( mono.both )                                                        % WTA for special case of dual inibitory and excitatory CCH
    eidx                    = mono.both( :, 2 ) > mono.both( :, 3 );
    mono.exc                = [ mono.exc; mono.both( eidx, 1 ) ];
    iidx                    = mono.both( :, 2 ) < mono.both( :, 3 );
    mono.inh                = [ mono.inh; mono.both( iidx, 1 ) ];
end
if ~isequal( mono.shankclu, shankclu( :, 1 : 2  ) )
    error( 'mismatch!!1' )
end

%--------------------------------------------------------------------%
% the light activation data:
%--------------------------------------------------------------------%
fprintf( '%s: Checking optical data (*par, *stm, *val*)...\n', upper( mfilename ) )

vals                        = LoadStims( filebase );

% base activation on single-shank PULSE stimulation.
didx                        = false( nclu, 1 ); % decided upon (tagged or not)
tidx                        = false( nclu, 1 ); % tested cell

% if slevel not enforced externally, run selectIntensity for each channel separately
nchans                      = length( stimchans );
nlevels                     = size( slevel, 1 );
if nchans == nlevels
elseif nlevels == 1
    slevel                  = ones( nchans, 1 ) * slevel;            
else
    slevel                  = NaN * ones( nchans, 2 );
end

% prepare blackout and spikes
blackout                    = vals( :, 1 : 2 );
fprintf( '%s: Loading spikes...\n', upper( mfilename ) )
spk                         = load_spikes( filebase, [], ilevel );

% optical tagging for each channel separately
allPeriods                  = cell( nchans, 1 );
allCorenames                = cell( nchans, 1 );
nperiods                    = zeros( nchans, 1 );
pstats                      = NaN * ones( nclu * nchans, 7 );
pinfo                       = NaN * ones( nclu * nchans, 6 );
peth                        = NaN * ones( nclu, 2 * nbins + 1, nchans );
bins                        = NaN * ones( 1, 2 * nbins + 1, nchans );
for i = 1 : nchans
    
    achan                   = stimchans( i );
    
    % check if slevel was provided, if not, run selectIntensity
    if any( isnan( slevel( i, : ) ) )
        fprintf( '%s: Selecting intensity...\n', upper( mfilename ) )
        slevel( i, : )      = selectIntensity( filebase, 'ilevel', ilevel...
            , 'stimchans', achan, 'wavRange', wavRange, 'sourcetypes', stimsources ...
            , 'durRange', durRange ...
            , 'celltype', celltype, 'supFlag', supFlag  ...
            , 'minSpikes', minSpikes, 'pTH', pTH ...
            , 'Overwrite', -2 );
    end
    
    % run get_triggers for that level
    params                  = { 'types', stimType, 'durs', durRange, 'vals', slevel( i, : ) };
    [ ~, tims, durs ]       = get_triggers( filebase, achan, uflag, multi, simOnly, params );
    periods                 = [ tims tims + round( durs * spkFs ) - 1 ];
    allPeriods{ i }         = periods;
    nperiods( i )           = size( periods, 1 );
    allCorenames{ i }       = sprintf( '%s.%s_c%s.celltypeClassification' ...
        , filename, sprintf( 'T%d', achan ), ilevel );

    % check if any periods, if not empty, run opticalTagging
    if nperiods( i ) <= 10
        continue
    end
    periodsI            = allPeriods{ i };
    corenameI           = allCorenames{ i };
    [ pstatsI, pinfoI, pethI, binsI ] = opticalTagging( spk.clu, spk.res, spk.map, periodsI...
        , 'blackout', blackout, 'ispyr', spk.shankclu( :, 3 )...
        , 'spkFs', spkFs, 'nbins', nbins, 'corename', corenameI  ...
        , 'filebase', filebase, 'Overwrite', Overwrite, 'graphics', graphics );
    idx                 = ( 1 : nclu ) + nclu * ( i - 1 );
    pstats( idx, : )    = pstatsI;
    pinfo( idx, : )     = [ pinfoI ones( nclu, 1 ) * achan ];
    peth( :, :, i )     = pethI';
    bins( 1, :, i )     = binsI;
    
end

% detemine whether light tagged:
pact                        = pstats( :, 5 );
psup                        = pstats( :, 6 );
gain                        = calc_gain( pstats( :, [ 4 3 ] ) );

kidx                        = nperiods > 0;
pTHcorr                     = pTH / sum( kidx );
if supFlag
    sidx                    = psup <= pTHcorr;                               % suppressed by Poisson
    slct                    = sidx & gain < 1;
else
    sidx                    = pact <= pTHcorr;
    ridx                    = pstats( :, 2 ) ./ pinfo( :, 4 ) > fTrials;       % enough spiking trials
    cidx                    = pstats( :, 2 ) > minSpikes;                      % enough spikes
    slct                    = sidx & gain > 1 & cidx & ridx;
end

% check if all units were tested/tagged
tested                  = ismember( spk.shankclu( :, 1 : 2 ), pinfo( :, 1 : 2 ), 'rows' );
tagged                  = ismember( spk.shankclu( :, 1 : 2 ), pinfo( slct, 1 : 2 ), 'rows' );
didx( tagged )          = 1;
tidx( tested )          = 1;

if sum( ~tidx ) > 0
    fprintf( 1, '\n\n\n\nNot all units identified by determine_units.m were tested by opticalTagging!!\n\n\n\n' )
end

%--------------------------------------------------------------------%
% summarize:
%--------------------------------------------------------------------%

inh                         = false( nclu, 1 );
exc                         = false( nclu, 1 );
inhibited                   = false( nclu, 1 );
excited                     = false( nclu, 1 );

inh( mono.inh )             = 1;
exc( mono.exc )             = 1;
inhibited( mono.inhibited ) = 1;
excited( mono.excited )     = 1;

% session-level summary
stats.filename              = { filename };
stats.ncchs                 = nclu * ( nclu - 1 ) / 2;
stats.ncchExc               = size( mono.pairsExc, 1 );
stats.ncchInh               = size( mono.pairsInh, 1 );

% unit-level summary
s.filename                  = repmat( { filename }, [ nclu, 1 ] );
s.opsinType                 = repmat( { opsinType }, [ nclu, 1 ] );
s.shankclu                  = shankclu;
if isempty( slevel )
    slevel                  = NaN;
end
s.slevel                    = repmat( slevel, [ nclu 1 ] );

% waveform data
s.t                         = t;                        % width [ms]
s.d                         = d;                        % trough-peak time [ms]
s.acom                      = acom;                     % ACH COM [ms]
s.hw                        = hw;                       % half-width
s.asy                       = asy;                      % assymetry
s.pidx                      = pidx;                     % PYR or INT 

% connectivity data (from check_mono)
s.exc                       = exc;                      % excitatory
s.inh                       = inh;                      % inhibitory
s.excited                   = excited;                  % by others
s.inhibited                 = inhibited;                % by others

% from optical tagging
s.pact                      = reshape( pact, [ nclu nchans ] );
s.psup                      = reshape( psup, [ nclu nchans ] );
s.gain                      = reshape( gain, [ nclu nchans ] );
s.act                       = didx;                     % tagged by light

% waveform data
s.depth                     = sst.depth;

% from optical tagging
s.peth                      = peth;
s.bins                      = repmat( bins, [ nclu 1 ] );
s.pinfo                     = permute( reshape( pinfo', [ 6 nclu nchans ] ), [ 2 1 3 ] );
s.pstats                    = permute( reshape( pstats', [ 7 nclu nchans ] ), [ 2 1 3 ] );

if length( s.pidx ) ~= length( s.pact )
    error( 'Not all units identified by determine_units.m were tested by opticalTagging' )
end
    
%--------------------------------------------------------------------%
s.pinfo                     = permute( reshape( pinfo', [ 6 nclu nchans ] ), [ 2 1 3 ] );
s.pstats                    = permute( reshape( pstats', [ 7 nclu nchans ] ), [ 2 1 3 ] );

if length( s.pidx ) ~= length( s.pact )
    error( 'Not all units identified by determine_units.m were tested by opticalTagging' )
end
    
%--------------------------------------------------------------------%
% save
%--------------------------------------------------------------------%
if ( Overwrite == 1 || ~isempty( savename ) && ~exist( savename, 'file' ) && Overwrite ~= -1 )
    fprintf( '%s: saving %s..\n', upper( mfilename ), savename )
    save( savename, 's', 'stats', '-v6' );
end

%--------------------------------------------------------------------%
% plot
%--------------------------------------------------------------------%
if ~graphics
    return
end

fig                         = figure;
subplot( 2, 2, 1 )
classify_waveform_plot( s, [], jitflag, 1, 0 );
[ ~, ~, fsep ]              = classify_waveform( [ s.d s.t ], 0 );
fh                          = fimplicit( fsep, [ xlim ylim ] );
set( fh, 'color', [ 0 0 0 ] );
title( '' )
ylabel( 'Duration [ms]' )
xlabel( 'Trough-to-peak [ms]' )

subplot( 2, 2, 2 )
shat                        = s;
shat.t                      = shat.acom;
classify_waveform_plot( shat, [], jitflag, 1, 0 );
ylims                       = ylim;
ylim( [ ylims( 1 ) 35 ] )
ylabel( 'ACH COM [ms]' )

textf( 0.5, 0.975, replacetok( corename, '\_', '_' ) );

if ~isempty( savetype ) && ( isa( savetype, 'cell' ) || ~all( isnan( savetype ) ) ) && ~isempty( figname )
    if ~isa( savetype, 'cell' )
        savetype            = { savetype };
    end
    for j                   = 1 : length( savetype )
        fignameI        = [ figname '.' savetype{ j } ];
        fprintf( '%s: Saving figure %s\n', upper( mfilename ), fignameI )
        fig_out( fig, 1, fignameI, savetype{ j } );
    end
end
    
return

% EOF

% blue LED, ChR2
s = celltypeClassification( filebase, 'sourceType', 'LED', 'Overwrite', 1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype', 0, 'supFlag', 0 );
% blue LD, ChR2
s = celltypeClassification( filebase, 'sourceType', 'LD', 'Overwrite', 1, 'wavRange', [ 400 500 ], 'durRange', [ 0.04 0.08 ], 'celltype', 0, 'supFlag', 0 );
% red LD, Jaws
s = celltypeClassification( filebase, 'sourceType', 'LD', 'Overwrite', 1, 'wavRange', [ 550 700 ], 'durRange', [ 0.15 0.25 ], 'celltype', 0, 'supFlag', 1 );
