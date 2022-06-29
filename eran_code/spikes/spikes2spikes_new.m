% spikes2spikes             compute s2s structure with ASG and STC information
%
% call                      s2s = spikes2spikes( filebase )
%
% gets                      filebase
% 
% returns                   s2s structure
%
% optional arguments (given as name/value pairs):
%
%               CCH computation parameters:
%
%                           isDeadTime          {1}         assume deadtime in same-shank detection
%                           deadTimeMS          {0.4}       [ms], ignored if isDeadTime is 0
%                           BinSizeMS           {1}         [ms], CCH bin size 
%                           halfWidthMS         {50}        [ms], CCH half-width 
%
%                           jitWindowMS         {5}         [ms], filtering half-width
%                           roiMS               {[NaN 5]}   [ms], for significance testing 
%                           convType            {'gauss'}   argument to cch_conv
%                           alfa                {0.001}     significance testing (global bands)
% 
%                           supportEdges        {[-1 3]}    jitWindowMS multiples i.e. keep
%                                                               h(t) only for the bins -5 : 15 
%                           asgMode             {1}         0   detects global maxima in ROI (always exists, may be non-causal)
%                                                           1   detects local maxima in ROI (may not exist, but preserves causality)
%                           roiHardness         {[1 0]}     hard/soft edges (limits/does not limit ASG base by ROI)
%
%               Flow control:
%                           verbose             {1}         
%                           suffix              {'s2s'}     of file on disk
%                           Overwrite           {-2}        1: compute and overwrite 
%                                                           0: only compute
%                                                           -1: load/compute but do not write
%                                                           -2: load/compute and write
%               Data selection parameters:
%
%                           shanknums           {[]}
%                           periods             {-1}        2-column matrix     periods, in samples, to keep spikes from
%                                                           -1                  flag, remove spikes during stims (segments specified in val files)
%                                                           -2                  flag, remove spikes during HFOs
%                                                           -3                  flag, remove spikes during stims and during HFOs 
%
% calls                     CCG, LoadXml                                        (blab)
%                           ParseArgPairs                                       (general)
%                           calc_stc_asg, cch_conv, cch_deconv, load_spikes     (spikes)
%                           LoadStims                                           (formats)
%                           get_hfo_times                                       (lfp) 
%                           inranges, resampleranges, sortranges, uniteranges   (sets)
%
% Note:
%               Short epochs, interval jitter, and other finesse not supported in this routine
%               For such analyses, call calc_asg with a single pair at a time
%
% see also                  calc_asg

% 03-mar-21 ES

% revisions
% 05-mar-21 matrix version (~100 times faster)
% 06-mar-21 added dead time support
% 07-mar-21 added ASG + STC computations
% 16-mar-21 (1) bug fix for nidx: was t_ROI, changed to t_ROID
%           (2) bug fix for nidx: gbUpper0/gbLower0 assigned 
%           (3) special treatment for a-causal bin (t==0)
%           (4) added new fields gbUpperZ, gbLowerZ
% 24-mar-21 cleaned up and renamed as spikes2spikes

% to do:
% (1) run for a few sessions, then rename this routine 'spikes2spikes.m'

function s2s = spikes2spikes( filebase, varargin )

% constants
expansionFactor                 = 4;                                      	% expansion factor for dual deconvolution
dcMethod                        = 'fft';

% argument handling
nargs                           = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ isDeadTime ...
    , deadTimeMS, BinSizeMS, halfWidthMS, jitWindowMS, roiMS ...
    , convType, alfa ...
    , supportEdges, asgMode, roiHardness ...
    , shanknums, ilevel, periods ...
    , verbose, suffix, Overwrite ]            = ParseArgPairs( ...
    { 'isDeadTime' ...
    , 'deadTimeMS', 'BinSizeMS', 'halfWidthMS', 'jitWindowMS', 'roiMS' ...
    , 'convType', 'alfa' ...
    , 'supportEdges', 'asgMode', 'roiHardness' ...
    , 'shanknums', 'ilevel', 'periods' ...
    , 'verbose', 'suffix', 'Overwrite' } ...
    , { 1 ...
    , 0.4, 1, 50, 5, [ NaN 5 ] ...
    , 'median', 0.001 ...
    , [ -1 3 ], 1, [ 1 0 ] ...
    , [], 'E', -1 ...
    , 1, '', -2 }...
    , varargin{ : } );

% i/o
if isempty( suffix )
    filename                 	= [ filebase '.s2s' ];
else
    filename                	= [ filebase '.' suffix ];
end
if verbose
    mfname                      = upper( mfilename );
end
if exist( filename, 'file' ) && Overwrite < 0 && nargout > 0
    if verbose
        fprintf( 1, 'Loading %s\n', filename )
    end
    load( filename, '-mat', 's2s' );
    return
end
if ~exist( fileparts( filebase ), 'dir' )
    fprintf( 1, '%s: missing path %s, exiting\n', mfname, fileparts( filebase ) )
    s2s                         = [];
    return
end

%------------------------------------------------------------------------
% (1) load spikes

if verbose
    t0                          = clock;
    fprintf( 1, '%s: Loading spikes and selecting periods... ', mfname )
end

% load the clu/res/map
spk                             = load_spikes( filebase, shanknums, ilevel );  % use
par                             = LoadXml( filebase );
SpikesFs                        = par.SampleRate;

% keep only the relevant spikes
if periods( 1 ) < 0
    
    % remove spikes during stimulation
    [ Vals, Trigs ]             = LoadStims( filebase );
    if periods( 1 ) == -1 || periods( 1 ) == -3                             % any trigger
        tidx                    = true( size( Trigs ) );
    elseif periods( 2 ) == -2                                               % specific subset
        tidx                    = ismember( Trigs, abs( periods ) );
    else
        tidx                    = [];
    end
    if isempty( tidx )
        uvals1                  = [];
    else
        uvals1                  = uniteranges( Vals( tidx, 1 : 2 ) );       % combine all the segments
    end
    if periods( 1 ) == -2 || periods( 1 ) == -3                             % remove HFO times
        hperiods                = get_hfo_times( filebase );
        hperiods                = sortranges( hperiods );
        eegFs                   = par.lfpSampleRate;
        uvals2                  = resampleranges( hperiods, SpikesFs, eegFs );
    else
        uvals2                  = [];
    end
    uvals                       = uniteranges( uvals1, uvals2 );
    if ~isempty( uvals )
        ridx                    = inranges( spk.res, uvals );               % remove any spike that is in any segment
        spk.res( ridx )         = [];
        spk.clu( ridx )         = [];
    end
    Bsec                        = sum( diff( uvals, [], 2 ) + 1 ) / SpikesFs;
    Rsec                        = spk.res( end ) / SpikesFs;
    Tsec                        = Rsec - Bsec;

elseif ~isempty( periods ) && size( periods, 2 ) == 2
    
    % keep only spikes during specified periods
    kidx                        = inranges( spk.res, periods );
    spk.res                     = spk.res( kidx );
    spk.clu                     = spk.clu( kidx );
    Tsec                        = sum( diff( periods, [], 2 ) + 1 ) / SpikesFs;
    
end

uClu                            = unique( spk.map( :, 1 ) );                % in case a unit only spikes during masked epochs
if ~isempty( setdiff( spk.map( :, 1 ), uClu ) )
    ridx                        = ~ismember( spk.map( :, 1 ), uClu );
    spk.shankclu( ridx, : )     = [];
    spk.map( ridx, : )          = [];
end
n                               = length( uClu );

if verbose
    et                          = etime( clock, t0 );
    fprintf( 1, ' done (%0.2g s)\n', et );
end

%------------------------------------------------------------------------
% (2) compute ACHs and CCHs

if verbose
    t0                          = clock;
    fprintf( 1, '%s: Computing count CCHs... ', mfname )
end

% prepare parameters for CCH time lag vector
nBins                           = halfWidthMS / BinSizeMS;               	% [bins] number of bins on each side of the CCH
t                               = ( -nBins : nBins )' * BinSizeMS;      	% [ms]
dt                              = diff( t( 1 : 2 ) ) / 1000;              	% [s]
m                               = length( t );                              % [bins] total number of CCH bins
mef                             = 2 * expansionFactor * nBins + 1;          % [bins] length of expanded CCH
BinSize                         = round( BinSizeMS * SpikesFs / 1000 );  	% [samples]

% prepare parameters for STC support
idx                             = ( supportEdges( 1 ) * jitWindowMS ) : ( supportEdges( 2 ) * jitWindowMS );
sidx                            = find( ismember( t, idx ) );
nsupport                        = length( sidx );                           % [bins]

% counts spikes for each unit
nspks0                          = hist( spk.clu, uClu )';                   % spike count

% compute raw CCH and ACH for all pairs and units
ccg                             = CCG( spk.res, spk.clu, BinSize, expansionFactor * nBins, SpikesFs, uClu, 'count' ); 
% ccg is organized s.t. ccg( :, n1, n2 ) is the CCH from n1->n2

if verbose
    et                          = etime( clock, t0 );
    fprintf( 1, ' done (%0.2g s)\n', et );
end

%------------------------------------------------------------------------
% (3) deconvolve ACHs from CCH

if verbose
    t0                          = clock;
    fprintf( 1, '%s: Deconvolving ACHs from CCHs... ', mfname )
end

% reshape cch and extract ach to prepare for deconvolution
cch0                            = reshape( ccg, [ mef n * n ] );            % columns are [ 11 12 ... 1n 21 22 ... 2n ... ]

% source column indices for ACH (n x n matrix organized as an 1 x n^2 vector): 
aidx                            = ( 0 : ( n - 1 ) ) * n + ( 1 : n );
aidx                            = ones( n, 1 ) * aidx;
aidx                            = aidx( : );

% target column indices for ACH1:
aidx1                           = reshape( ( 1 : n^2 ), [ n n ] )';
aidx1                           = aidx1( : );
ach1                            = zeros( mef, n * n );
ach1( :, aidx1 )                = cch0( :, aidx );

% target column indices for ACH2:
aidx2                           = ( 1 : n^2 )';
ach2                            = zeros( mef, n * n );
ach2( :, aidx2 )                = cch0( :, aidx );

% spike count indices
nidx                            = ( 1 : n )' * ones( 1, n );
nidx                            = nidx( : );
nspks1                          = zeros( 1, n * n );
nspks2                          = zeros( 1, n * n );
nspks1( aidx2 )                 = nspks0( nidx );
nspks2( aidx1 )                 = nspks0( nidx );

% actually run the deconvolution
[ cch, ~, kidx ]                = cch_deconv( cch0, ach1, nspks1, ach2, nspks2, dcMethod, mef );

% cch matrix includes deconvolved ACHs, replace those with the raw ones
aidx                            = ( 0 : ( n - 1 ) ) * n + ( 1 : n );
cch( :, aidx )                  = cch0( kidx, aidx );

% keep only central part of the raw count cch
cch0                            = cch0( kidx, : );
sem                             = sqrt( cch0 );

if verbose
    et                          = etime( clock, t0 );
    fprintf( 1, ' done (%0.2g s)\n', et );
end

%------------------------------------------------------------------------
% (4) compute predictor and p-values

if verbose
    t0                          = clock;
    fprintf( 1, '%s: Computing predictors and p-values... ', mfname )
end

% dead time considerations
shank                           = spk.map( :, 2 );
idx0                            = shank( ( 1 : n )' * ones( 1, n ) ) == shank( ones( n, 1 ) * ( 1 : n ) );
idx0                            = reshape( idx0, [ 1 n * n ] );             % same shank (including same unit)
if isDeadTime
    nidx                        = t <= deadTimeMS & t >= -deadTimeMS;       % exclude these bins from all subsequent steps
    idx1                        = t < -deadTimeMS;
    idx2                        = sum( idx1 ) + 1 : sum( ~nidx );           % indices of the excluded vector
    idx1                        = 1 : sum( idx1 );
else
    nidx                        = false( length( t ), 1 );
end

% convolution window
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % optional: causality imposed
end
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );        % causality not necessarily imposed
nBonf                           = sum( t_ROI );
W                               = 2 * ceil( jitWindowMS / BinSizeMS ) + 1;  % [samples]
tZ                              = t == 0;                                   % give the zero-lag bin special treatment

% same shank pairs - compute without the dead-time bins, then replace with NaNs
if sum( nidx )
    
    t_ROID                      = t_ROI & ~nidx;                            % dead time imposed
    nBonfD                      = sum( t_ROID );
    
    cchD                        = cch( :, idx0 );                           % same-shank pairs (incl. ACH)
    cchD                        = cchD( ~nidx, : );                         % ignore under-estimated bins d.t. dead time
    
    % significance (point-wise probabilities; predictor)
    [ pvalsUpper, pred, pvalsLower ]    = cch_conv( cchD, W, convType ); 	% linear/non-linear filtering of the CCH
    
    % global limits, with a conservative Bonferroni correction, on the tested range:
    gbUpper                     = poissinv( 1 - alfa / nBonfD, max( pred( t_ROID, : ), [], 1 ) );
    gbLower                     = poissinv( alfa / nBonfD, min( pred( t_ROID, : ), [], 1 ) );
    gbUpperZ                    = poissinv( 1 - alfa, max( pred( tZ, : ), [], 1 ) );
    gbLowerZ                    = poissinv( alfa, min( pred( tZ, : ), [], 1 ) );
    
    % detect the bins with high/low global counts
    hiBins                      = false( size( cchD ) );
    loBins                      = false( size( cchD ) );
    hiBins( t_ROID, : )         = bsxfun( @gt, cchD( t_ROID, : ), gbUpper ) & ( cchD( t_ROID, : ) > 0 );
    loBins( t_ROID, : )         = bsxfun( @lt, cchD( t_ROID, : ), gbLower ) & ( ones( nBonfD, 1 ) * gbLower ) > 0;
    hiBins( tZ, : )             = bsxfun( @gt, cchD( tZ, : ), gbUpperZ ) & ( cchD( tZ, : ) > 0 );
    loBins( tZ, : )             = bsxfun( @lt, cchD( tZ, : ), gbLowerZ ) & ( gbLowerZ > 0 );
    nans                        = NaN( sum( nidx ), sum( idx0 ) );
    
    % determine if any bin in t_ROI is globally sig. (binary decision)
    act0                       	= any( cchD( t_ROID, : ) >= ones( sum( t_ROID ), 1 ) * gbUpper, 1 );
    sil0                        = any( cchD( t_ROID, : ) <= ones( sum( t_ROID ), 1 ) * gbLower, 1 );
    
    pred0                       = [ pred( idx1, : ); nans; pred( idx2, : ) ];
    pvalsUpper0                 = [ pvalsUpper( idx1, : ); nans; pvalsUpper( idx2, : ) ];
    pvalsLower0                 = [ pvalsLower( idx1, : ); nans; pvalsLower( idx2, : ) ];
    hiBins0                     = [ hiBins( idx1, : ); false( 1, sum( idx0 ) ); hiBins( idx2, : ) ];
    loBins0                     = [ loBins( idx1, : ); false( 1, sum( idx0 ) ); loBins( idx2, : ) ];
    gbUpper0                    = gbUpper;
    gbLower0                    = gbLower;
    gbUpperZ0                   = gbUpperZ;
    gbLowerZ0                   = gbLowerZ;
    
end

% significance (point-wise probabilities; predictor)
[ pvalsUpper, pred, pvalsLower ] = cch_conv( cch, W, convType );            % linear/non-linear filtering of the CCH

% global limits, with a conservative Bonferroni correction, on the tested range:
gbUpper                         = poissinv( 1 - alfa / nBonf, max( pred( t_ROI, : ), [], 1 ) );
gbLower                         = poissinv( alfa / nBonf, min( pred( t_ROI, : ), [], 1 ) );
gbUpperZ                        = poissinv( 1 - alfa, max( pred( tZ, : ), [], 1 ) );
gbLowerZ                        = poissinv( alfa, min( pred( tZ, : ), [], 1 ) );

% detect the bins with high/low global counts
hiBins                          = false( size( cch ) );
loBins                          = false( size( cch ) );
hiBins( t_ROI, : )              = bsxfun( @gt, cch( t_ROI, : ), gbUpper ) & ( cch( t_ROI, : ) > 0 );
loBins( t_ROI, : )              = bsxfun( @lt, cch( t_ROI, : ), gbLower ) & ( ones( nBonf, 1 ) * gbLower ) > 0;
hiBins( tZ, : )                 = bsxfun( @gt, cch( tZ, : ), gbUpperZ ) & ( cch( tZ, : ) > 0 );
loBins( tZ, : )                 = bsxfun( @lt, cch( tZ, : ), gbLowerZ ) & ( gbLowerZ > 0 );

% determine if any bin in t_ROI is globally sig. (binary decision)
act                             = any( cch( t_ROI, : ) >= ones( sum( t_ROI ), 1 ) * gbUpper, 1 );
sil                             = any( cch( t_ROI, : ) <= ones( sum( t_ROI ), 1 ) * gbLower, 1 );

% replace the results for the same-shank pairs
if sum( nidx )
    pred( :, idx0 )          	= pred0;
    pvalsUpper( :, idx0 )   	= pvalsUpper0;
    pvalsLower( :, idx0 )    	= pvalsLower0;
    hiBins( :, idx0 )         	= hiBins0;
    loBins( :, idx0 )         	= loBins0;
    gbUpper( :, idx0 )          = gbUpper0;
    gbLower( :, idx0 )          = gbLower0;
    act( :, idx0 )              = act0;
    sil( :, idx0 )              = sil0;
    gbUpperZ( :, idx0 )         = gbUpperZ0;
    gbLowerZ( :, idx0 )         = gbLowerZ0;
end

if verbose
    et                          = etime( clock, t0 );
    fprintf( 1, ' done (%0.2g s)\n', et );
end

%------------------------------------------------------------------------
% (5) scale to determine difference CCH, rate ('gain') cch, STC, and ASG
%-------------------------------------------------------------------------

if verbose
    t0                          = clock;
    fprintf( 1, '%s: Scaling and computing STCs and ASGs ... ', mfname )
end

% compute difference cch
dcch                            = cch - pred;                               % [counts]

% compute rate ('gain') cch
gcch                            = dcch ./ ( ones( m, 1 ) * nspks1 * dt );   % this will always be assymetric
gsem                            = sem ./ ( ones( m, 1 ) * nspks1 * dt );    % [spks/s]
if sum( nidx )
    gsemD                       = gsem( :, idx0 );                       	% same-shank pairs (incl. ACH)
    gsemD( nidx, : )            = NaN;
    gsem( :, idx0 )             = gsemD;
end

% allocate space
g1mat                           = NaN( 1, n * n );
g2mat                           = NaN( 1, n * n );
g1tidx                          = NaN( nsupport, n * n );
g1gcch                          = NaN( nsupport, n * n );
g2tidx                          = NaN( nsupport, n * n );
g2gcch                          = NaN( nsupport, n * n );

% compute ASG and STC
for i = 1 : size( gcch, 2 )
    
    % compute ASG and STC for a single CCH
    if idx0( i ) && sum( nidx )
        t_roi                   = t_ROID;
    else
        t_roi                   = t_ROI;
    end
    [ g1, g2, g1base, g2base ]  = calc_stc_asg( gcch( :, i ), t_roi, dt, asgMode, roiHardness );
    
    % collect ASG
    g1mat( :, i  )              = g1;
    g2mat( :, i )               = g2;
    
    % determine the support of the ASGs
    idx1t                       = ismember( sidx, g1base );
    idx2t                       = ismember( sidx, g2base );
    [ ~, idx1s ]                = intersect( g1base, sidx );
    [ ~, idx2s ]                = intersect( g2base, sidx );
    
    % assign the STC for ASGe and ASGi
    g1tidx( idx1t, i )          = g1base( idx1s );
    g1gcch( idx1t, i )          = gcch( g1base( idx1s ), i );
    g2tidx( idx2t, i )          = g2base( idx2s );
    g2gcch( idx2t, i )          = gcch( g2base( idx2s ), i );

end

if verbose
    et                          = etime( clock, t0 );
    fprintf( 1, ' done (%0.2g s)\n', et );
end

%------------------------------------------------------------------------
% (6) assign output

if verbose
    t0                          = clock;
    fprintf( 1, '%s: Organizing in structure ... ', mfname )
end

% reshape back
cch                             = reshape( cch, [ m n n ] );
pred                            = reshape( pred, [ m n n ] );
pvalsUpper                      = reshape( pvalsUpper, [ m n n ] );
pvalsLower                      = reshape( pvalsLower, [ m n n ] );
hiBins                          = reshape( hiBins, [ m n n ] );
loBins                          = reshape( loBins, [ m n n ] );
gbUpper                         = reshape( gbUpper, [ n n ] );
gbLower                         = reshape( gbLower, [ n n ] );
gbUpperZ                        = reshape( gbUpperZ, [ n n ] );
gbLowerZ                        = reshape( gbLowerZ, [ n n ] );

act                             = reshape( act, [ n n ] );
sil                             = reshape( sil, [ n n ] );
cch0                            = reshape( cch0, [ m n n ] );
sem                             = reshape( sem, [ m n n ] );
gcch                            = reshape( gcch, [ m n n ] );
gsem                            = reshape( gsem, [ m n n ] );

g1mat                           = reshape( g1mat, [ n n ] );
g2mat                           = reshape( g2mat, [ n n ] );
g1tidx                          = reshape( g1tidx, [ nsupport n n ] );
g1gcch                          = reshape( g1gcch, [ nsupport n n ] );
g2tidx                          = reshape( g2tidx, [ nsupport n n ] );
g2gcch                          = reshape( g2gcch, [ nsupport n n ] );

% 'standard' fields (conforms to old spikes2spikes output)
s2s.filebase                    = filebase;
s2s.shankclu                    = spk.shankclu( :, 1 : 2 );
s2s.nspks0                      = nspks0;
s2s.nspks                       = nspks0;
s2s.Tsec                        = Tsec;
s2s.ccg                         = cch;
s2s.t                           = t;
s2s.State                       = [];
s2s.cutoff                      = 0;
s2s.window                      = jitWindowMS;
s2s.convtype                    = convType;
s2s.alpha                       = alfa;
s2s.t_ROI                       = t_ROI;
s2s.periods                     = periods;
s2s.pred                        = pred;
s2s.pvalsUpper                  = pvalsUpper;
s2s.pvalsLower                  = pvalsLower;
s2s.hiBins                      = hiBins;
s2s.loBins                      = loBins;
s2s.gbUpper                     = gbUpper;
s2s.gbLower                     = gbLower;

% other relevant fields (new)
s2s.act                         = act;
s2s.sil                         = sil;
s2s.cch0                        = cch0;                                     % [count], raw CCH
s2s.sem                         = sem;                                      % [count], SEM of the deconvolved CCH ('ccg')
s2s.gcch                        = gcch;                                     % [spks/s], gain CCH (derived from the dcCCH)
s2s.gsem                        = gsem;                                     % [spks/s], SEM of the gain CCH
s2s.gbUpperZ                    = gbUpperZ;
s2s.gbLowerZ                    = gbLowerZ;

% STC and ASG fields (new)
s2s.g1mat                       = g1mat;
s2s.g2mat                       = g2mat;
s2s.g1tidx                      = g1tidx;
s2s.g1gcch                      = g1gcch;
s2s.g2tidx                      = g2tidx;
s2s.g2gcch                      = g2gcch;
if verbose
    et                          = etime( clock, t0 );
    fprintf( 1, ' done (%0.2g s)\n', et );
end

%------------------------------------------------------------------------
% (7) save to disk

if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( filename, 'file' ) )
    if verbose
        fprintf( 1, 'Saving %s\n', filename )
    end
    save( filename, 's2s', '-v6' );
end

return 

% EOF

% to compare run times:
tic, s2s_new = spikes2spikes_new( filebase, 'Overwrite', 0 ); toc
tic, s2s_old = spikes2spikes_old( filebase, 'Overwrite', 0 ); toc

% to compare CCHs:
figure, plot_s2s( s2s_old, [ 3 9; 3 7 ] );
figure, plot_s2s( s2s_new, [ 3 9; 3 7 ] );

% to compare plot_ss:
load( [ filebase '.sst' ], '-mat' );
figure, plot_ss( sst, [ 3 9 ], 's2s', s2s_old )
figure, plot_ss( sst, [ 3 9 ], 's2s', s2s_new )
