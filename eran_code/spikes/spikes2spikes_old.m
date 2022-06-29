% spikes2spikes         Computes all pair-wise CCGs and determines significance
%
% call          s2s = spike2spikes( filebase, BinSizeMS, halfWidthMS, jitHalfWindowMS, roiMS...
%                   , convType, alpha, CutOffMS, clustr, Overwrite, periods, Electrodes )
%
% BinSizeMS, halfWidthMS        {1,50} for CCH computation
% jitHalfWindowMS               {5} half width for jittering/convolution
% roiMS                         {[-5 5 ]} region of interest for global significance estimation
% convType                      {'gauss'}, or 'rect', 'triang', or 'jitter'
% alpha                         {0.001} p-value for detection of extreme values (this is bonferonni corrected in the ROI)
%                                NOTE that for jittering, the number of jitters is automatically
%                                determined to be 5 times 1/alpha (to be enable to estimate tails
%                                properly)
% CutOffMS                      {5} to burst filter trains before CCH computation
% clustr                        {'clc'} which clu file to load
% Overwrite                     {1} to recompute and overwrite; -1 to just recompute 
% periods                       {-1}. periods (in samples) to keep spikes from. overloaded: 
%                               -1 indicates to remove spikes during segments specified in val files; 
%                               -2 also removes spikes during HFOs (see get_hfo_times.m)
%                               -3 removes from segments indicates in val and in HFOs 
%                               -XX indicates to remove only during segments in val.tXX (can be a vector for multiple triggers).
%                               right now the combination of the two is not suppported (i.e. keep spikes in certain periods but not during
%                                   stimulation etc). also not relevant for merged files 
% 
% calls         FileExists, ParseArgPairs, LoadXml, LoadCluRes, LoadStims
%               get_hfo_times, sortranges, resampleranges, uniteranges, inranges
%               CCG, CCG_jitter, cch_conv
%
% run time      6 sec for a 1 hour recording of 1,000,000 spikes from 70 units
%               (most time spent on loading from HD.. less than 3 sec for
%               computations; jittering was 2 hours)
%
% nspks         may be lower than clu d.t. (1) burst filtering (2) trigger
%               filtering (3) clc filtering (all are optional, though)

% 11-jan-12 ES (based on Spikes2Spikes)

% revisions
% 15-jan-12 various modifications (care of zero-lag bin)
% 16-jan-12 additional modifications
% 24-jan-12 load/save decision
% 27-jan-12 remove spikes during stimulation (according to *.val.* files)
% 02-mar-12 changed so writes if no file OR overwrite flag; 
%               changed default overwrite to 0; 
%               added module to support cell array (datestr/fnum) filebase format
% 13-may-12 (1) changed so Overwrite -1 does not write anymore
%           (2) fixed bug in shankclu in case of periods-removed clusters
% 09-jun-12 (1) changed format to support r2009b
%           (2) handled case of no spikes at all
% 17-jul-12 (1) modification of input case cell array, non-merged file 
% 19-may-13 (1) datenum2filebase integrated
%           (2) suffix external specification optional
%           (3) LoadVals -> LoadStims
%           (4) get_hfo_times integrated (periods == -2; -3 excludes both stimulus and hfo times)
%           (5) Overwrite logic extended to -2/-1/0/1 
% 09-oct-13 (1) W rounded to support arbitrary sampling rates
% 01-aug-19 (1) changed DefaultArgs to ParseArgPairs
% 17-aug-19 cleaned up

% there are multiple problems with the jittering procedure:
% (1) too slow..
% (2) straight-forward jittering violates same-shank dead-time
% (3) globalband bonferonni corrects by n(cch) which is way too conservative
% which is why we use the convolution method.
% (1) is solved automatically
% (2) is solved by ignoring the zero-lag bin from same-shank CCH 
% (3) is solved by taking the global bands only in the region of interest
% (still bonf. corrected)

function s2s = spikes2spikes( filebase, varargin )

verbose                 = 1;

[ BinSizeMS, halfWidthMS, jitWindowMS, roiMS, convType, alpha, CutOffMS ...
    , clustr, Overwrite, periods, shanknums, suffix ] = ParseArgPairs( ...
    { 'BinSizeMS', 'halfWidthMS', 'jitWindowMS', 'roiMS', 'convType', 'alpha', 'CutOffMS' ...
    , 'clustr', 'Overwrite', 'periods', 'shanknums', 'suffix' } ...
    , { 1, 50, 5, [ -5 5 ], 'gauss', 0.001, 5 ...
    , 'clu', -2, -1, [], '' }...
    , varargin{ : } );

% load/save
if isa( filebase, 'cell' ) && length( filebase ) == 2
    filebase            = datenum2filebase( filebase );
end
if isempty( suffix )
    SaveFn              = [ filebase '.s2s' ];
    if strcmp( convType, 'jitter' )
        SaveFn          = [ SaveFn '.jit' ];
    end
else
    SaveFn              = [ filebase '.' suffix ];
end
if FileExists( SaveFn ) && Overwrite < 0 && nargout > 0
    if verbose
        fprintf( 1, 'Loading %s\n', SaveFn )
    end
    load(SaveFn,'-MAT');
    return
end

% get the parameters
if ~exist( 'par', 'var' )
    par                 = LoadXml( sprintf( '%s.xml', filebase ) );
end
if isempty( shanknums )
    shanknums           = 1 : par.nElecGps;
end
SpikesFs                = par.SampleRate; % sampling rate (samples/s)
nBins                   = halfWidthMS/BinSizeMS; % number of bins on each side of the CCH
BinSize                 = round( BinSizeMS * SpikesFs / 1000 ); % ms -> samples
CutOff                  = round( CutOffMS * SpikesFs / 1000 ); % ms -> samples
W                       = ceil( 2 * jitWindowMS / ( BinSize / SpikesFs * 1000 ) + 1 ); % ms of half-window -> nbins for full window

% load 
[ Res, Clu, Map ]       = LoadCluRes( filebase, shanknums, [], clustr );
if isempty( Clu )
    fprintf( 1, 'No spikes in %s, shanks %s\n', filebase, num2str( shanknums ) )
    s2s                 = [];
    return
end

% select spikes
if periods( 1 ) < 0
    % remove spikes during stimulation
    [ Vals, Trigs ]     = LoadStims( filebase );
    if periods( 1 ) == -1 || periods( 1 ) == -3 % any trigger
        tidx            = true( size( Trigs ) );
    elseif periods( 2 ) == -2 % specific subset
        tidx            = ismember( Trigs, abs( periods ) );
    else 
        tidx            = [];
    end
    if isempty( tidx )
        uvals1          = [];
    else
        uvals1          = uniteranges( Vals( tidx, 1 : 2 ) ); % combine all the segments
    end
    if periods( 1 ) == -2 || periods( 1 ) == -3 % remove HFO times
        hperiods        = get_hfo_times( filebase );
        hperiods        = sortranges( hperiods );
        eegFs           = par.lfpSampleRate; 
        uvals2          = resampleranges( hperiods, SpikesFs, eegFs );
    else
        uvals2          = [];
    end
    uvals               = uniteranges( uvals1, uvals2 );
    if ~isempty( uvals )
        ridx            = inranges( Res, uvals ); % remove any spike that is in any segment
        Res( ridx )     = [];
        Clu( ridx )     = [];
    end
elseif ~isempty( periods ) && size( periods, 2 ) == 2
    % keep only spikes during specified periods
    kidx                = inranges( Res, periods );
    Res                 = Res( kidx );
    Clu                 = Clu( kidx );
end
if isempty( Clu )
    fprintf( 1, 'No spikes in %s, shanks %s\n', filebase, num2str( shanknums ) )
    s2s                 = [];
    return
end

% compute raw ACHs, CCHs
uClu                    = unique( Map( :, 1 ) );% in case unit only spikes during masked epochs
[ ccg, t ]              = CCG( Res, Clu, BinSize, nBins, SpikesFs, uClu, 'count' ); % ACH
t                       = t( : );
t_ROI                   = t >= roiMS(1) & t <= roiMS(2);
nBonf                   = sum( t_ROI );
[ m, n, p ]             = size( ccg );

% dilute trains and recompute CCHs
nspks0                  = hist( Clu, uClu )';
if CutOff > 0
    rmv                 = [];
    for i               = uClu( : ).'
        idx             = find( Clu == i );
        isis            = diff( Res( idx ) );
        rmvidx          = find( isis <= CutOff ) + 1;
        rmv             = [ rmv; idx( rmvidx ) ];
    end
    Res( rmv )          = [];
    Clu( rmv )          = [];
    ccg0                = ccg;
    ccg                 = CCG( Res, Clu, BinSize, nBins, SpikesFs, uClu, 'count' ); % ACH
    for i               = 1 : n % keep the ACH for the non-dilulted trains
        ccg( :, i, i )  = ccg0( :, i, i );
    end
end

% jitter/convolve for significance
if strcmp( convType, 'jitter' )
    
    pred                = zeros( 2 * nBins + 1, n, p );
    hiBins              = pred;
    loBins              = pred;
    gbUpper             = zeros( p, n );
    gbLower             = zeros( p, n );
    for i               = 1 : n
        for j           = i : p % do no compute twice..
            if verbose
                fprintf( 1, '%dx%d\n', uClu( i ), uClu( j ) )
            end
            [ ~, ~, hiBins( :, i, j ), loBins( :, i, j ), ccgjMtx ] = CCG_jitter( Res...
                , Clu, uClu( i ), uClu( j ), BinSize, nBins, jitWindowMS, 5/alpha, alpha, 0 );
            pred( :, i, j )     = mean( ccgjMtx, 2 );
            smax                = sort( max( ccgjMtx ) ); 
            smin                = sort( min( ccgjMtx ) ); 
            gidx                = ( alpha * [ 1 -1 ] + [ 0 1 ] ) * size( ccgjMtx, 2 ); 
            gbUpper( i, j )     = smax( gidx( 2 ) );
            gbLower( i, j )     = smin( gidx( 1 ) );
            % fill the other half
            pred( :, j, i )     = flipud( pred( :, i, j ) );
            hiBins( :, j, i )   = hiBins( :, i, j );
            loBins( :, j, i )   = loBins( :, i, j );
            gbUpper( j, i )     = gbUpper( i, j );
            gbLower( j, i )     = gbLower( i, j );
        end
    end
    pvalsUpper                  = [];
    pvalsLower                  = [];
    
else
    
    % reshape
    cch                         = reshape( ccg, [ m n * p ] );
    
    % identify same-shank pairs
    shank                       = Map( :, 2 );
    idx0                        = shank( ( 1 : n )' * ones( 1, p ) ) == shank( ones( n, 1 ) * ( 1 : p ) );
    idx0                        = reshape( idx0, [ 1 n * p ] ); % same shank (including same unit)
    
    % fix zero-lag bins
    if 1
        % compute without the zero-lag bins, then replace with NaNs
        cch0                                = cch( :, idx0 );
        cch0( nBins + 1, : )                = [];
        [ pvalsUpper, pred, pvalsLower ]    = cch_conv( cch0, W, convType );
        gbUpper                             = poissinv( 1 - alpha / nBonf, max( pred( t_ROI, : ), [], 1 ) );
        gbLower                             = poissinv( alpha / nBonf, min( pred( t_ROI, : ), [], 1 ) );
        hiBins                              = false( size( cch0 ) );
        loBins                              = false( size( cch0 ) );
        hiBins( t_ROI, : )                  = bsxfun( @gt, cch0( t_ROI, : ), gbUpper ) & ( cch0( t_ROI, : ) > 0 );
        loBins( t_ROI, : )                  = bsxfun( @lt, cch0( t_ROI, : ), gbLower ) & ( ones( nBonf, 1 ) * gbLower ) > 0;
        nans                                = NaN * ones( 1, sum( idx0 ) );
        pvalsUpper0                         = [ pvalsUpper( 1 : nBins, : ); nans; pvalsUpper( nBins + 1 : 2 * nBins, : ) ];
        pred0                               = [ pred( 1 : nBins, : ); nans; pred( nBins + 1 : 2 * nBins, : ) ];
        pvalsLower0                         = [ pvalsLower( 1 : nBins, : ); nans; pvalsLower( nBins + 1 : 2 * nBins, : ) ];
        hiBins0                             = [ hiBins( 1 : nBins, : ); false( 1, sum( idx0 ) ); hiBins( nBins + 1 : 2 * nBins, : ) ];
        loBins0                             = [ loBins( 1 : nBins, : ); false( 1, sum( idx0 ) ); loBins( nBins + 1 : 2 * nBins, : ) ];
        
    else
        % average peri-zero bins for same-shank (affects the conv)
        cch( nBins + 1, idx0 )              = round( mean( cch( nBins + [ 0 2 ], idx0 ), 1 ) );
    end
    
    % significance (point-wise probabilities; predictor)
    [ pvalsUpper, pred, pvalsLower ]        = cch_conv( cch, W, convType );
    
    % global limits, with a conservative Bonferroni correction, on the tested range:
    gbUpper                                 = poissinv( 1 - alpha / nBonf, max( pred( t_ROI, : ), [], 1 ) );
    gbLower                                 = poissinv( alpha / nBonf, min( pred( t_ROI, : ), [], 1 ) );
    
    % detect the bins with high/low global counts
    hiBins                                  = false( size( cch ) );
    loBins                                  = false( size( cch ) );
    hiBins( t_ROI, : )                      = bsxfun( @gt, cch( t_ROI, : ), gbUpper ) & ( cch( t_ROI, : ) > 0 );
    loBins( t_ROI, : )                      = bsxfun( @lt, cch( t_ROI, : ), gbLower ) & ( ones( nBonf, 1 ) * gbLower ) > 0;
    
    % replace the results for the same-shank pairs
    if 1
        pvalsUpper( :, idx0 )               = pvalsUpper0;
        pred( :, idx0 )                     = pred0;
        pvalsLower( :, idx0 )               = pvalsLower0;
        hiBins( :, idx0 )                   = hiBins0;
        loBins( :, idx0 )                   = loBins0;
    end
    
    % reshape
    pred                                    = reshape( pred, [ m n p ] );
    pvalsUpper                              = reshape( pvalsUpper, [ m n p ] );
    pvalsLower                              = reshape( pvalsLower, [ m n p ] );
    hiBins                                  = reshape( hiBins, [ m n p ] );
    loBins                                  = reshape( loBins, [ m n p ] );
    gbUpper                                 = reshape( gbUpper, [ n p ] );
    gbLower                                 = reshape( gbLower, [ n p ] );
    
end

% store
s2s.filebase                    = filebase;
s2s.shankclu                    = Map( uClu, 2 : 3 );
s2s.nspks0                      = nspks0;
s2s.nspks                       = hist( Clu, uClu )';
s2s.Tsec                        = ( Res( end ) - Res( 1 ) ) / SpikesFs;
s2s.ccg                         = ccg;
s2s.t                           = t;
s2s.State                       = [];
s2s.cutoff                      = CutOffMS;
s2s.window                      = jitWindowMS;
s2s.convtype                    = convType;
s2s.alpha                       = alpha;
s2s.t_ROI                       = t_ROI;
s2s.periods                     = periods;
s2s.pred                        = pred;
s2s.pvalsUpper                  = pvalsUpper;
s2s.pvalsLower                  = pvalsLower;
s2s.hiBins                      = hiBins;
s2s.loBins                      = loBins;
s2s.gbUpper                     = gbUpper;
s2s.gbLower                     = gbLower;

if Overwrite == 1 || ( Overwrite ~= -1 && ~exist( SaveFn, 'file' ) )
    if verbose
        fprintf( 1, 'Saving %s\n', SaveFn )
    end
    save( SaveFn, 's2s', '-v6' );
end


return

% EOF

