% check_mono         check short-term synchrony/mono-synaptic functional connectivity
%
% get a list of pairs, each with the time lag and p-value of the sign.
%
% call              mono = check_mono( filebase, gidx, minCounts, nInhBins, Optional, convType )
%                        , prePrune, excludeSameShank, suffix )
%
% ARGUMENTS
% filebase - can receive either a filebase or an s2s structure
%
% gidx - logical vector, output of check_cluster_quality, for pre/post-pruning
%       alternatively isolation level class (default, {'B'} )
% minCounts - mean counts/CCH bin in the mono-synaptic ROI. defaults to {5}
% nInhBins - minimal number of consecutive inh bins to define a CCH as inh; {2}
%           excitatory connectivity 
% Optional - [ 1 1 0 ]. 
%   First element - prevent same-shank inh. pairs to be mutually-inhibitory (could be ref. period if over-sorted)
%   Second element - prevent a pair from begin exc and inh on same side of CCH (could be for HFOs)
%   Third element - prevent a cell from being exc and inh (could be if HFO or poorly isolated)
% convType - {'gauss'}; alternatives include 'jitter' (assuming an
%               *s2s.jit file exists); alternate convolution windows (e.g.
%               'rect'), and 'none'; in the latter case, significance is
%               computed locally based on a predictor generated from the
%               peri-mono bins (default is <-10 and >10 ms)
% prePrune - flag; {1}; Post-pruning is always done
% excludeSameShank - {0};
% suffix - {''}; if not empty, replaces the 's2s' (or 's2s.jit') suffix
% ccTH - 
%
% DOES
% pairs are classified as either exc, inh, sync, or de-sync:
% 1. exc - the reference unit in a CCH with a short-latency pos.lag peak 
% 2. inh - the reference unit in a CCH with a short-latency pos.lag trough
% 3. both - the reference unit in a CCH with a short-latency pos.lag peak+trough
% 4. unc - unclassified units (no peaks/referenced unit/zero-lag peak/trough)
% 
% then the lists are pruned as follows:
% 1. pairs with sync + exc/inh are classified as only sync (they are not strictly exc/inh)
% 2. pairs with de-sync + inh are classified as only de-sync (may be wide)
% 
%
% see also spikes2spikes, plot_s2s, check_cluster_quality

% 15-jan-12 ES

% revision
% 24-jan-12 pruning algorithms modified
% 30-jan-12 gidx input added
% 16-may-13 (1) optional prePruning (default - 1)
%           (2) local significance computed by CCH baseline
%           (3) optional same-shank exclusion (default - 0)
%           (4) optional call to determine_units (default - 'B')
%           (5) suffix added
%           (6) flexible argument handling
%           (7) excited and inhibited lists added

% called by:
% -plot_ss
% -wnAnalysisCelltypeClassification

function mono = check_mono( filebase, gidx, varargin )

% FileExists 

if nargin < 2 || isempty( gidx )
    gidx                        = 'B';
end
[ minCounts, nInhBins, Optional, convType, prePrune, excludeSameShank, suffix, ccTH ] = ParseArgPairs(...
    { 'minCounts', 'nInhBins', 'Optional', 'convType', 'prePrune', 'excludeSameShank', 'suffix', 'ccTH' }...
    , { 5, 2, [ 1 1 0 ], 'gauss', 1, 0, '', 0.9 }...
    , varargin{ : } );

if isa( filebase, 'char' )
    if isempty( suffix )
        LoadFn                  = [ filebase '.s2s' ];
        if strcmp( convType, 'jitter' )
            LoadFn              = [ LoadFn '.jit' ];
        end
    else
        LoadFn                  = [ filebase '.' suffix ];
    end
    if FileExists( LoadFn )
        load( LoadFn, '-mat' )
    else
        % recompute..
        if isempty( suffix )
            s2s                 = spikes2spikes( filebase );
        elseif isequal( suffix, 's2s.nohfo' )
            s2s                 = spikes2spikes( filebase, 'clustr', 'clu', 'Overwrite', -2, 'periods', -3, 'suffix', 's2s.nohfo' );
        end
    end
elseif isa( filebase, 'struct' )
    s2s                         = filebase;
    clear filebase
end


% accomodate old formats
if isfield( s2s, 'ElClu' )
    s2s.shankclu = s2s.ElClu;
    s2s = rmfield( s2s, 'ElClu' );
end
if isfield( s2s, 'FileBase' )
    s2s.filebase = s2s.FileBase;
    s2s = rmfield( s2s, 'FileBase' );
end
if ~isfield( s2s, 'nspks0' )
    s2s.nspks0 = s2s.nspks;
end
if isempty( s2s.pvalsUpper )
    s2s.pvalsUpper = NaN * ones( size( s2s.ccg ) );
end
if isempty( s2s.pvalsLower )
    s2s.pvalsLower = NaN * ones( size( s2s.ccg ) );
end
if size( s2s.hiBins, 1 ) == sum( s2s.t_ROI )
    hiBins = zeros( size( s2s.ccg ) );
    hiBins( s2s.t_ROI, :, : ) = s2s.hiBins;
    s2s.hiBins = hiBins;
end
if size( s2s.loBins, 1 ) == sum( s2s.t_ROI )
    loBins = zeros( size( s2s.ccg ) );
    loBins( s2s.t_ROI, :, : ) = s2s.loBins;
    s2s.loBins = loBins;
end

% select well-isolated cells
if isa( gidx, 'char' ) && length( gidx ) == 1
    ilevel = gidx;
    if exist( 'filebase', 'var' )
        shankclu = determine_units( filebase, [], ilevel );
        shankclu0 = s2s.shankclu;
        gidx = ismember( shankclu0( :, 1 : 2 ), shankclu( :, 1 : 2 ), 'rows' );
    else
        gidx = true( size( s2s.shankclu, 1 ), 1 );
    end
end

% pre-prune
if prePrune && ~isempty( gidx )
    if length( gidx ) ~= size( s2s.shankclu, 1 )
        % this is a problem. hack by n (might be way off..)
        ok = 0;
        if ~ok
            shankcluD = determine_units( filebase, [], 'B' );
            if size( shankcluD, 1 ) == length( gidx )
                ok = 1;
            end
        end
        if ~ok
            shankcluD = determine_units( filebase, [], 'C' );
            if size( shankcluD, 1 ) == length( gidx )
                ok = 1;
            end
        end
        if ~ok
            shankcluD = determine_units( filebase, [], 'D' );
            if size( shankcluD, 1 ) == length( gidx )
                ok = 1;
            end
        end
        if ~ok
            error( 'mismatch..' )
        end
        gidx = ismember( s2s.shankclu, shankcluD( gidx, 1 : 2 ), 'rows' );
    end
    s2s.shankclu = s2s.shankclu( gidx, : ); 
    s2s.nspks = s2s.nspks( gidx, : ); 
    s2s.nspks0 = s2s.nspks0( gidx, : ); 
    s2s.ccg = s2s.ccg( :, gidx, gidx ); 
    s2s.pred = s2s.pred( :, gidx, gidx ); 
    s2s.pvalsUpper = s2s.pvalsUpper( :, gidx, gidx ); 
    s2s.pvalsLower = s2s.pvalsLower( :, gidx, gidx ); 
    s2s.hiBins = s2s.hiBins( :, gidx, gidx ); 
    s2s.loBins = s2s.loBins( :, gidx, gidx ); 
    gidx = true( sum( gidx ), 1 );
end
[ m n p ] = size( s2s.ccg );

% excitation/synchrony/inhibition range
%idx1 = s2s.t_ROI & s2s.t < 0;
idx2 = s2s.t_ROI & s2s.t == 0; % sync/de-sync
idx3 = s2s.t_ROI & s2s.t > 0;  % exc/inh

% compute significance locally (use the mean -50:-10, 10:50 as a predictor)
sameshank = ( s2s.shankclu( :, 1 ) * ones( 1, n ) ) == ( ones( n, 1 ) * s2s.shankclu( :, 1 )' );
sameshank = reshape( sameshank, [ 1 n * p ] );
if strcmp( convType, 'none' )
    ccg = reshape( s2s.ccg, [ m n * p ] );
    t_ROI = s2s.t_ROI;
    nBonf = sum( t_ROI );
    sameshank = ( s2s.shankclu( :, 1 ) * ones( 1, n ) ) == ( ones( n, 1 ) * s2s.shankclu( :, 1 )' );
    sameshank = reshape( sameshank, [ 1 n * p ] );
    cidx = s2s.t == 0;
    ccg( cidx, sameshank ) = NaN;
    TH = 2 * max( abs( s2s.t( t_ROI ) ) );
    bidx = s2s.t < -TH | s2s.t > TH;
    mm = ones( length( s2s.t ), 1 ) * mean( ccg( bidx, : ) );
    gbUpper = poissinv( 1 - s2s.alpha / nBonf, mm( 1, : ) );
    gbLower = poissinv( s2s.alpha / nBonf, mm( 1, : ) );
    hiBins = false( size( ccg ) );
    loBins = false( size( ccg ) );
    hiBins( t_ROI, : ) = bsxfun( @gt, ccg( t_ROI, : ), gbUpper ) & ( ccg( t_ROI, : ) > 0 );
    loBins( t_ROI, : ) = bsxfun( @lt, ccg( t_ROI, : ), gbLower ) & ( ones( nBonf, 1 ) * gbLower ) > 0;
    pvalsUpper = 1 - poisscdf( ccg - 1, mm ) - poisspdf( ccg, mm ) * 0.5;
    pvalsLower = 1 - pvalsUpper;
else
    % the results of Spikes2SpikesConv
    hiBins = reshape( s2s.hiBins, [ m n * p ] );
    loBins = reshape( s2s.loBins, [ m n * p ] );
    pvalsUpper = reshape( s2s.pvalsUpper, [ m n * p ] );
    pvalsLower = reshape( s2s.pvalsLower, [ m n * p ] );
end

% pair members
u1 = [ 1 : n ]' * ones( 1, p );
u2 = ones( n, 1 ) * [ 1 : p ];
pairs = [ u1( : ) u2( : ) ];

% valid pairs - sufficient counts in the relevant range + not ACH
validpairs = squeeze( sum( s2s.ccg( s2s.t_ROI, :, : ), 1 ) ) / sum( s2s.t_ROI ) > minCounts;
validpairs = reshape( validpairs, [ 1 n * p ] ) & [ pairs( :, 1 ) ~= pairs( :, 2 ) ]';
if excludeSameShank
    validpairs( sameshank ) = 0;
end

% % compare auto and cross for same-shank
% if ~excludeSameShank
%     ccg = reshape( s2s.ccg, [ m n * p ] );
%     didx = diagidx( zeros( n, p ) );
%     ach = ccg( :, didx );
% end

% BASIC SOLUTIONS
% because ccg contains each pair twice, don't need to check idx1 and all lags are zero/positive
% if a pair has two peaks (or two dips), then it will simply be listed twice
% u1 excites u2 (for any number of bins; take the lag as the minimal p-value):
idxExc = find( sum( hiBins( idx3, : ), 1 ) > 0 & validpairs );
pairsExc = zeros( length( idxExc ), 4 );
j = 0;
for i = idxExc
    j = j + 1;
    fidx = find( hiBins( :, i ) & idx3 );
    [ pval idx ] = min( pvalsUpper( fidx, i ) );
    lag = s2s.t( fidx( idx ) );
    pairsExc( j, : ) = [ pairs( i, : ) pval lag ];
end
% u1 inhibits u2 (for at least nInhBins consecutive bins; take the lag as the minimal p-value in the longest stretch):
idxInh = find( sum( loBins( idx3, : ), 1 ) > 0 & validpairs );
pairsInh = zeros( length( idxInh ), 4 );
j = 0;
for i = idxInh
    fidx = find( loBins( :, i ) & idx3 );
    mat = parse( fidx );
    [ nb, midx ] = max( diff( mat, [], 2 ) + 1 );
    if nb >= nInhBins
        fidx = mat( midx, 1 ) : mat( midx, 2 );
    else
        continue
    end
    [ pval idx ] = min( pvalsLower( fidx, i ) );
    lag = s2s.t( fidx( idx ) );
    j = j + 1;
    pairsInh( j, : ) = [ pairs( i, : ) pval lag ];
end
pairsInh( j + 1 : length( idxInh ), : ) = [];
% u1,u2 synchronized:
idxSync = find( sum( hiBins( idx2, : ), 1 ) > 0 );
pairsSync = zeros( length( idxSync ), 4 );
j = 0;
for i = idxSync
    j = j + 1;
    fidx = find( hiBins( :, i ) & idx2 );
    [ pval idx ] = min( pvalsUpper( fidx, i ) );
    lag = s2s.t( fidx( idx ) );
    pairsSync( j, : ) = [ pairs( i, : ) pval lag ];
end
% u1,u2 desynchronized (may indicate splitting a well-isolated unit)
idxDesync = find( sum( loBins( idx2, : ), 1 ) > 0 & validpairs );
pairsDesync = zeros( length( idxDesync ), 4 );
j = 0;
for i = idxDesync
    j = j + 1;
    fidx = find( loBins( :, i ) & idx2 );
    [ pval idx ] = min( pvalsLower( fidx, i ) );
    lag = s2s.t( fidx( idx ) );
    pairsDesync( j, : ) = [ pairs( i, : ) pval lag ];
end

% PRUNE:
% (0) same-shank pairs cannot have CCH too similar to the ACHs
if ~excludeSameShank
    shanks = s2s.shankclu( :, 1 );
    nbins = floor( length( s2s.t ) / 2 );
    
    sameExc = find( diff( shanks( pairsExc( :, 1 : 2 ) ), 1, 2 ) == 0 );
    cc = zeros( size( sameExc, 1 ), 4 );
    for i = 1 : length( sameExc )
        ij = pairsExc( sameExc( i ), 1 : 2 );
        cch = s2s.ccg( :, ij( 1 ), ij( 2 ) ) * ones( 1, 2 );
        achs = [ s2s.ccg( :, ij( 1 ), ij( 1 ) ) s2s.ccg( :, ij( 2 ), ij( 2 ) ) ];
        %cc( i, : ) = calc_pearson( cch, achs );
        cch = [ cch( 1 : nbins, : ) cch( nbins + 2 : end, : ) ];
        achs = [ achs( 1 : nbins, : ) achs( nbins + 2 : end, : ) ];
        cc( i, : ) = calc_pearson( cch, achs );
        %     [ cc2( i, 1 ) pval( i, 1 ) ] = calc_spearman( cch( :, 1), achs( :, 1 ), 1000 );
        %     [ cc2( i, 2 ) pval( i, 2 ) ] = calc_spearman( cch( :, 2), achs( :, 2 ), 1000 );
    end
    i1 = any( cc > ccTH, 2 );
    pairsExc( sameExc( i1 ), : ) = [];
    
    sameInh = find( diff( shanks( pairsInh( :, 1 : 2 ) ), 1, 2 ) == 0 );
    cc = zeros( size( sameInh, 1 ), 4 );
    for i = 1 : length( sameInh )
        ij = pairsInh( sameInh( i ), 1 : 2 );
        cch = s2s.ccg( :, ij( 1 ), ij( 2 ) ) * ones( 1, 2 );
        achs = [ s2s.ccg( :, ij( 1 ), ij( 1 ) ) s2s.ccg( :, ij( 2 ), ij( 2 ) ) ];
        cch = [ cch( 1 : nbins, : ) cch( nbins + 2 : end, : ) ];
        achs = [ achs( 1 : nbins, : ) achs( nbins + 2 : end, : ) ];
        cc( i, : ) = calc_pearson( cch, achs );
    end
    i1 = any( cc > ccTH, 2 );
    pairsInh( sameInh( i1 ), : ) = [];
    
    sameSync = find( diff( shanks( pairsSync( :, 1 : 2 ) ), 1, 2 ) == 0 );
    cc = zeros( size( sameSync, 1 ), 4 );
    for i = 1 : length( sameSync )
        ij = pairsSync( sameSync( i ), 1 : 2 );
        cch = s2s.ccg( :, ij( 1 ), ij( 2 ) ) * ones( 1, 2 );
        achs = [ s2s.ccg( :, ij( 1 ), ij( 1 ) ) s2s.ccg( :, ij( 2 ), ij( 2 ) ) ];
        cch = [ cch( 1 : nbins, : ) cch( nbins + 2 : end, : ) ];
        achs = [ achs( 1 : nbins, : ) achs( nbins + 2 : end, : ) ];
        cc( i, : ) = calc_pearson( cch, achs );
    end
    i1 = any( cc > ccTH, 2 );
    pairsSync( sameSync( i1 ), : ) = [];
    
    sameDesync = find( diff( shanks( pairsDesync( :, 1 : 2 ) ), 1, 2 ) == 0 );
    cc = zeros( size( sameDesync, 1 ), 4 );
    for i = 1 : length( sameDesync )
        ij = pairsDesync( sameDesync( i ), 1 : 2 );
        cch = s2s.ccg( :, ij( 1 ), ij( 2 ) ) * ones( 1, 2 );
        achs = [ s2s.ccg( :, ij( 1 ), ij( 1 ) ) s2s.ccg( :, ij( 2 ), ij( 2 ) ) ];
        cch = [ cch( 1 : nbins, : ) cch( nbins + 2 : end, : ) ];
        achs = [ achs( 1 : nbins, : ) achs( nbins + 2 : end, : ) ];
        cc( i, : ) = calc_pearson( cch, achs );
    end
    i1 = any( cc > ccTH, 2 );
    pairsDesync( sameDesync( i1 ), : ) = [];
end


% (1.1) Exc pairs cannot be sync (on same side of the CCH; another way is to keep the larger effect)
[ ~, i1 ] = intersect( pairsExc( :, 1 ) + 1e3 * pairsExc( :, 2 ), pairsSync( :, 1 ) + 1e3 * pairsSync( :, 2 ) );
pairsExc( i1, : ) = [];
% (1.2) Inh pairs cannot be sync (on smae side of the CCH)
[ ~, i1 ] = intersect( pairsInh( :, 1 ) + 1e3 * pairsInh( :, 2 ), pairsSync( :, 1 ) + 1e3 * pairsSync( :, 2 ) );
pairsInh( i1, : ) = [];
% (2) Inh pairs cannot be de-sync
[ ~, i1 ] = intersect( pairsInh( :, 1 ) + 1e3 * pairsInh( :, 2 ), pairsDesync( :, 1 ) + 1e3 * pairsDesync( :, 2 ) );
pairsInh( i1, : ) = [];
% (3) sync pairs are all listed twice
% change this to a proper method..
i = 1;
while i < size( pairsSync, 1 )
    idx = find( pairsSync( :, 2 ) == pairsSync( i, 1 ) & pairsSync( :, 1 ) == pairsSync( i, 2 ) );
    if ~isempty( idx )
        pairsSync( idx, : ) = [];
    end
    i = i + 1;
end
% (4) de-sync pairs are also listed twice
i = 1;
while i < size( pairsDesync, 1 )
    idx = find( pairsDesync( :, 2 ) == pairsDesync( i, 1 ) & pairsDesync( :, 1 ) == pairsDesync( i, 2 ) );
    if ~isempty( idx )
        pairsDesync( idx, : ) = [];
    end
    i = i + 1;
end
if Optional( 1 )
    % (5.1) Same-shank Inh pairs cannot be mutually inhibitory (may indicate over-sorting)
    [ ~, i1 ] = intersect( pairsInh( :, 1 ) + 1e3 * pairsInh( :, 2 ), pairsInh( :, 2 ) + 1e3 * pairsInh( :, 1 ) );
    shanks = s2s.shankclu( :, 1 );
    inhshanks = shanks( pairsInh( :, 1 : 2 ) );
    if numel( inhshanks ) == 2
        inhshanks = inhshanks( : ).';
    end
    if ~isempty( inhshanks )
        ix = intersect( find( inhshanks( :, 1 ) == inhshanks( :, 2 ) ), i1 );
        pairsInh( ix, : ) = [];
    end
    % (5.2) Same-shank Exc pairs cannot be mutually excitatory (may indicate over-sorting)
    [ ~, i1 ] = intersect( pairsExc( :, 1 ) + 1e3 * pairsExc( :, 2 ), pairsExc( :, 2 ) + 1e3 * pairsExc( :, 1 ) );
    shanks = s2s.shankclu( :, 1 );
    excshanks = shanks( pairsExc( :, 1 : 2 ) );
    if numel( excshanks ) == 2
        excshanks = excshanks( : ).';
    end
    if ~isempty( excshanks )
        ix = intersect( find( excshanks( :, 1 ) == excshanks( :, 2 ) ), i1 );
        pairsExc( ix, : ) = [];
    end
end
if Optional( 2 )
    % (6) The same pair cannot be exc & inh at the same side of the CCH (NOTE: may be in the case of HFO)
    [ ~, i1 i2 ] = intersect( pairsExc( :, 1 ) + 1e3 * pairsExc( :, 2 ), pairsInh( :, 1 ) + 1e3 * pairsInh( :, 2 ) );
    if ~isempty( i1 )
%         fprintf( '%s: check this!\n', upper( mfilename ) )
%         disp( i1 )
%         keyboard
        rmvExc = [];
        rmvInh = [];
        for i = 1 : length( i1 )
            pidx = pairsExc( i1( i ), 1 : 2 );
            idx = ( pidx( 2 ) - 1 ) * n + pidx( 1 );
            hi = find( hiBins( :, idx ) );
            lo = find( loBins( :, idx ) );
            hi( hi < m / 2 ) = [];
            lo( lo < m / 2 ) = [];
            if min( hi ) < min( lo )
                rmvInh = [ rmvInh; i2( i ) ];
            else
                rmvExc = [ rmvExc; i1( i ) ];
            end
        end
        pairsExc( rmvExc, : ) = [];
        pairsInh( rmvInh, : ) = [];
    end
%     pairsExc( i1, : ) = [];
%     pairsInh( i2, : ) = [];
    %pairsInh( union( i2, ix ), : ) = [];
    % right now both are pruned. this is not the best strategy, because 
    % another way to go here is to choose the first/larger effect
end
if Optional( 3 )
    % (7) Exc cells cannot be Inh (NOTE: may be in the case of HFO/poor isolation)
    ix = intersect( unique( pairsExc( :, 1 ) ), unique( pairsInh( :, 1 ) ) );
    pairsExc( ismember( pairsExc( :, 1 ), ix ), : ) = [];
    pairsInh( ismember( pairsInh( :, 1 ), ix ), : ) = [];
end

% SUMMARIZE
mono.filebase = s2s.filebase;
mono.shankclu = s2s.shankclu;
mono.minCounts = minCounts;
mono.nInhBins = nInhBins;
mono.Optional = Optional;
mono.pairsExc = pairsExc; 
mono.pairsInh = pairsInh; 
mono.pairsSync = pairsSync; 
mono.pairsDesync = pairsDesync;
exc = unique( mono.pairsExc( :, 1 ) );
inh = unique( mono.pairsInh( :, 1 ) );
mono.exc = setdiff( exc, inh ); 
mono.inh = setdiff( inh, exc ); 
mono.both = intersect( exc, inh );
mono.unc = setdiff( ( 1: n )', [ exc; inh ] );

% post-prune
if ~prePrune && ~isempty( gidx )
    fidx = find( gidx );
    mono.pairsExc = mono.pairsExc( ismember( mono.pairsExc( :, 1 ), fidx ) & ismember( mono.pairsExc( :, 2 ), fidx ), : );
    mono.pairsInh = mono.pairsInh( ismember( mono.pairsInh( :, 1 ), fidx ) & ismember( mono.pairsInh( :, 2 ), fidx ), : );
    mono.pairsSync = mono.pairsSync( ismember( mono.pairsSync( :, 1 ), fidx ) & ismember( mono.pairsSync( :, 2 ), fidx ), : );
    mono.pairsDesync = mono.pairsDesync( ismember( mono.pairsDesync( :, 1 ), fidx ) & ismember( mono.pairsDesync( :, 2 ), fidx ), : );
    exc = unique( mono.pairsExc( :, 1 ) );
    inh = unique( mono.pairsInh( :, 1 ) );
    mono.exc = setdiff( exc, inh );
    mono.inh = setdiff( inh, exc );
    mono.both = intersect( exc, inh );
    mono.unc = setdiff( fidx, [ exc; inh ] );
end

mono.excited = unique( mono.pairsExc( :, 2 ) ); 
mono.inhibited = unique( mono.pairsInh( :, 2 ) );

% count instances for dual cells
cmat = zeros( length( mono.both ), 2 );
for i = 1 : length( mono.both )
    cmat( i, : ) = [ sum( mono.pairsExc( :, 1 ) == mono.both( i ) ) sum( mono.pairsInh( :, 1 ) == mono.both( i ) ) ];
end
mono.both = [ mono.both cmat ];

return
% EOF

% NOTE:
% check the pruning protocol: on 26nov11_1, 12-6 are exc-inh pair but 6
% does not appear at all in the pairsInh list!

% NOTE: that most exc+inh cells (except for 3.11) were actually multi-units
% with HFOs in the ACH as well (most were later removed by spikes_stats quality checks)

% However:
% (1) a given pair can be both inh and exc - but not in the same direction. if that
% happens at monosynaptic delays, then this probably means that the CCH has
% HFOs (e.g. es25nov11_4, 4.3x4.5 or many other pairs). From such CCHs one
% cannot determine mono-synaptic connectivity, so they should be pruned
% from BOTH lists
% a sync pair can also have exc
% a desync pair can also have inh

shankclu = determine_units( filebase, [], 'B' );
load( [ filebase '.sst' ], '-mat' );
shankclu0 = sst.shankclu;
gidx = ismember( shankclu0( :, 1 : 2 ), shankclu( :, 1 : 2 ), 'rows' );
mono = check_mono( filebase, gidx );

% consider using only different-shank pairs for the classification
% another issue - some of the pvalues are NEGATIVE!

poissinv( 1 - s2s.alpha / sum( s2s.t_ROI ), max( s2s.pred( s2s.t_ROI, :, : ), [], 1 ) ); 
poissinv( s2s.alpha / sum( s2s.t_ROI ), min( s2s.pred( s2s.t_ROI, :, : ), [], 1 ) );


i1 = 20; i2 = 18; 

figure
bh = bar( s2s.t, s2s.ccg( :, i1, i2 ) ); set( bh, 'EdgeColor', [ 0 0 1 ], 'FaceColor', [ 0 0 1 ] ); xlim( s2s.t( [ 1 end ] ) )
gHi = poissinv( 1 - s2s.alpha / sum( s2s.t_ROI ), max( s2s.pred( s2s.t_ROI, i1, i2 ), [], 1 ) );
gLo = poissinv( s2s.alpha / sum( s2s.t_ROI ), min( s2s.pred( s2s.t_ROI, i1, i2 ), [], 1 ) );

% gHi = poissinv( 1 - s2s.alpha / sum( s2s.t_ROI ), mm );
% gLo = poissinv( s2s.alpha / sum( s2s.t_ROI ), mm );
% 
alines( [ gLo gHi ], 'y', 'color', [ 1 0 0 ] );
alines( minmax( s2s.t( s2s.t_ROI ) ) + [ -0.5 0.5 ], 'x', 'color', [ 1 0 0 ] );


pvals = 1 - poisscdf( s2s.ccg( :, i1, i2 ) - 1, mm ) - poisspdf( s2s.ccg( :, i1, i2 ), mm ) * 0.5

% In principle:
% An excitatory cell: 
% -Is presynaptic in excitatory connections
% -May be postsynaptic in excitatory or inhibitory connections (although in
% CA1 it is most likely to be postsynaptic to inhibitory connections)
% -May be synchronized with other cells. This is possible due to
% common input to a pair of cells, and monosynaptic connectivity 
% between one of them and a third.
% 
% An inhibitory cell:
% -Is presynaptic in inhibitory connections
% -May be postsynaptic in excitatory or inhibitory connections (although in
% CA1 it is most likely to be postsynaptic to excitatory connections)
% -May be synchronized with other cells, e.g. due to gap junctions or
% common input
% -Unlikely to be desynchronized with other cells, although that is in
% principle possible (if there is precise excitation to one and
% simultaneous inhibition to the other; but I don't know what
% cellular/network mechanism can cause that)
%
% Single-pair issues:
% -a given CCH (i.e. a given pair of units) may only have one type of
% connectivity (if on the same side). Thus if a CCH has both inh and exc
% bins on the same side, this is most likley due to network events (e.g. HFOs), 
% poor isolation, or third party local circuits. 
%       -> The conservative solution is to treat the pair as unconnected; 
%       -> the more logical is to take the first (or larger) event and ignore the other
% -a given pair may have both sync and exc, but then it is not purely
% mono-synaptically presynaptic. It could be classified as sync or exc. 
%       -> One option is to ignore the excitatory connection and treat as only sync; 
%       -> the other is to classify according to the largest peak
% -a same-shank CCH cannot resemble the ACH of either unit too much. 
% If it does, this most likely indicates oversorting (e.g. bursting cells 
% treated as excitatory pairs, wide ACH treated as inhibitory pairs) and 
% the pair should not be treated as connected 
%       -> ignore the CCH
% 
% Algorithm:
% -go over CCHs, look only at the postsynaptic side. 
% -check for extreme bins in the monosynaptic range. This comprises the basic solutions. 
% -classify individual CCHs. still look only at postsynaptic side
%       -if same-shank, compare CCH to ACHs. if too similar (>0.9 CC), ignore
%       the CCH
%       -if only exc/inh/sync/desync, leave as is
%       -if two types (sync/exc, exc/inh, .. ), prune it: either use the
%       first (default), the larger, or remove from both classifications
%       ('conservative')

% this may still yield wierd cases due to marginal significance (e.g.
% m258r1_42, 3.22 x 2.21, which are both nice inhibitory interneurons
% also note e.g. 3.31 x 2.21, which is clearly driven (oscillatory CCHs,
% probably due to 3.31 bursting). 
% some "inhibitory" cells are only sync. with other inhibitory cells, and
% have no clear suppression of spiking. There is nothing much 

% e.g. 25nov11_3: 1.3 and 1.14, 1.22, 1.27, and 1.34 are probably oversorted
