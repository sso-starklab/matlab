% calc_asg          estimate ASG from spike trains
%
% call              [ g1, g2, act, sil, s, fig, cchmat ] = calc_asg( st1, st2 );
% 
% gets              st1             spike train, given as spike times (or a sparse vector)
%                   st2             second spike train
%
% does              (1) computes CCH, ACH1, ACH2, SEM, and other basics (or derives from the input)
%                   (2) deconvolves the ACH from the CCH
%                   (3) computes a predictor
%                   (4) computes difference CCH
%                   (5) computes rate ('gain') CCH
%                   (6) determines significance (in roiMS)
%                   (7) computes the ASG (in roiMS, regardless of the significance)
%                   (8) computes the SG50 and CA (for g1 only)
%                   (9) organizes results in a structure s
%                   (10) plots
%
% alternative input configurations are defined using an optional argument cmode: 
%
%                   cmode 1:        calc_asg( s2s, n12 )            assumes s2s has fields: ccg, shankclu, t, nspks, pred, filebase, alpha
%                   cmode 2:        calc_asg( cch, nspks1 )         receives the actual CCH and nspks1 (optional parameters: dt, ACHs)
%                   {cmode 3}:      calc_asg( st1, st2 )            receives the spike trains directly
%
%                   s2s             structure (see spikes2spikes)
%                   n12             [ idx1 idx2 ], indices of s2s
%                   cch             cross-correlation histogram
%                   nspks1          number of spikes in st1
%
% optional arguments (given as name/value pairs) include:
%
%                   cmode               {3}                 computational mode: 
%                                                               1   gets CCH, ACH, etc. from s2s structure
%                                                               2   gets CCH and nspks1 etc. from input arguments
%                                                               3   gets spike trains and computes everything internally
%                   dcmode              {1}                 de-convolution mode: 
%                                                               0   does not deconvolve
%                                                               1   deconvolves ACH1 then ACH2
%                                                              -1   deconvolves ACH1
%                                                              -2   deconvolves ACH2
%                                                              -3   deconvolves ACH1*ACH2
%                                                              -4   deconvolves ACH1 then ACH2
%                   jmode               {-1}                jittering mode: 
%                                                              -2   global mean
%                                                              -1   uses median filtering (non-linear)
%                                                               0   uses convolution (linear filtering)
%                                                               1   spike time jitter (equivalent to convolution)
%                                                               2   interval jitter
%                                                           jittering is supported only if (1) 'periods' is [] and (2) 'approxSEM' is 1
%                                                           (if periods required, remove before calling this routine)
%
%                   ach1                {[]}                ignored if cmode is 1 or 3
%                   ach2                {[]}                ignored if cmode is 1 or 3
%                   t                   {[]}                ignored if cmode is 3
%
%                   BinSizeMS           {1}         [ms]    CCH bin size
%                   halfWidthMS         {50}        [ms]    half-window for CCH calculations (ignored if cmode is 3)
%                   SpikesFs            {20000}     [Hz]    sampling rate (ignored if cmode is 1 or 2)
%                   deadTimeMS          {0}         [ms]    may be set to 0.4 for spike trains recorded on the same shank
%
%                   approxSEM           {1}                 approximate SEM if periods not given
%                   periods             {[]}                ignored if jmode is 1 or 2 (in that case, exclude the irrelevant spikes before calling this routine)
%                   correctPeriods      {1}                 correct periods (excludes trigger spike but not reference)
%
%                   CutOffMS            {0}         [ms]    burst filtering; should be set to 0 to enable proper devonvolution (ignored/warning if cmode is 1 or 2)
%                   dcmethod            {'fft'}             argument 6 to cch_deconv; 'wiener' also supported
%
%                   jitWindowMS         {5}         [ms]    half-window for jittering/filtering window
%                   convType            {'gauss'}   text    argument 3 to cch_conv ('triang' is equivalent to jittering both trains in a boxcar )
%                   hollowF             {[]}                argument 4 to cch_conv (1 means completely hollow)
%                   njreps              {1000}              argument to cch_jit
%
%                   roiMS               {[NaN 5]}   [ms]    should be causal; NaN in lower edge sets to first bin above zero lag
%                   asgMode             {1}                 0   detects global maxima in ROI (always exist, may be non-causal)
%                                                           1   detects local maxima in ROI (may not exist, but preserve causality)
%                   roiHardness         {[1 0]}             hard/soft edges (limits/does not limit ASG base by ROI)
%                   alpha0              {0.001}             only influences act/sil output arguments 
%
%                   graphics            {0}
%                   plottype            {'bar'}             'patch' plots lines with SEM bands; 'barlines' plots bars with SEM dashed lines
%
% returns           g1                  ASG of global positive extremum in the ROI
%                   g2                  ASG of global negative extremum in the ROI
%                   act                 f, 1 if any CCH bin in the ROI is sig. above the predictor
%                   sil                 flag, 1 if any CCH bin in the ROI is sig. below the predictor
%                   s                   structure with fields:
%                       cch0            [counts]    before deconvolution
%                       cch             [counts]    after deconvolution 
%                       sem             [counts]    SEM of the raw CCH
%                       t               [ms]        time vector
%                       ach1            [counts]    zero-lag bin set to zero
%                       ach2            [counts]    zero-lag bin set to zero
%                       sem1            [counts]    SEM of ACH1
%                       sem2            [counts]    SEM of ACH2
%                       nspks1                      total number of spikes in st1
%                       nspks2                      total number of spikes in st2
%                       pred                        jiitered cch
%                       dcch            [counts]    difference CCH 
%                       gcch            [spks/s]    rate ('gain') CCH 
%                       gsem            [spks/s]    SEM of the gain CCH 
%                       pvals                       p-values
%                       g1base          [samples]   temporal support for positive extremum 
%                       g2base          [samples]   temporal support for negative extremum
%                       theta           [SD]        [ theta0 theta1 ]; spiking threshold (SDs of membrane potential, standard normal PDF)
%                       sg50            []          [ sg50 CA ]
%                       asg             []          [ g1 g2 ]
%                       dt              [s]         bin size 
%                   fig                 handle to figure
%                   cchmat              cch matrix: one row (cch) for each relevant spike in st1
%
% calls             CCG                                         (blab)
%                   ParseArgPairs, replacetok                   (general)
%                   alines, patch_band                          (graph)
%                   inranges, sortranges                        (sets)
%                   cch_conv, cch_deconv, cch_jit, cch_mat      (spikes)
%                   calc_sem                                    (stats)
%
% see also          spikes2spikes

% 14-may-20 ES

% revisions
% 15-may-20 (1) expanded i/o to support various formats (partial)
% 31-may-20 (1) added dt to s output
% 08-jul-20 (1) added support for direct computations from spike trains (cmode 3)
%           (2) fixed ASG computation (nansum->nanmean)
%           (3) cmode 2 supported
% 02-oct-20 (1) implemented devconvolution of ACH from CCH
%           (2) modified nidx to support deadTimeMS of zero
% 03-oct-20 (1) implemented dual deconvolution of the ACHs from the CCH
%           (2) flow control, documentation, default arguments
%           (3) added blue vertical lines for g2
% 05-oct-20 (1) added computation and plotting of SEM bands (cmode 3)
%           (2) added computation of SEM bands to ACH (cmode 3)
%           (3) implemented periods (cmode 3)
% 08-oct-20 (1) added computation of SG50, CA
%           (2) applied kidx to ach1 and ach2 after deconvolution
%           (3) fix instead of ceil in make_cchmat
%           (4) extended support of cmode 1
% 09-oct-20 (1) theta0,1 defaulted to NaN
%           (2) added plots of ACH1, ACH2
%           (3) approxSEM added, defaulted to 1 (exact if BinSizeMS < CutOffMS)
% 10-oct-20 (1) fixed issue with convolution of ACH1*ACH2: first convolve,
%               then set zero-lag to 1 (delta function) - otherwise can
%               amplify a single frequency
%           (2) imposed dead time constraints (relevant for synthetic data)
%           (3) added roiHardness, defaulted to [ 1 0 ], preventing
%               a-causal ASGs
% 14-oct-20 (1) formulation of W simplified for all cmodes
%           (2) jmode added to compute predictor using interval jitter (cch_jit)
%           (3) deconvolution entirely by cch_deconv
%           (4) added ASG to the output structure
% 15-oct-20 (1) full support of jmode 1, 2
%           (2) excluded cch bins within deadTimeMS from pred estimation
%           (3) excluded nidx (deadTimeMS) in pred (influences significance estimation and ASG support)
%           (4) jmode -1 added (median filtering)
% 16-oct-20 (1) deconvolution of jittered CCH with jittered ACH
%           (2) modified default support (roiMS) to exclude zero-lag bin
% 17-oct-20 (1) fixed bug in roiMS( 1 )
%           (2) set nidx of gsem to NaN
%           (3) compute SEMs also for cmodes 1 and 2
%           (4) hollowF (cch_conv 4th argument) added 
%           (5) improved flow control and error checking (periods/cmode/dcmode/jmode)
%           (6) compute predictor also for cmodes 1 and 2
%           (7) make_cchmat moved to an external routine cch_mat
%           (8) njreps added (default, 1000)
%           (9) default values changed: dcmode=1; jmode=-1; CutOffMS=0
%           (10) W always odd
%           (11) pvals added to s
%           (12) median filtering done in cch_conv
%           (13) wiener filtering supported in cch_deconv and deconvfft
% 20-oct-20 (1) plottype added, defaulted to 'bar'
% 21-oct-20 (1) dcmode 3 added
% 29-oct-20 (1) dcmode 3 supported in argument check
%           (2) jmode -2 added and supported (global mean)
% 24-nov-20 (1) dcmode 4 added
%           (2) jmode 1/2 will not run for dcmode > 1
% 03-dec-20 (1) passed sequential deconvolution to cch_deconv
%           (2) nfft fixed for all deconvolution methods
% 06-dec-20 (1) jmode -2 modified to exclude ROI and a negative version of it
% 11-dec-20 (1) bug in ASG of negative peak corrected
% 12-dec-20 (1) modified computation of ASG s.t. only local extrema are
%               considered; this solves the issue of de/synchrony being
%               interpreted as negative/positive ASG
% 13-dec-20 (1) asgMode added, defaults to 1 (new behavior from 12-dec-20)
%           (2) asgMode x modified
%           (3) expansionFactor defined also for cmode 1
%           (4) full support of deconvolution also of jittered spike trains
%           (5) asgMode documented
% 14-dec-20 (1) asgMode 1 refined
% 16-dec-20 (1) BinSizeMS hard-limited by sampling rate
% 17-dec-20 (1) refinement of asgMode 1: case of consecutive identical values handled
% 12-jan-21 (1) added pvalue estimation for jmode -2
% 19-jan-21 (1) refinement of asgMode 1: if no local maximum, but the global maximum is at the far end of the
%               ROI (far from zero-lag), then use asgMode 0 (and the mirror for the asg2)
% 24-jan-21 (1) bug fix (asgMode 1)
% 03-mar-21 (1) bug fix (asgMode 1)

function [ g1, g2, act, sil, s, fig, cchmat ] = calc_asg( st1, st2, varargin )

%-------------------------------------------------------------------------
% 0. preparations

% constants
ARP                             = 0.002;                                    % [s], ARP for SG50 and CA computations
mfname                          = mfilename;

% initialize outputs
cchmat                          = [];

% arguments
[ BinSizeMS, halfWidthMS, jitWindowMS, SpikesFs ...
    , roiMS, roiHardness, convType, hollowF, alpha0, deadTimeMS, CutOffMS ...
    , approxSEM, correctPeriods ... 
    , ach1, ach2, t, periods, cmode, dcmode, dcmethod, jmode, njreps, asgMode ...
    , graphics, plottype ]      = ParseArgPairs( ...
    { 'BinSizeMS', 'halfWidthMS', 'jitWindowMS', 'SpikesFs' ...
    , 'roiMS', 'roiHardness', 'convType', 'hollowF', 'alpha0', 'deadTimeMS', 'CutOffMS' ...
    , 'approxSEM', 'correctPeriods' ... 
    , 'ach1', 'ach2', 't', 'periods', 'cmode', 'dcmode', 'dcmethod', 'jmode', 'njreps', 'asgMode' ...
    , 'graphics', 'plottype' } ...
    , { 1, 50, 5, 20000 ...
    , [ NaN 5 ], [ 1 0 ], 'gauss', [], 0.001, 0, 0 ...
    , 1, 1 ...
    , [], [], [], [], 3, 1, 'fft', -1, 1000, 1 ...
    , 0, 'bar' }...
    , varargin{ : } );

% check arguments 
if ~ismember( cmode, [ 1 2 3 ] )
    error( '%s: cmode %d not supported\n', mfname, cmode )
end
if ~ismember( dcmode, [ 0 1 -1 -2 -3 -4 ] )
    error( '%s: dcmode %d not supported\n', mfname, dcmode )
end
if ~ismember( jmode, [ -2 -1 0 1 2 ] )
    error( '%s: jmode %d not supported\n', mfname, jmode )
end
if jmode > 0 && dcmode > 1
    error( '%s: jmode %d not supported when dcmode is %d\n', mfname, jmode, dcmode )
end
if dcmode ~= 0
    if cmode == 3
        if dcmode == 1 || dcmode == -4                                      % two consecutive deconvolutions
            expansionFactor  	= 4;
        elseif dcmode == -1 || dcmode == -2 || dcmode == -3                 % single deconvolution step
            expansionFactor    	= 2;
        end
        halfWidthMS          	= expansionFactor * halfWidthMS;         	% expand edges for deconvolution
    else
        expansionFactor         = 1;
    end
end
if CutOffMS > 0 && dcmode ~= 0
    warning( '%s: Will not deconvolve fully when CutOff = %0.2g ms\n', mfname, CutOffMS )
end
if ~isempty( periods )
    periods                     = sortranges( periods );
    if ismember( cmode, [ 1 2 ] )
        error( '%s: Cannot apply periods to an existing CCH\n', mfname )
    end
    if ismember( jmode, [ 1 2 ] )
        error( '%s: Jittering with periods not supported; apply periods before calling this routine\n', mfname )
    end
end
if ismember( jmode, [ 1 2 ] ) && ismember( cmode, [ 1 2 ] )
    error( '%s: Cannot jitter without spike trains\n', mfname )
end
if ismember( jmode, [ 1 2 ] ) && approxSEM ~= 1
    error( '%s: The compbination of jittering and exact SEM estimation is not supported\n', mfname )
end
if ~ismember( plottype, { 'bar', 'patch', 'barlines' } )
    error( '%s: plottype %d not supported\n', mfname, plottype )
end
if ~ismember( asgMode, [ 0 1 ] )
    error( '%s: asgMode %d not supported\n', mfname, asgMode )
end
if BinSizeMS < ( 1000 / SpikesFs )
    error( '%s: bin size of %0.2g not supported when sampling rate is %0.2g', mfname, BinSizeMS, SpikesFs )
end

% derive W
W                               = 2 * ceil( jitWindowMS / BinSizeMS ) + 1;  % [samples]

% predictor name
switch jmode
    case -2
        pred_type               = 'global mean';
    case -1
        pred_type               = 'median filter';
        convType                = 'median';
    case 0
        pred_type               = 'convolution';
    case 1
        pred_type               = 'spike jitter';
    case 2
        pred_type               = 'interval jitter';
end

%-------------------------------------------------------------------------
% 1. compute cch, ach1, ach2 (and cchbins, nspks1, dt, alfa, t_ROI)
str                             = '';
switch cmode
    
    case 1 % use the data in s2s (cch, pred, ach)
        
        if ~isa( st1, 'struct' ) || ~isa( st2, 'numeric' )
            error( '%s: cmode %d requires first argument to be an s2s structure and second argument to be numeric', mfname, cmode )
        end
        s2s                     = st1;
        n12                     = st2;                                   	% [ shankclu1; shankclu2 ];
        n1                      = find( s2s.shankclu( :, 1 ) == n12( 1, 1 ) & s2s.shankclu( :, 2 ) == n12( 1, 2 ) );
        n2                      = find( s2s.shankclu( :, 1 ) == n12( 2, 1 ) & s2s.shankclu( :, 2 ) == n12( 2, 2 ) );
        
        t                       = s2s.t;                                    % [ms]
        nBins                   = ( length( t ) - 1 ) / 2; 
        
        cch                     = s2s.ccg( :, n1, n2 );
        ach1                    = s2s.ccg( :, n1, n1 );
        ach2                    = s2s.ccg( :, n2, n2 );
        nspks1                  = s2s.nspks( n1 );
        nspks2                  = s2s.nspks( n2 );
        
        [ ~, basename ]         = fileparts( s2s.filebase );
        str                     = sprintf( '%s: %d.%d x %d.%d; ' ...
            , basename, n12( 1, 1 ), n12( 1, 2 ) ...
            , n12( 2, 1 ), n12( 2, 2 ) );

        % compute SEM
        if approxSEM
            sem                 = sqrt( cch );
            sem1                = sqrt( ach1 );
            sem2                = sqrt( ach2 );
        else
            sem                 = zeros( size( cch ) );
            sem1                = zeros( size( ach1 ) );
            sem2                = zeros( size( ach2 ) );
        end
        
    case 2 % use the cch as given (and if given, t and ach)
        
        % check inputs
        if ~isa( st1, 'numeric' ) || ~isa( st2, 'numeric' )
            error( '%s: cmode %d requires first and second arguments to be numeric', mfname, cmode )
        end
        cch                     = st1;
        nspks1                  = st2;
        nBins                   = ( length( cch ) - 1 ) / 2; 
        if isempty( t )                                                  	% assume default bin size (e.g. 1 ms)
            t                   = ( -nBins : nBins )' * BinSizeMS;
        else
            BinSizeMS         	= diff( t( 1 : 2 ) );                       % assume fixed bin size
        end
        if isempty( ach1 )
            ach1                = NaN( length( cch ), 1 );
            if dcmode ~= 0
                error( '%s: Cannot deconvolve ACH1 if ACH1 not given\n', mfname )
            end
        end
        if isempty( ach2 )
            ach2                = NaN( length( cch ), 1 );
            if ismember( dcmode, [ -4 -3 1 ] )
                error( '%s: Cannot deconvolve both ACH1 and ACH2 if ACH2 not given\n', mfname )
            end
        end
        nspks2                  = NaN;
        
        % compute SEM
        if approxSEM
            sem                 = sqrt( cch );
            sem1                = sqrt( ach1 );
            sem2                = sqrt( ach2 );
        else
            sem                 = zeros( size( cch ) );
            sem1                = zeros( size( ach1 ) );
            sem2                = zeros( size( ach2 ) );
        end
        
    case 3 % compute raw cch from spike trains (diluted according to periods)
        
        % check inputs
        if ~isa( st1, 'numeric' ) || ~isa( st2, 'numeric' )
            error( '%s: cmode %d requires first and second arguments to be numeric', mfname, cmode )
        end
        if issparse( st1 ) || issparse( st2 )
            if size( st1( : ) ) ~= size( st2( : ) )
                error( 'input size mistmatch: if st1 or st2 is sparse, both must have the same dimensions' )
            end
            st1 = find( st1( : ) );
            st2 = find( st2( : ) );
        end
        
        % get some basic parameters
        nBins                   = halfWidthMS / BinSizeMS;               	% number of bins on each side of the CCH
        BinSize                 = round( BinSizeMS * SpikesFs / 1000 );   	% [samples]
        CutOff                  = round( CutOffMS * SpikesFs / 1000 );   	% [samples]
        t                       = ( -nBins : nBins )' * BinSizeMS;       	% [ms]
        
        % impose dead time (detection) constraints
        if deadTimeMS > 0
            deadTime            = round( deadTimeMS * SpikesFs / 1000 );    % [samples]
            st                  = [ [ st1 ones( size( st1, 1 ), 1 ) ]; [ st2 2 * ones( size( st2, 1 ), 1 ) ] ];
            st                  = sortrows( st, 1 );
            isis                = diff( st( :, 1 ) );
            ridx                = find( isis <= deadTime ) + 1;
            st( ridx, : )       = [];
            idx1                = st( :, 2 ) == 1;
            st1                 = st( idx1, 1 );
            st2                 = st( ~idx1, 1 );
        end
        
        % burst-filter each train
        if CutOff > 0
            isis1               = diff( st1 ); 
            isis2               = diff( st2 ); 
            idx1                = find( isis1 <= CutOff ) + 1; 
            idx2                = find( isis2 <= CutOff ) + 1; 
            st1( idx1 )         = [];
            st2( idx2 )         = [];
        end
        
        % compute the CCH and the ACH
        nspks1                  = length( st1 );
        nspks2                  = length( st2 );
        res                     = [ st1; st2 ];
        clu                     = [ ones( nspks1, 1 ); 2 * ones( nspks2, 1 ) ];
        [ ~, sidx ]             = sort( res );
        clu                     = clu( sidx );
        res                     = res( sidx );
        [ ccg, ~, pairs ]       = CCG( res, clu, BinSize, nBins, SpikesFs, [ 1 2 ], 'count' );
        
        % attend to periods and compute SEM
        if isempty( periods ) && ( approxSEM || CutOff > BinSize )
            
            % approximate SEM by square root of the count
            % (is exact if all cchmat counts are 0/1, i.e. if CutOff > BinSize )
            cch                 = ccg( :, 1, 2 );
            ach1                = ccg( :, 1, 1 );
            ach2                = ccg( :, 2, 2 );
            sem                 = sqrt( cch );
            sem1                = sqrt( ach1 );
            sem2                = sqrt( ach2 );
            
            if jmode > 0
                [ pvals, ~, mcch, ~, mach1, mach2 ]  = cch_jit( st1, st2 ...
                    , 'BinSizeMS', BinSizeMS, 'halfWidthMS', halfWidthMS ...
                    , 'SpikesFs', SpikesFs, 'jitWindowMS', jitWindowMS ...
                    , 'jmode', jmode, 'nreps', njreps );
            end
            
        else
            
            % compute an exact SEM by averaging over all in-period trigger spikes
            
            % compute CCH matrices (one row for each relevant spike)
            [ cchmat, ust1c ]   = cch_mat( clu, res, BinSize, nBins, pairs, [ 1 2 ] );
            [ achmat1, ust1 ]   = cch_mat( clu, res, BinSize, nBins, pairs, [ 1 1 ] );
            [ achmat2, ust2 ]   = cch_mat( clu, res, BinSize, nBins, pairs, [ 2 2 ] );
            
            % keep only trigger spikes in periods (referred spikes may be out of the periods)
            if ~isempty( periods )
                
                if correctPeriods                                           % reduce periods by halfwin from each side (then, referred spikes will also be in the periods)
                    osp         = ones( size( periods, 1 ), 1 ) * [ 1 -1 ] * nBins;
                    periods     = periods + osp;
                end
                
                idx1            = inranges( st1, periods );
                idx2            = inranges( st2, periods );
                st1p            = st1( idx1 );
                st2p            = st2( idx2 );
                nspks1          = length( st1p );
                nspks2          = length( st2p );
                [ ~, i1 ]       = intersect( ust1c, st1p );
                cchmat          = cchmat( i1, : );
                [ ~, i1 ]       = intersect( ust1, st1p );
                achmat1         = achmat1( i1, : );
                [ ~, i1 ]       = intersect( ust2, st2p );
                achmat2         = achmat2( i1, : );
            end
            
            % compute CCH, ACH, and SEM by averaging over all relevant trigger spikes
            cch                 = sum( cchmat, 1 )';
            ach1                = sum( achmat1, 1 )';
            ach2                = sum( achmat2, 1 )';
            sem                 = calc_sem( cchmat, 1 )' * size( cchmat, 1 );
            sem1                = calc_sem( achmat1, 1 )' * size( achmat1, 1 );
            sem2                = calc_sem( achmat2, 1 )' * size( achmat2, 1 );
            
        end

end

if nBins ~= round( nBins )
    error( '%s: Check CCH length (should be odd)', mfname )
end
if ~isequal( size( cch ), size( ach1 ), size( ach2 ), size( t ) )
    error( '%s: Check correspondence between cch, ach, and t lengths (should be the same)', mfname )
end

dt                              = diff( t( 1 : 2 ) ) / 1000;              	% [s]
n                               = 2 * nBins + 1;
cch0                            = cch;

%-------------------------------------------------------------------------
% 2. deconvolve the ACH(s) from the CCH
if dcmode ~= 0
    
    % deconvolve
    nfft                        = length( cch );
    switch dcmode 
        case 1
            [ cch, bn, kidx ]   = cch_deconv( cch, ach1, nspks1, ach2, nspks2, dcmethod, nfft );
            cch0                = cch0( kidx, : );
            sem                 = sem( kidx, : );
            ach1               	= ach1( kidx, : );
            ach2             	= ach2( kidx, : );
            sem1               	= sem1( kidx, : );
            sem2              	= sem2( kidx, : );
        case -1
            [ cch, bn, kidx ]   = cch_deconv( cch, ach1, nspks1, [], [], dcmethod, nfft );
        case -2
            [ cch, bn, kidx ]   = cch_deconv( cch, ach2, nspks2, [], [], dcmethod, nfft );
        case -3
            [ cch, bn, kidx ]   = cch_deconv( cch, ach1, nspks1, ach2, nspks2, dcmethod, nfft, 3 );
        case -4
            [ ccht, ~, kidxt ]  = cch_deconv( cch, ach1, nspks1, [], [], dcmethod, nfft );
            cch0                = cch0( kidxt, : );
            sem                 = sem( kidxt, : );
            ach1               	= ach1( kidxt, : );
            ach2             	= ach2( kidxt, : );
            sem1               	= sem1( kidxt, : );
            sem2              	= sem2( kidxt, : );
            [ ccht, bn, kidx ] 	= cch_deconv( flipud( ccht ), ach2, nspks2, [], [], dcmethod, nfft );
            cch                 = flipud( ccht );
    end
    if jmode > 0                                                            % deconvolve also the jitter-based cch
        switch dcmode
            case 1
                mcch            = cch_deconv( mcch, mach1, nspks1, mach2, nspks2, dcmethod, nfft );
            case -1
                mcch            = cch_deconv( mcch, mach1, nspks1, [], [], dcmethod, nfft );
            case -2
                mcch            = cch_deconv( mcch, mach2, nspks2, [], [], dcmethod, nfft );
            case -3
                mcch            = cch_deconv( mcch, mach1, nspks1, mach2, nspks2, dcmethod, nfft, 3 );
        end
    end

    % assign
    t                           = bn * BinSizeMS;
    n                           = 2 * floor( nBins / expansionFactor ) + 1;
    if ismember( dcmode, [ -1 -2 -3 -4 ] )
        cch0                 	= cch0( kidx, : );
        sem                   	= sem( kidx, : );
        ach1                  	= ach1( kidx, : );
        ach2                	= ach2( kidx, : );
        sem1                 	= sem1( kidx, : );
        sem2                  	= sem2( kidx, : );
    end

end

%-------------------------------------------------------------------------
% 3. compute the predictor
if deadTimeMS ~= 0
    nidx                        = t <= deadTimeMS & t >= -deadTimeMS;       % exclude these bins from all subsequent steps
else
    nidx                        = false( length( t ), 1 );
end
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % optional: causality imposed
end
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );        % causality not necessarily imposed
t_ROI                           = t_ROI & ~nidx;                            % dead time imposed
ft_ROI                          = find( t_ROI );
if cmode == 1                                                               % get the alpha level from s2s
    alfa                        = s2s.alpha;
else
    alfa                        = alpha0;
end
if ( W * 3 + 2 ) * sqrt( 2 ) > length( cch )
    error( '%s: Mismatch between CCH length (%d) and jitWindowMS (%d): elongate CCH or shorten window' ...
        , mfname, length( cch ), jitWindowMS )
end
if ismember( jmode, [ 1 2 ] )
    pred                        = mean( mcch, 2 );                        	% jitter-based predictor (ignores dead time, refractoriness, etc.)
elseif jmode == -2
    eval                        = max( abs( t( ft_ROI ) ) );
    tidx                        = t < -eval | t > eval;
    pred                        = ones( size( cch, 1 ), 1 ) * mean( cch( tidx, : ), 1 );
    pvals                       = 1 - poisscdf( cch - 1, pred ) - poisspdf( cch, pred ) * 0.5;
else
    if sum( nidx )
        cchN                    = cch( ~nidx, : );                          % ignore under-estimated bins d.t. dead time
        [ pvalsN, predN ]       = cch_conv( cchN, W, convType, hollowF );   % linear/non-linear filtering of the CCH
        pred                    = NaN( size( cch ) );
        pred( ~nidx, : )        = predN;
        pvals                   = NaN( size( cch ) );
        pvals( ~nidx, : )       = pvalsN;
    else
        [ pvals, pred ]         = cch_conv( cch, W, convType, hollowF );   % linear/non-linear filtering of the CCH
    end
end
pred( nidx, : )                 = NaN;
pvals( nidx, : )                = NaN;

%-------------------------------------------------------------------------
% 4. compute difference cch
dcch                            = cch - pred;                               % [counts]

%-------------------------------------------------------------------------
% 5. compute rate ('gain') cch
gcch                            = dcch / ( nspks1 * dt );                   % this will always be assymetric
gsem                            = sem / ( nspks1 * dt );                    % [spks/s]
gsem( nidx, : )                 = NaN;

%-------------------------------------------------------------------------
% 6. determine if any bin in t_ROI is significant
nBonf                           = sum( t_ROI );                             % dead-time considered
gbUpper                         = poissinv( 1 - alfa / nBonf, max( pred( t_ROI), [], 1 ) );
gbLower                         = poissinv( alfa / nBonf, min( pred( t_ROI ), [], 1 ) );
act                             = false( 1 );
sil                             = false( 1 );
if any( cch( t_ROI ) >= gbUpper )
    act                         = 1;
end
if any( cch( t_ROI ) <= gbLower )
    sil                         = 1;
end

%-------------------------------------------------------------------------
% 7. compute the ASG

% find the largest local extremum in the ROI
% then find the first zero crossings before and after that extremum 
% compute the intergral (can be negative or positive): mean divided by span
% 
% do this:
% -for both the maxima (positive ASG) and for the minima (negative)
% -for the extrema in the ROI, regardless of significance (step 6)

% support of positive extremum
switch asgMode
    case 0
        x                       = cch( t_ROI );
        [ ~, maxidx ]           = max( x );
    case 1
        pidx                    = ft_ROI( [ 1 end ] ) + [ -1 1 ]';
        xidx                    = [ pidx( 1 ); ft_ROI; pidx( 2 ) ];
        nROI                    = length( ft_ROI ) + 1;
        x                       = cch( xidx );
        sx                      = find( diff( x ) == 0 );                   % same values of x
        ux                      = x;
        ix                      = ( 1 : length( x ) )';
        ix( sx )                = [];
        ux( sx )                = [];
        tmpidx                  = find( diff( sign( diff( ux ) ) ) < -1 ) + 1;
        maxidx                  = ix( tmpidx ) - 1;
        if ~isempty( sx ) && sx( end ) == nROI && length( ux ) > 1 && ux( end ) > ux( end - 1 )
            maxidx              = sx( end ) - 1;
        end
end
if length( maxidx ) > 1
    [ ~, subidx ]               = max( x( maxidx ) );
    maxidx                      = maxidx( subidx );
end
if isempty( maxidx )
    sidx                        = [];
else
    pidx                    	= ft_ROI( maxidx );
    si                        	= 1;                                    	% first preceding negative value
    for i                     	= pidx : -1 : 1
        if gcch( i ) < 0
            si                 	= i + 1;
            break
        end
    end
    ei                         	= n;                                    	% first proceeding negative value
    for i                      	= pidx : n
        if gcch( i ) < 0
            ei                 	= i - 1;
            break
        end
    end
    sidx                       	= si : ei;
end
if roiHardness( 1 )
    rmv                         = sidx < ft_ROI( 1 );
    sidx( rmv )                 = [];
end
if roiHardness( 2 )
    rmv                         = sidx > ft_ROI( end );
    sidx( rmv )                 = [];
end
g1base                          = sidx;
% compute the integral
a1                              = nanmean( gcch( sidx ) );
b1                              = sum( ~isnan( gcch( sidx ) ) ) * dt;
g1                              = a1 * b1;

% support of negative extremum
switch asgMode
    case 0
        [ ~, minidx ]           = min( x );
    case 1
        tmpidx                  = find( diff( sign( diff( ux ) ) ) > 1 ) + 1;
        minidx                  = ix( tmpidx ) - 1;
        if ~isempty( sx ) && sx( end ) == nROI && length( ux ) > 1 && ux( end ) < ux( end - 1 )
            minidx              = sx( end ) - 1;
        end
end
if length( minidx ) > 1
    [ ~, subidx ]               = min( x( minidx ) );
    minidx                      = minidx( subidx );
end
if isempty( minidx )
    sidx                        = [];
else
    pidx                     	= ft_ROI( minidx );
    si                         	= 1;                                       	% first preceding positive value
    for i                      	= pidx : -1 : 1
        if gcch( i ) > 0
            si                	= i + 1;
            break
        end
    end
    ei                       	= n;                                    	% first proceeding positive value
    for i                      	= pidx : n
        if gcch( i ) > 0
            ei               	= i - 1;
            break
        end
    end
    sidx                     	= si : ei;
end
if roiHardness( 1 )
    rmv                         = sidx < ft_ROI( 1 );
    sidx( rmv )                 = [];
end
if roiHardness( 2 )
    rmv                         = sidx > ft_ROI( end );
    sidx( rmv )                 = [];
end
g2base                          = sidx;
% compute the integral
a2                              = nanmean( gcch( sidx ) );
b2                              = sum( ~isnan( gcch( sidx ) ) ) * dt;
g2                              = a2 * b2;

%-------------------------------------------------------------------------
% 8. compute SG50
if isempty( g1base )
    theta0                      = NaN;
    theta1                      = NaN;
else
    lambda1                     = max( cch( g1base ) ) / ( nspks1 * dt );
    lambda0                     = mean( pred( g1base ) ) / ( nspks1 * dt );
    theta0                      = -norminv( lambda0 * ARP );                % depends on ARP
    theta1                      = -norminv( lambda1 * ARP );
end
dtheta                          = theta0 - theta1;
sg50                            = dtheta / theta0;
ca                              = sg50 / ( 2 * g1 );                        % coincidence advantage

%-------------------------------------------------------------------------
% 9. organize results in a structure
s.cch0                          = cch0;                                     % before deconvolution, [counts]
s.cch                           = cch;                                      % after deconvolution, [counts]
s.sem                           = sem;                                      % SEM of the raw CCH, [counts]
s.cchbins                       = t;                                        % time vector, [ms]
s.ach1                          = ach1;                                     % zero-lag bin set to zero, [counts]
s.ach2                          = ach2;                                     % zero-lag bin set to zero, [counts]
s.sem1                          = sem1;                                     % SEM of ACH1, [counts]
s.sem2                          = sem2;                                     % SEM of ACH2, [counts]
s.nspks1                        = nspks1;                                   % total number of spikes in st1
s.nspks2                        = nspks2;                                   % total number of spikes in st2
s.pred                          = pred;                                     % jiitered cch
s.dcch                          = dcch;                                     % cch - pred, [counts]
s.gcch                          = gcch;                                     % dcch / ( dt * nspks1 ) [spks/s]
s.gsem                          = gsem;                                     % SEM of the gain CCH, [spks/s]
s.pvals                         = pvals;                                    % p-values
s.g1base                        = g1base;                                   % support for positive extremum [samples]
s.g2base                        = g2base;                                   % support for negative extremum [samples]
s.theta                         = [ theta0 theta1 ];                        % theta (SD of standard normal distribution)
s.asg                           = [ g1 g2 ];                                % ASG 
s.sg50                          = [ sg50 ca ];                              % SG50 and CA
s.dt                            = dt;                                       % bin size [s]

%-------------------------------------------------------------------------
% 10. plot if requested
if ~graphics
    fig                         = NaN;
    return
end
fig                             = figure;
for spn                         = 1 : 6
    subplot( 2, 3, spn )
    switch spn
        case 1
            switch plottype
                case 'bar'
                    bar( s.cchbins, s.ach1, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                case 'barlines'
                    bar( s.cchbins, s.ach1, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                    line( s.cchbins, s.ach1 + s.sem1, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    line( s.cchbins, s.ach1 - s.sem1, 'color', [ 0 0 0 ], 'linestyle', '--' );
                case 'patch'
                    patch_band( s.cchbins, s.ach1, s.sem1, [ 0 0 0 ], [ 1 1 1 ] * 0.7 );
            end
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            ylabel( 'Counts/bin' )
            title( 'Raw ACH1' )
        case 2
            switch plottype
                case 'bar'
                    bar( s.cchbins, s.ach2, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                case 'barlines'
                    bar( s.cchbins, s.ach2, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                    line( s.cchbins, s.ach2 + s.sem2, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    line( s.cchbins, s.ach2 - s.sem2, 'color', [ 0 0 0 ], 'linestyle', '--' );
                case 'patch'
                    patch_band( s.cchbins, s.ach2, s.sem2, [ 0 0 0 ], [ 1 1 1 ] * 0.7 );
            end
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            ylabel( 'Counts/bin' )
            title( 'Raw ACH2' )
        case 3
            switch plottype
                case 'bar'
                    bar( s.cchbins, s.cch0, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                case 'barlines'
                    bar( s.cchbins, s.cch0, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                    line( s.cchbins, s.cch0 + s.sem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    line( s.cchbins, s.cch0 - s.sem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    ylim( [ 0 max( ylim ) ] )
                case 'patch'
                    patch_band( s.cchbins, s.cch0, s.sem, [ 0 0 0 ], [ 1 1 1 ] * 0.7 );
            end
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            ylabel( 'Counts/bin' )
            title( 'Raw CCH' )
        case 4
            switch plottype
                case 'bar'
                    bar( s.cchbins, s.cch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                case 'barlines'
                    bar( s.cchbins, s.cch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                    line( s.cchbins, s.cch + s.sem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    line( s.cchbins, s.cch - s.sem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    ylim( [ 0 max( ylim ) ] )
                case 'patch'
                    patch_band( s.cchbins, s.cch, s.sem, [ 0 0 0 ], [ 1 1 1 ] * 0.7 );
            end
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            line( s.cchbins, s.pred, 'linewidth', 2, 'color', [ 1 0 0 ] )
            ylabel( 'Counts/bin' )
            if dcmode ~= 0
                title( sprintf( 'Deconvolved CCH (%d) with %s predictor (%0.2g ms)', dcmode, pred_type, jitWindowMS ) )
            else
                title( sprintf( 'CCH with %s predictor (%0.2g ms)', pred_type, jitWindowMS ) )
            end
        case 5
            switch plottype
                case 'bar'
                    bar( s.cchbins, s.dcch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                case 'barlines'
                    bar( s.cchbins, s.dcch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                    line( s.cchbins, s.dcch + s.sem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    line( s.cchbins, s.dcch - s.sem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    ylim( [ 0 max( ylim ) ] )
                case 'patch'
                    patch_band( s.cchbins, s.dcch, s.sem, [ 0 0 0 ], [ 1 1 1 ] * 0.7 );
            end
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
            ylabel( 'Counts/bin' )
            title( 'Difference CCH' )
        case 6
            switch plottype
                case 'bar'
                    bar( s.cchbins, s.gcch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                case 'barlines'
                    bar( s.cchbins, s.gcch, 1, 'facecolor', 'k', 'edgecolor', 'k' );
                    line( s.cchbins, s.gcch + s.gsem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                    line( s.cchbins, s.gcch - s.gsem, 'color', [ 0 0 0 ], 'linestyle', '--' );
                case 'patch'
                    patch_band( s.cchbins, s.gcch, s.gsem, [ 0 0 0 ], [ 1 1 1 ] * 0.7 );
            end
            ylabel( 'Conditional rate [spikes/s]' )
            alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
            alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
            xlims               = s.cchbins( minmax( s.g1base ) )';
            if length( xlims ) == 2
                alines( xlims + [ -0.5 0.5 ] * dt * 1000, 'x', 'linestyle', '--', 'color', [ 1 0 0 ] );
            end
            xlims               = s.cchbins( minmax( s.g2base ) )';
            if length( xlims ) == 2
                alines( xlims + [ -0.5 0.5 ] * dt * 1000, 'x', 'linestyle', '--', 'color', [ 0 0 1 ] );
            end
            tstr                = sprintf( 'Rate CCH; %sASG=%0.3g, %0.3g', str, g1, g2 );
            title( replacetok( tstr, '\_', '_' ) )
    end
    set( gca, 'box', 'off', 'tickdir', 'out' )
end

return

% EOF

%------------------------------------------------------------
% example (same shank PYR-INT (mono))
filebase                        = filebaseLookup( 'mC41', -33 );
shankclu1                       = [ 3 9 ];
shankclu2                       = [ 3 7 ];

%------------------------------------------------------------
% mode #1 (from the s2s file)


% prepare for call
load( [ filebase '.s2s' ], '-mat', 's2s' )
n12                             = [ shankclu1; shankclu2 ];

% call 
[ g1, g2, act, sil, s ]         = calc_asg( s2s, n12, 'graphics', 1, 'cmode', 1 );                      % w/ deconvolution

%------------------------------------------------------------
% mode #2 (from the ach and/or the cch)

% prepare for call
n1                              = find( s2s.shankclu( :, 1 ) == shankclu1( 1 ) & s2s.shankclu( :, 2 ) == shankclu1( 2 ) );
n2                              = find( s2s.shankclu( :, 1 ) == shankclu2( 1 ) & s2s.shankclu( :, 2 ) == shankclu2( 2 ) );
cch                             = s2s.ccg( :, n1, n2 );
ach1                            = s2s.ccg( :, n1, n1 );
nspks1                          = s2s.nspks( n1 );

% call 
[ g1, g2, act, sil, s ]         = calc_asg( cch, nspks1, 'cmode', 2, 'graphics', 1, 'dcmode', 0 );      % w/o deconvolution
[ g1, g2, act, sil, s ]         = calc_asg( cch, nspks1, 'ach1', ach1, 'cmode', 2, 'graphics', 1 );     % w/ ACH and deconvolution

%------------------------------------------------------------
% mode #3 (from the raw spike trains):

% load all spikes (of B and better units)
spk                             = load_spikes( filebase );
par                             = LoadXml( filebase );
SpikesFs                        = par.SampleRate;

% remove spikes during stimuli
vals                            = LoadStims( filebase );
uvals                           = uniteranges( vals( :, 1 : 2 ) );
ridx                            = inranges( spk.res, uvals );
spk.clu( ridx )                 = [];
spk.res( ridx )                 = [];

% keep only the spike times of the two relevant units
clunum1                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu1, 'rows' ), 1 );
clunum2                         = spk.map( ismember( spk.map( :, 2 : 3 ), shankclu2, 'rows' ), 1 );
st1                             = spk.res( spk.clu == clunum1 );
st2                             = spk.res( spk.clu == clunum2 );

% call 
[ g1, g2, act, sil, s ]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 1, 'dcmode', 0 );     % w/o deconvolution
[ g1, g2, act, sil, s ]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 1 );                  % w/ deconvolution
[ g1, g2, act, sil, s ]         = calc_asg( st1, st2, 'cmode', 3, 'graphics', 1, 'BinSizeMS', 0.1 );% w/ deconvolution and a smaller bin size
