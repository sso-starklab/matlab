% CALC_SPATIAL_WAVEFORM_FEATURES        bi-polarity and wavefront
%
% call                  [ vA, vT, bpi, tV, vB, pA, nA, pT, nT, upol ] = calc_spatial_waveform_features( w, TH, spkFs, graphics )
%
% gets                  w           waveform over multiple channels of a single unit, [ nsites x nsamples ]
%                       TH          {40}; threshold of SNR, micro-volts 
%                                            could be a matrix of sd over multiple channels, same size as w
%                       spkFs       {20000}; sampling rate [Hz]
%                       graphics    {0}; plot
%                       THflag      {0}; if 0, thank TH is calculated as
%                                        2SDs to the BL. if 1, 1SD to a point N samples to the 
%                                        left for the positive peak. 
%
% returns               vA          [uV] vector of ( nsites x 1 ), signed amplitude of extremum on each site
%                                       with the condition of peak above TH 
%                       vT          [ms] same, Time lag of local extermum relative to global extremum
%                                       with the condition of peak above TH 
%                       bpi         global BPI, the vB of the extremal channel, 
%                                       the closer to 1, the more symmetric ('bipolar')
%                       tV          [ms] SD of vT
%                       vB          vector of ( nsites x 1), channel level BPI, calculated as: 
%                                       ( p - abs( n ) ) / ( p + abs( n ) )
%                                   BPI is scaled from -1 to 1 where:
%                                       0     := symmetric (perfectly bipolar)
%                                       1     := perfect Pspike
%                                       -1    := perfect Nspike
%                                     with the condition of peak appearing before trough
%                       pA          [uV] vector of ( nsites x 1 ), peak value on each site
%                       nA          [uV] vector of ( nsites x 1 ), trough value on each site
%                       pT          [ms] same, Time lag of peak relative to global extremum
%                       nT          [ms] same, Time lag of trough relative to global extremum
%                       upol        an index which determine whether all valid extrema are 
%                                       homogenous or heterogenous 
%                                       1     := Pspikes                       (homogenous)
%                                       -1    := Nspikes                       (homogenous)
%                                       0     := some Pspikes and some Nspikes (heterogenous)
%
% calls                 alines, plot_spk_waveforms, imupsample, myjet, calibration
%                       local_max

% 08-dec-19 SSo & ES

% revisions
% 27-dec-19 set vT to NaN for uninformative amplitudes (vA==0)
% 23-feb-21 (1) extended vA and vT to {pA,nA} and {pT,nT}
%           (2) modified graphics to include {pA,nA} and {pT,nT}
%           (3) changed BPI to be channel-based, and include a condition on
%                   the sequence (positive then negative)
%           (4) added channel level BPI (vB)
%           (5) compute tV (SD of vT) - unit level metric
% 01-mar-21 updated documentation
% 30-mar-21 commented
% 06-apr-21 modified BP definition to rely on local extrema rather than
%               global extrema 
% 07-apr-21 (1) added delta parameter: "significant" difference between amplitudes of different samples
%           (2) added BPI conditions: global maximum is larger by delta from last sample
%                                     global minimum is smaller by delta from first sample
% 23-jun-21 (1) removed delta; changed TH to 40 uV (same as in check_cluster_quality)
% 28-jun-21 (1) changed everything completely: the criteria of a �Bipolar unit�: 
%                                                   (1) one positive and one negative local extremum per channel
%                                                   (2) both extrema are detectable (above threshold of 40 ?V)
%                                                   (3) positive extremum precedes negative extremum
%           (2) including the graphics
%           (3) added upol for determining whether all valid extrema are 
%               homogenous or heterogenous (Punit, Nunit, Dunit)
% 12-sep-21 (1) fix graphics to support m>10
%           (2) changed vB to iclude data for all units (also for units that do not answer
%               the conditions of Bipolar units)
% 13-sep-21 changed graphics to include CSD data 
% 15-feb-22 added recieving TH as a matrix of sd per channel (e.g. sst.sd)
% 17-feb-22 conseptual changes:
%         (1) We first corrected every mean spike according to a baseline calculated as you do in wfeatures function (upsample and mean).
%         (2)Then the definition is:
%             existence of a positive peak
%             the positive peak happens before the negative peak
%             the corrected negative peak is smaller by 2sd (according to sst.sd) than the baseline 
%             a) calculating the number of samples between the positive and the negative peaks (pre-index) 
%             b) the corrected positive peak would be 1sd larger than the value on the pre-index sample before
% 19-jun-22 (1) actually using the TH kmat correctly
%           (2) adding abvTH logical flag for units with absolute value of at least one channel above TH

% to do (28-jun-21):
% (1) update documentation (help)
% (2) consider extending BPI to AUC (compute pA and nA using AUC)

function [ vA, vT, bpi, tV, vB, pA, nA, pT, nT, upol, tp, rA, abvTH ] = calc_spatial_waveform_features( w, TH, spkFs, graphics,varargin)

[ THflag ] = ParseArgPairs( ...
    { 'THflag'} ...
    , { 0}...
    , varargin{ : } );

nargs                           = nargin;
if nargs < 1 || isempty( w )
    return
end
if nargs < 2 || isempty( TH )
    TH                          = 40;                                       % same units as w (e.g. micro-V)
end
if nargs < 3 || isempty( spkFs )
    spkFs                       = 20000;
end
if nargs < 4 || isempty( graphics )
    graphics                    = 0;
end

if ~ismatrix( w )
    return
end
[ m, n ]                        = size( w );
abvTH = false;
if size( TH,1 ) >1 
    [ msd, nsd ]                        = size( TH );
    if ~isequal([m,n],[msd,nsd])
        error('TH and w not the same size')
    end        
end
%---------------------------------------------------------------
% general channel specific properties (nA, pA, vA, nT, pT, vT)
%---------------------------------------------------------------
tp = NaN;
% upsample waveforms
w2                           = fft_upsample( w', 4, 1 );

% compute
[ nn, mm ]                    = size( w2 );                                % number of samples for baseline estimation
bvals                       = mean( w2( 1 : floor( nn / 3 ), : ), 1 );    % baseline

mean_amp = bvals';
% detect all local extrema
wt                              = (w-mean_amp)'; 
[ eidx, evals, etype ]       	= local_max( wt, 'ext' );
emat                            = [ eidx evals etype ];

% get the peak and trough values for each channel (value, sign, and lag)
% based on local extrema (in time) only:
nA                              = NaN( m, 1 );
pA                              = NaN( m, 1 );
nidx                            = NaN( m, 1 );
pidx                            = NaN( m, 1 );
for i = 1 : m
    imat                        = emat( emat( :, 2 ) == i, : );
    [ pvali, pidxi ]            = max( imat( :, 3 ) );
    pidx( i )                   = imat( pidxi, 1 );
    pA( i )                     = pvali;
    [ nvali, nidxi ]            = min( imat( :, 3 ) );
    nidx( i )                   = imat( nidxi, 1 );
    nA( i )                     = nvali;
end

% determine overall extremum (pos/neg) for each channel - value, sign, and lag
pnmat                           = [ pA nA ];
pnidx                           = [ pidx nidx ];
[ ~, midx ]                     = max( [ pA abs( nA ) ], [], 2 );
idx                             = ( midx - 1 ) * m + ( 1 : m )';
vA                              = pnmat( idx );
extidx                          = pnidx( idx );

% determine time lag of local exterma relative to global extremum
[ ~, maxidx ]                   = max( abs( vA ) );                         % [channels] get the origin
orig                            = extidx( maxidx );                         % [samples]
vT                              = extidx - orig;                         	% shift origin
vT                              = vT / spkFs * 1000;                        % convert to ms
nT                              = nidx - orig;
pT                              = pidx - orig;
nT                              = nT / spkFs * 1000;
pT                              = pT / spkFs * 1000;

% calculate peak-to-trough

%---------------------------------------------------------------
% channel specific BPI and indicator according to amplitude and ordinal considerations
%---------------------------------------------------------------
% vB: from -1 to 1
% 0     := symmetric (perfectly bipolar)
% 1     := perfect Pspike
% -1    := perfect Nspike

% note that to convert to a BPI one can 
% bpi = 1 - abs( vB )

% one positive and one negative local extremum per channel
kmat                            = [ pA nA ];
kmat1                           = kmat; % for general bpi (for every unit,no conditions)

% find the kmat values in the TH matrix
if size(TH,1)>1
    TH1 = TH';
    if ~THflag
    for i = 1: length(pA)
       idxthPA = find(ismember(wt(:,i),pA(i)));
       idxthNA = find(ismember(wt(:,i),nA(i)));
%        if mean_amp(i)>0
%           THmat (i,:) = [TH1(idxthPA(1),i)+mean_amp(i) TH1(idxthNA(1),i)+mean_amp(i)];
%        else
%           THmat (i,:) = [TH1(idxthPA(1),i)-mean_amp(i) TH1(idxthNA(1),i)-mean_amp(i)];
%        end
       THmat (i,:) = [TH1(idxthPA(1),i) TH1(idxthNA(1),i)];
       % find the tp if the positive peak comes before the negetive peak
       if i == maxidx && idxthNA>idxthPA
           tp = abs(idxthPA - idxthNA); % in samples
           tp = tp / spkFs * 1000;      % in ms
       end
    end
    
    elseif THflag
        for i = 1: length(pA)       
           idxthPA = find(ismember(wt(:,i),pA(i)));
           idxthNA = find(ismember(wt(:,i),nA(i)));
           tp_temp = abs(idxthPA(1)-idxthNA(1));
           idx_pre = idxthPA(1)-tp_temp;
           if idx_pre<1
               idx_pre=1;
           end
           THmat (i,:) = [TH1(idx_pre,i)+wt(idx_pre,i) 2*TH1(idxthNA(1),i)];
           % find the tp if the positive peak comes before the negetive peak
           if i == maxidx && idxthNA>idxthPA
               tp = abs(idxthPA - idxthNA); % in samples
               tp = tp / spkFs * 1000;      % in ms        
           end
        end
    end
    %kmat_temp = [w1(idxthPA) w1(idxthNA)];
    TH = THmat;
else
       idxthPA = find(ismember(wt(:,maxidx),pA(maxidx)));
       idxthNA = find(ismember(wt(:,maxidx),nA(maxidx))); 
       tp = abs(idxthPA - idxthNA); % in samples
       tp = tp / spkFs * 1000;      % in ms
end
    rA = pA(maxidx)-nA(maxidx);


% keep only cases with both extrema above TH
kmat( abs( kmat ) < TH )        = NaN;
if sum(sum(~isnan(kmat)))>0
    abvTH = true;
end
nans                            = sum( isnan( kmat ), 2 ) > 0;
kmat( nans, : )                 = NaN;

% keep only cases in which positive extremum precedes negative
kmat( pT > nT, : )              = NaN;

% derive BPI for each channel
% ( a - abs( b ) ) / ( a + abs( b ) )
a                               = kmat( :, 1 );
b                               = kmat( :, 2 );
vB                              = ( a - abs( b ) ) ./ ( a + abs( b ) );
%cbpi                            = 1 - abs( vB );

% determine which extremum is larger
pol                             = -2 * midx + 3; % 1: positive, -1: negative 

%---------------------------------------------------------------
% unit level metrics
%---------------------------------------------------------------
if size (TH,1)>1
    vidx                            = abs( vA ) >= TH(:,1);  
else
    vidx                            = abs( vA ) >= TH;
end
% determine whether all valid extrema are homogenous (Pspikes: 1; Nspikes: -1) 
% or heterogenous (some Pspikes and some Nspikes: 0)
upol                            = unique( pol( vidx ) );
if isempty( upol )
    upol                        = NaN;
elseif length( upol ) > 1
    upol                        = 0;
end

% determine the BPI (vB) according to the channel with the maximal p2p for bipolar units
bpi                             = vB( maxidx );

% determine the vB for all units
a1                               = kmat1( :, 1 );
b1                               = kmat1( :, 2 );
vB2                              = ( a1 - abs( b1 ) ) ./ ( a1 + abs( b1 ) );

% estimate the probability of multi-compartmental recording by the variance in the peak lag
tV                              = nanstd( vT( vidx ), 0, 1 );            	

%---------------------------------------------------------------
% graphics
%---------------------------------------------------------------

if ~graphics
    return
end

newplot

% raw waveform
xx                          = ( ( -n / 2 ) : 1 : ( n / 2 - 1 ) ) / spkFs * 1000; % [ms]
ylims                       = [ min( nA ) max( pA ) ];
xlims                       = xx( [ 1 n ] );
for i = 1 : m
    ch                      = m - i + 1;            % bottom channel is 1
    if m > 10
        ff = m/9;
        
        axes( 'position', [ 0.13 0.1/ff * ch 0.21341 0.1/ff ] )
    else
        axes( 'position', [ 0.13 0.1 * ch 0.21341 0.1 ] )
    end
    %subplot( m, 3, ( i - 1 ) * 3 + 1 )
    plot( xx, w( ch, : ), 'k', nT( ch ), nA( ch ), '.b', pT( ch ), pA( ch ), '.r' )
    set( gca, 'tickdir', 'out', 'box', 'off', 'xlim', xlims, 'ylim', ylims );
    alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
    axis off
    if i == m
        calibration( [ 0.5 100 ], [], [], { 'ms', '\muV' } );
    end

end

% same, as an upsampled color map:
USF                         = 50;
x0                          = ( ( -n / 2 ) : ( n / 2 - 1 ) ) / spkFs * 1000;
y0                          = 1 : m;
%[ x1, y1, w1 ]              = imupsample( x0, y0, w, USF );
[ csd0 ]                    = csd_mat (w);
[ x1, y1, csd1 ]            = imupsample( x0, y0, csd0, USF );
subplot( 1, 3, 2 )
imagesc( x1, y1, csd1 )
axis xy
colormap( myjet )
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 'CSD') )

% plot vA
subplot( 3, 3, 3 )
plot( vA, y0, '.-k' )
hold on
ph = plot( vA( ~vidx ), y0( ~vidx ), '.' );
set( ph, 'color', [ 1 1 1 ] * 0.8 )
ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
xlim( ylims + [ -1 1 ] * 10 )
alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Amplitude [{\mu}V]' )
title( sprintf( '%0.3g \\muV', vA( maxidx ) ) )

% plot vT
subplot( 3, 3, 6 )
plot( vT, y0, '.-k' )
hold on
ph = plot( vT( ~vidx ), y0( ~vidx ), '.' );
set( ph, 'color', [ 1 1 1 ] * 0.8 )
ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
xlim( xlims )
alines( 0, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
set( gca, 'tickdir', 'out', 'box', 'off' );
xlabel( 'Time lag [ms]' )
title( sprintf( '%0.3g ms, SD = %0.3g', vT( maxidx ), tV ) )

% plot vB
subplot( 3, 3, 9 )
plot( vB, y0( : ), '.-k' )
ylim( y0( [ 1 end ] ) + [ -1 1 ] * 0.5 )
xlim( [ -1 1 ] )
alines( bpi, 'x', 'linestyle', '--', 'color', [ 0 0 0 ] );
alines( maxidx, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );
set( gca, 'tickdir', 'out', 'box', 'off' );
title( sprintf( 'BPI=%0.3g', bpi ) )

return

% EOF

