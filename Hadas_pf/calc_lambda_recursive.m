% calc_lambda      calculates the baseline mean firing rate on track after outlier removal  
%
% CALL             [lambda, min_lambda, ratemap, occumap, orig_lambda] = calc_lambda( spk, pos, periods )
%
% GETS              
%                  spk            spike times, sampled at spkFs
%                  pos            position of the mouse, sampled at movFs
%                  periods        time of the begging and end of each trial, sampled at spkFs
%
% OPTIONAL
%                  spd            {[]} movement speed, passed to calc_field.m
%                  minspd         {[5]} minimum speed, to prune trials longer than it [cm/s], passed to calc_field.m
%                  binSize        {[5]} size of spatial bins [cm]
%                  spkFs          {[20000]} sampling frequency of spikes [Hz]
%                  movFs          {[39.0625]} sampling frequency of position and speed [Hz]
%                  xSmooth        {[1]} spatial smoothing, passed to calc_field.m
%                  LT_lims        {[80 220]} xlimits of the linear track 
%                  sd_mult        {[1]} multipications of sd above which values are prunned from the mean, inf = no prunning
%                  graphics       {[0]}
% STEPS
%                   (1) calculate the rate and occupancy maps
%                   (2) calculate mean firing rate 
%                   (3) find values sd*sd_mult above and below the weighted mean in the mean rate vector and change them to NaN
%                   (4) re-calculate weighted mean
%
% RETURNS           
%                   lambda           weighted mean firing rate during periods, with outliers removed
%                   min_lammbda      lower precentile of firing rate during periods
%                   ratemap          rate map of the unit
%                   occumap          occupancy map of the mouse
%                   original_lambda  weighted mean firing rate during periods, without removing outliers
%
% CALLS             ParseArgPairs, calc_field, calc_com, calc_lambda_outliers,  patch_band, alines, calc_sem

% written by        HS + ES 14-Jan-20
% modified          HS      21-Jan-20 added min_lambda as output and precentile as optional input
%                   HS      03-Jun-20 added prunning of extreme values from mean, added graphics
%                   HS + ES 07-jun-20 further modified 
%                   HS + ES 08-Jun-20 cleaned up
%                   HS + ES 18-jun-20 (1) renamed as recursive
%                                     (2) changed algorithmic details for the recursion
%                                     (3) changed occupancy vector to be sum over trials (not mean)
%                   HS      19-Jul-20 changed periods values to be in spkFs instead of movFs
%                   HS      26-Aug-20 (1) added smoothing of occupancy map
%                                     (2) changed mean calculation in graphics to weighthed mean (like analyze_LT_dir)
%                   HS      29-Sep-20 (1) added prunning of occurances with speed lower than minspd before running calc_field
%                                     (2) updated help

function [lambda_final, ratemap, occumap] = calc_lambda_recursive( spk, pos, periods, varargin )

% constants
ndirs                           = 2;

% arguments
[ spd, dirs, binSize, spkFs, movFs, minspd, xSmooth, LT_lims ...
    , sd_mult, alpha, doqbins, nrounds ...
    , graphics]        = ParseArgPairs (...
    { 'spd', 'dirs', 'binSize', 'spkFs', 'movFs', 'minspd', 'xSmooth', 'LT_lims' ...
    , 'sd_mult', 'alpha', 'doqbins', 'nrounds' ...
    , 'graphics'} ...
    ,{ [], [], 5, 20000, 39.0625, 10, 5, [80 220] ...
    , 1.5, 0.05, 1, 5 ...
    , 0 }, varargin{ : } );


%--------------------------------------------------------------------%
% check inputs
%--------------------------------------------------------------------%
nargs                     = nargin;

if nargs < 3 || isempty( spk ) || isempty( pos ) || isempty( periods )
    error( 'missing arguments' )
end

if ~isempty( dirs )
    if size( dirs, 1 ) ~= size( periods, 1 )
        error( 'input size mismatch' )
    end
    if setdiff( unique( dirs ), 1 : ndirs )
        error( 'unexpected dirs' )
    end
end

%--------------------------------------------------------------------%
% calculate rate and occupancy maps
%--------------------------------------------------------------------%

edges                              = ( LT_lims( 1 ) : binSize : LT_lims( 2 ) )';
binC                               = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;

% get periods and spikes in sec
s1                                 = spk / spkFs;
periods_sec                        = periods / spkFs;

% intiate
nxbins                             = length(binC);
ntrials                            = size(periods_sec,1);
ratemap                            = NaN(nxbins,ntrials);
occumap                            = NaN(nxbins,ntrials);

% prune occurances with low speed
pos2                               = pos(:,1);
pos2(spd < minspd )                = 1000;

% calculate place field per trial
for tr = 1 : ntrials
    
    [ratemap(:,tr), ~, ~, occumap(:,tr)] = calc_field( s1, pos2, edges, movFs, periods_sec( tr, : ), minspd, [], [], spd, xSmooth );
    
end

% smooth occupancy map
SDx                                             = xSmooth / binSize;
nx                                              = ceil( 6 * SDx ) + mod( ceil( 6 * SDx ) + 1, 2 );
w                                               = gausskernel( SDx, 0, nx, 1 );
w                                               = w / sum( w );
occumap                                         = firfilt( occumap, w );


%--------------------------------------------------------------------%
% calculate mean lambda per each directions
%--------------------------------------------------------------------%
nbins                           = size( ratemap, 1 );   % number of spatial bins
lambda                          = zeros( nrounds, ndirs );
rate_vec_dirs                   = NaN( nbins, ndirs );
occu_vec_dirs                   = NaN( nbins, ndirs );

for j = 1 : ndirs % loop over directions
   
    % select the relevant trials and compute firing rate vectors
    ratemap_dir                 = ratemap( :, dirs == j );
    occumap_dir                 = occumap( :, dirs == j );
    
    rate_vec                    = calc_com( ratemap_dir', occumap_dir' ); % weighted mean firing rate
    rate_vec                    = rate_vec( : );
    occu_vec                    = nansum( occumap_dir, 2 );
    
    % initialize the vectors to be pruned
    rate_vec_running            = rate_vec;
    occu_vec_running            = occu_vec;

    % run the recursion
    for i = 1 : nrounds
        % compute the lambda without outliers
        lambda(i,j)                 = calc_lambda_outliers( rate_vec_running, occu_vec_running, sd_mult );
        % detect significant bins (or continuous stretches) and define fields
        limFields                   = field_poisson_test( rate_vec, occu_vec, lambda(i,j) ...
            , 'alpha', alpha, 'FRthr_pq', [ NaN NaN ], 'doqbins', doqbins );
        % apply pruning within fields
        if ~isempty(limFields)
            idx                     = enumerate( limFields ); 
            rate_vec_running( idx ) = NaN;
            occu_vec_running( idx ) = NaN;
        end
    end
  
    rate_vec_dirs( :, j )           = rate_vec_running;
    occu_vec_dirs( :, j )           = occu_vec_running;
    
end

   
lambda_final                        = calc_lambda_outliers( rate_vec_dirs( : ), occu_vec_dirs( : ), sd_mult );    


%--------------------------------------------------------------------%
% graphics
%--------------------------------------------------------------------%
if graphics
    figure
    ylims                           = zeros( ndirs, 2 );
    for j                           = 1 : ndirs
        
        % compute
        [ rate_vec, sd_rate ]                      = calc_com( ratemap( :, dirs == j )', occumap( :, dirs == j )' ); % weighted mean firing rate
        rate_map_sem                               = sd_rate / sqrt( size( ratemap( :, dirs == j ), 2 ) );

        % plot
        subplot( 1, ndirs, j );
        ph                          = patch_band( binC, rate_vec, rate_map_sem );
        hold on
        ph                          = plot( binC, rate_vec, '.' ); 
        set( ph, 'MarkerSize', 20, 'color', [ 0 0 1 ] );
        a1 = alines( lambda_final, 'y','color', 'magenta','LineStyle','-','LineWidth',2 );
        a2 = alines( lambda( 1, j ), 'y','color', 'red','LineStyle','--','LineWidth',2 );
        a3 = alines( lambda( 2, j ), 'y','color', 'green','LineStyle','--','LineWidth',2 );
        legend( [a1 a2 a3], {'\lambda', '\lambda1', '\lambda2'} );
        xlabel('Position [cm]');
        ylabel('Firing rate [spks/s]');
        title(sprintf('Direction %d',j));
        ylims( j, : ) = ylim;
    end

    ylims = [ 0 max( ylims( : ) ) ];
    
    for j                           = 1 : ndirs
        % compute
        rate_vec                    = calc_com( ratemap( :, dirs == j )', occumap( :, dirs == j )' ); % weighted mean firing rate
        rate_vec                    = rate_vec(:);
        limFields                   = field_poisson_test( rate_vec, occu_vec, lambda_final ...
            , 'alpha', alpha, 'FRthr_pq', lambda_final * [ 1 1 ], 'doqbins', 0 );
        
        % expand limFields by half bin towards edges
        if ~isempty(limFields)
            mat0                        = reshape( binC( limFields ), size( limFields ) );
            mat1                        = repmat( diff( binC( 1 : 2 ) ) / 2 * [ -1 1 ], [ size( limFields, 1 ) 1 ] );
            limFieldsPos                = mat0 + mat1;
        else
            limFieldsPos                = [];
        end
        
        % plot
        subplot( 1, ndirs, j )
        set( gca, 'tickdir', 'out', 'box', 'off', 'ylim', ylims );
        alines( limFieldsPos, 'x', 'color', [ 0 0 0 ], 'linestyle', '--' );
    end    
    
end

return;

%EOF