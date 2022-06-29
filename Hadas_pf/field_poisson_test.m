% field_poisson_test   Finds significant fields within a vector and outputs their limits and p value
%
% CALL                 [ limFields, fields_stat ] = field_poisson_test(rate_map, occupancy_map, lambda)
%
%
% GETS                 rate_map            rate in independent variable bins [spikes/s]
%                      occupancy_map       the time the mouse spent in each bin [s]
%                      lambda              weighted mean firing rate across all bins (if not entered, calculated within the function)
%
% OPTIONAL ARGUMENTS
%                      alpha               {0.05}
%                      doqbins             {[0]} calculate negative fields in addition to positive fields
%                      max_gain            {inf}    maximum possible gain (larger brain is clipped to this value 
%                      FRthr_pq            {[lambda lambda]} threshold for field limits for positive and negative fields
%                      graphics            {0}
% STEPS
%                      (1) calculate the cumulative Poisson propability to get
%                          the rate in each bin in comparison to lambda, by calling
%                          poissonTest
%                      (2) perform benferroni correction for alpha
%                          according to effective indepndent measurments
%                      (3) determine whether there is a place field according to log likelhood 
%                      (4) find field limits
%                      (5) calculate field statistics
%
% RETURNS
%                     pbins               p value per bin ( length = nbins )
%                     limFields           a matrix of the bins in which each significant field segment begins and ends (#F x 2 matrix)
%                     field_stat          a matrix of #F x 7, with the statistics of each field: 
%                                         [mean FR, SD of FR, mean pval, area, FR outside the field, peakFR, field_size] 
%
% CALLS             poissonTest, effective_independence, parse, isoverlap, nangeomean, patch_band
%

% written by        HS + ES 12-12-19
% modified          HS      17-12-19 fixed lambda size error
%                   HS      18-12-19 fixed field_stat and added a fourh statistic (gain)
%                   HS      29-12-19 changed gain to be calculated  according to lambda_hat (rate outside field) instead of lambda
%                   HS      30-12-19 changed alpha correction from alpha/bins to alpha/indepndent bins using  autocorrelation (see effective_indepndence)
%                   HS      05-02-20 added peakFR to outputs of field_stat, and added clipping of gain 
%                   HS      18-02-20 added field size to outputs of field_stat
%                   HS      17-06-20 (1) changed lambda to an obligatory input instead of optional
%                                    (1) changed the method of finding fields to log likelihood
%                                    (2) changed the methods of finding field limits to above FRthr
%                                    (3) added graphics
%                   HS + ES 18-06-20 (1) changed FR to be called SC (this is the spike count per bin, after smoothing)
%                                    (2) added qbins
%                   HS + ES 19-06-20 cleaned up
%                   HS      12-08-20 fixed an error in the calculation of firing rate outside field (and therfore gain calculations)
%                   HS      19-08-20 changed default of max_gain from 100 to inf 
%                   HS      11-11-20 changed the output of fields_stat(:,5) from field gain = mean_field / mean_outside_field to mean_outside_field. Can calculate gain post hoc by fields_stat(1,6) / fields_stat(:,5) 

function [ limFields, fields_stat ] = field_poisson_test(rate_map, occupancy_map, lambda, varargin)

%--------------------------------------------------------------------%
% check input
%--------------------------------------------------------------------%

nargs = nargin;

if nargs < 3 || isempty( rate_map ) || isempty( occupancy_map ) || isempty( lambda )
    error( 'missing arguments' )
end


[ alpha, doqbins, FRthr_pq, graphics ] = ParseArgPairs (...
    { 'alpha', 'doqbins', 'FRthr_pq', 'graphics' } ...
    ,{ 0.05, 0, [ lambda lambda ], 0 ...
    },varargin {:} );

if length( rate_map ) ~= length( occupancy_map )
    error('rate map and occupancy map are not in the same length');
end

if numel( lambda ) > 1
    error('lambda should be a scalar');
end

% initialize outputs
limFields                           = [];
fields_stat                         = [];

%--------------------------------------------------------------------------------------%
% test whether there is a place field
%--------------------------------------------------------------------------------------%

% check Poisson propability per bin
nbins                               = size( rate_map, 1 );
lambda_FP                           = repmat( lambda, nbins, 1 );
SC                                  = rate_map .* occupancy_map;
[ pbins, qbins ]                    = poissonTest( lambda_FP, SC, occupancy_map );

% correct alpha according to bonferoni (effective indepndent measurments)
m                                   = effective_independence ( rate_map ); %get the number of effective indepndent measurment for bonferoni correction
cutoff                              = alpha / ( ( double( doqbins ) + 1 ) * m ); % if doqbins, 2 in case of doqbins since both pbins and qbins

% check significance 
sig_p                               = pbins <= cutoff;
if doqbins
    sig_q                           = qbins <= cutoff;
else
    sig_q                           = false( nbins, 1 );
end
sig_bins                            = sig_p | sig_q;
is_field                            = any(sig_bins);


%--------------------------------------------------------------------------------------%
% find and expand field limits
%--------------------------------------------------------------------------------------%

sig_limits                          = parse( find( sig_bins ) );


if is_field && all( ~isnan( FRthr_pq ) )

    thr_limits_p                    = parse( find( rate_map > FRthr_pq( 1 ) ) );
    bp                              = isoverlap( thr_limits_p, sig_limits );
    limFieldsP                      = thr_limits_p( bp, : );  

    if doqbins
        thr_limits_q                = parse( find( rate_map < FRthr_pq( 2 ) ) );
        bq                          = isoverlap( thr_limits_q, sig_limits );
        limFieldsQ                  = thr_limits_q( bq, : );  
    else
        limFieldsQ                  = [];
    end
    
    limFields                       = uniteranges( limFieldsP, limFieldsQ );
    
else
    
    limFields                       = sig_limits;
    limFieldsP                      = limFields;
    
end


%--------------------------------------------------------------------------------------%
% calculate lambda_hat, the mean firing rate outside of all place fields
%--------------------------------------------------------------------------------------%

if is_field 
    
    idx                                    = limFieldsP(1):limFieldsP(2);
    not_idx                                = true(nbins,1);
    not_idx(idx)                           = false;
    outlim_rate                            = rate_map ( not_idx );
    outlim_occu                            = occupancy_map ( not_idx );
    px                                     = outlim_occu / nansum( outlim_occu ) ;     % probability to be in a given spatial bin
    lambda_hat                             = nansum( outlim_rate .* px );              % weighted mean firing rate


%--------------------------------------------------------------------------------------%
% calculate mean rate and mean p value within each field
%--------------------------------------------------------------------------------------%

    nFields                             = size(limFieldsP,1);
    fields_stat                         = zeros(nFields,3);

    for f = 1 : nFields

        idx                             = limFieldsP( f, 1 ) : limFieldsP( f, 2 );
        fields_stat(f,1)                = mean ( rate_map( idx  ) );
        fields_stat(f,2)                = std( rate_map( idx  ) );
        fields_stat(f,3)                = nangeomean ( pbins( idx ) );
        fields_stat(f,4)                = sum( rate_map( idx  ) .* occupancy_map( idx ) ) / sum( occupancy_map( idx ) );
        fields_stat(f,5)                = lambda_hat;
        fields_stat(f,6)                = max (rate_map (idx) );
        fields_stat(f,7)                = limFieldsP( f, 2 ) - limFieldsP( f, 1 ) + 1;
    end
    
end


%--------------------------------------------------------------------------------------%
% graphics
%--------------------------------------------------------------------------------------%


if graphics
    
   figure;
    plot(rate_map);
    ylabel('mean firing rate [spk/s]')
    xlabel('spatial bins');
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines(limFieldsP,'x', 'color', 'red','LineStyle','--','LineWidth',1);
    if doqbins
       alines(limFieldsQ,'x', 'color', 'blue','LineStyle','--','LineWidth',1);
    end
    alines(lambda,'y', 'color', 'black','LineStyle','--','LineWidth',1);
    t                                   = sprintf ('Lambda = %.2f', lambda);
    title(t);
    
    
end

return

%EOF

