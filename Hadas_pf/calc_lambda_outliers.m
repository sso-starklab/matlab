function lambda = calc_lambda_outliers( mean_rate, mean_occupancy, sd_mult)

%--------------------------------------------------------------------%
% calculate mean lambda
%--------------------------------------------------------------------%
px                               = mean_occupancy / nansum( mean_occupancy ) ;        % probability to be in a given spatial bin
orig_lambda                      = nansum( mean_rate .* px );                        % weighted mean firing rate

%--------------------------------------------------------------------%
% prune extreme values and recalculate
%--------------------------------------------------------------------%
if sd_mult == inf
    lambda                       = orig_lambda;
    return
end
[ myu, sd ]                  = calc_com( mean_rate, px);
hi                           = mean_rate > ( myu + sd * sd_mult );
lo                           = mean_rate < ( myu - sd * sd_mult );
ratevec_pr                   = mean_rate;
ratevec_pr( hi | lo )        = NaN;
px( hi | lo )                = NaN;
px                           = px / nansum( px );
lambda                       = nansum( ratevec_pr .* px );

return;

%EOF
