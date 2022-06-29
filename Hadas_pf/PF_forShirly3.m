
% parameters
binSize = 2.5; %[cm]
minspd  = 10; %[cm/s]
xSmooth = 5; %[cm]

% load position
mov_name                        = sprintf('%s.mov.mat', filebase);
load( mov_name, '-mat', 'mov');
pos = mov.pos; %[cm]
spd = spd.pos;

% calculate trials
[ trials, tmov ]                = get_LT_trials( filebase );
periods                         = trials(:,1:2);
dir                             = trials(:,3);

% find track limits
LL1                             = max( pos(periods( rdir == 2 ,1) ) );
LL2                             = max( pos(periods( rdir == 1 ,2) ) );
HL1                             = min( pos(periods( rdir == 2 ,2) ) );
HL2                             = min( pos(periods( rdir == 1 ,1) ) );
LL                              = max([LL1,LL2]);
HL                              = min([HL1,HL2]);
pos                             = pos - LL;
LT_lims                         = [0, HL - LL];
pos_dirs                        = [ LT_lims( 2 ) + LT_lims( 1 ) - pos pos ];


s                                 = load_spikes( filebase, [], 'B' );
units                             = 1 : size( s.shankclu, 1 );
nunits                            = length(units);

for i = 1 : nunits

    % get spk of a specific unit
    ui                         = units(i);
    ashankclu                 = s.map( ui, 2 : 3 );
    aclunum                   = s.map( ui, 1 );
    atype                     = s.shankclu( ui, 3 );
    idx                       = s.clu == aclunum;
    spk                       = s.res( idx ); 

    % calculate baseline firing rate
    lambda                        = calc_lambda_recursive( spk, pos, periods,'dirs',dir,'spd',spd,'binSize',binSize);
    % spk - spike time (from s.res, for the relevant unit)
    % pos - position
    % periods - times of trials in spkFs

    % generate rate and occupancy maps
    s1                                = spk / spkFs;
    periods_sec                       = periods / spkFs;
        edgeMin                                         = ceil( LT_lims( 1 ) / binSize ) * binSize;
        edgeMax                                         = floor( ( LT_lims( 2 ) ) / binSize ) * binSize;
        edges                                           =  edgeMin : binSize : edgeMax;
        binC                                            = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
        nxbins                                          = length(binC);
        
   for d = 1 :2 
    
        my_periods                                      = periods_sec(dir == d, :);
        ntrials                                         = size(my_periods,1);
        trials_rate                                     = NaN(nxbins,ntrials);
        trials_occu                                     = NaN(nxbins,ntrials);
        trials_count                                    = NaN(nxbins,ntrials);

        % prune occurances with low speed
        pos2                                            = pos_dirs(:,d);
        pos2(spd <minspd )                              = 1000;

        % calculate place field per trial
        for tr = 1 : ntrials

            [trials_rate(:,tr), ~,trials_count(:,tr), trials_occu(:,tr)] = calc_field( s1, pos2, edges, movFs, my_periods( tr, : ), minspd, [], [], spd, xSmooth );

        end

        % smooth occupancy map
        SDx                                             = xSmooth / binSize;
        nx                                              = ceil( 6 * SDx ) + mod( ceil( 6 * SDx ) + 1, 2 );
        w                                               = gausskernel( SDx, 0, nx, 1 );
        w                                               = w / sum( w );
        trials_occu                                     = firfilt( trials_occu, w );

        % calculate mean firing rate
        [ m_rate, sd_rate ]                             = calc_com( trials_rate', trials_occu' ); % weighted mean firing rate
        sum_occu                                        = nansum(trials_occu,2);
        [limFields, fields_stat]                        = field_poisson_test(m_rate, sum_occu,lambda, 'FRthr_pq',[NaN NaN]);

        % limFields - limits of place fields in cm
        % field_stat- [mean FR, SD of FR, mean pval, area, FR outside the field, peakFR, field_size] 
   end
   
end