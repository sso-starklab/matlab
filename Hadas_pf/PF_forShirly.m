%filename
%clear all

        %stuff
        tot_fields = [];
        tot_pcells = [];
        tot_PF_TE = [];
        tot_not_PF_TE =[];
        tot_PC_TE = [];
        tot_PC_not_TE = [];
        num_PC_TE= [];
        num_PC_not_TE= [];
        
        
        mat_f= [];
        mat_p = [];
        
sess_names_be = {'mS234_02','mS234_03','mS234_04'};
sess_names_af = {'mS234_16','mS234_17','mS234_18','mS234_19','mS234_20'};

    for sess = 1 : length (sess_names_be)
        filename                        = sess_names_be{sess};
        filebase = filebase_lookup(filename);
        % parameters
        binSize = 2.5; %[cm]
        minspd  = 10; %[cm/s]
        xSmooth = 5; %[cm]
        


        % load position
        mov_name                        = sprintf('%s.mov.mat', filebase);
        load( mov_name, '-mat', 'mov');
        pos = mov.pos(:,1); %[cm]
        spd = mov.spd;
        movFs = mov.Fs;
        %load spike times
        spk = LoadRes([ filebase '.res.1' ]);

        par                 = LoadXml( filebase );
        spkFs               = par.SampleRate;
        % calculate trials
        trials = [];
        while size(trials,1) < 20
        [ trials, tmov ]                = get_LT_trials( filebase, 'graphics',1, 'overwrite',1,'alphaLevel',0.1 ,'extractMode','byPOS');
        end
        periods                         = trials(:,1:2);
        dir                             = trials(:,3);
        
%         S_RUFF=0;
%         if S_RUFF 
%             pos(pos<200) = 0;
%         end
        % find track limits
        LL1                             = max( pos(periods( dir == 2 ,1) ) );
        LL2                             = max( pos(periods( dir == 1 ,2) ) );
        HL1                             = min( pos(periods( dir == 2 ,2) ) );
        HL2                             = min( pos(periods( dir == 1 ,1) ) );
        LL                              = max([LL1,LL2]);
        HL                              = min([HL1,HL2]);
        LT_lims                         = [LL, HL];
        % pos_dirs                        = [ LT_lims( 2 ) + LT_lims( 1 ) - pos pos ];

        s                                 = load_spikes( filebase, [], 'B' );
        units                             = 1 : size( s.shankclu, 1 );
        nunits                            = length(units);
        periods1=trials(:,4:5);
        figtemp = figure;
        field_count = [0 0];
        
        % find the number of PF in the last 20 cm of the track (DW)
        RB = 20/binSize; %relevent_bins
        jj=0;
        jjj=0;
        PF_TE=[];
        PF_not_TE = [];
        field = [];
        for i = 1 : nunits

            % get spk of a specific unit
            ui                         = units(i);
            ashankclu                 = s.map( ui, 2 : 3 );
            aclunum                   = s.map( ui, 1 );
            atype                     = s.shankclu( ui, 3 );
            idx                       = s.clu == aclunum;
            spk                       = s.res( idx ); 

            % calculate baseline firing rate
            lambda                        = calc_lambda_recursive( spk, pos, periods1,'dirs',dir,'spd',spd,'binSize',binSize,'LT_lims', LT_lims );
            % spk - spike time (from s.res, for the relevant unit)
            % pos - position
            % periods - times of trials in spkFs

            % generate rate and occupancy maps
            s1                                = spk / spkFs;
            periods_sec                       = periods1 / spkFs;
            ntrials                           = size(periods1,1);
            edgeMin                                         = ceil( LT_lims( 1 ) / binSize ) * binSize;
            edgeMax                                         = floor( ( LT_lims( 2 ) ) / binSize ) * binSize;
            edges                                           =  edgeMin : binSize : edgeMax;
            binC                                            = ( edges( 1 : end - 1 ) + edges( 2 : end ) ) / 2;
            nxbins                                          = length(binC);


              for d = 1:2  
              subplot(1,2,d)
                my_periods                                      = periods_sec(dir == d, :);
                ntrials                                         = size(my_periods,1); 
                trials_rate                                     = NaN(nxbins,ntrials);
                trials_occu                                     = NaN(nxbins,ntrials);
                trials_count                                    = NaN(nxbins,ntrials);

                % prune occurances with low speed
                pos2                                            = pos;
                pos2(spd <minspd )                              = 1000;
                
                %simcki hack
%                 idx0 = pos2 > 200 & pos2 <1000;
%                 posdw = pos2;
%                 posdw(pos2 > 200 ) = 1000; 
%                 pos_dw = pos2(idx0);
%                 spd_dw = spd(idx0);
%                 my_periods_dw = my_periods(idx,:);
                
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
                [limFields, fields_stat]                        = field_poisson_test(m_rate', sum_occu,lambda, 'FRthr_pq',[NaN NaN]);
                if ~isempty(limFields)
                    field{i,d}.lim = limFields;
                    for l = 1: size(limFields,1)
                    plot(limFields(l,1):limFields(l,2),fields_stat(l,1)*ones(fields_stat(l,7),1))
                     hold on;
                    text(limFields(l,2),fields_stat(l,1),sprintf('%d %d',ashankclu))
                    if limFields(l,1) >= (nxbins - RB)
                       jj= jj+1;
                       PF_TE (jj,1:2) = ashankclu;
                    else
                        jjj = jjj+1;
                       PF_not_TE (jjj,1:2) = ashankclu;
                    end
        %             xlim([0 51])
                    end
                            field_count (d) = field_count (d)+l;

                else
                    field{i,d}.lim = NaN;
                end
                   field{i,3} = ashankclu; 
                % limFields - limits of place fields in cm
                % field_stat- [mean FR, SD of FR, mean pval, area, FR outside the field, peakFR, field_size]

              end
        end
        
        
        
        n1=0;
        n2=0;
        n3=0;
        field_num = length(field);
        for i = 1:field_num
              if ~isnan (field{i,1}.lim) & ~isnan (field{i,2}.lim)
                  field{i,4} = 3;
                  n3=n3+1;
              elseif ~isnan (field{i,1}.lim)
                  field{i,4} = 1;
                  n1=n1+1;
              elseif ~isnan (field{i,2}.lim)
                  field{i,4} = 2;
                  n2=n2+1;
              end
        end
        
        
        %to do
        % find out if the number of place fieald in the are of inretest is
        % larger when considering that the area of interest is small
        %done?

        tot_fields = sum(field_count);
        tot_pcells(sess) = sum(n1+n2+n3);
        tot_PF_TE(sess) = size(PF_TE,1);
        tot_not_PF_TE(sess) = size(PF_not_TE,1);
        
        uniqe_PC = unique(PF_TE,'rows');
        tot_PC_TE = [tot_PC_TE;uniqe_PC];
        
        num_PC_TE(sess) = size(uniqe_PC,1);
       
       
        
   
        uniqe_PC2 = unique(PF_not_TE,'rows');
        tot_PC_not_TE = [tot_PC_not_TE;uniqe_PC2];       
        num_PC_not_TE(sess) = size(uniqe_PC2,1);
        mat_f{sess} = [field_count(1) tot_fields; field_count(2) tot_fields];
        mat_p{sess} = [n1+n3 tot_pcells; n2+n3 tot_pcells];

        pv(sess) = gtest(mat_f{sess});
    end
    
    % DO NOT RUN
    
    
    PC1 = [tot_PC_not_TE tot_PC_TE];
     % for place field
    vec_after = [1/1 6/17 2/14 6/17 4/25];
    vec_before = [5/36 2/9 3/15];
    
    %for place cells
    vec_PC_before = [num_PC_TE./(num_PC_not_TE+num_PC_TE)];
    vec_PC_before =[0.2083    0.2857    0.2000]
    vec_PC_after = [ 1.0000    0.4545    0.2857    0.3636    0.2667]
    
    
    %plot
    figbox = figure;
    subplot(1,2,1)
dw = vec_before';
pw = vec_after';
g = [1 1 1 2 2 2 2 2];
 boxplot([dw;pw],g)
  
    subplot(1,2,2)
dw = vec_PC_before';
pw = vec_PC_after';
g = [1 1 1 2 2 2 2 2];
 boxplot([dw;pw],g)

chance_lvl = 8/51;
%PC
pv = myBinomTest(16,45,chance_lvl,'one')
pv = myBinomTest(9,41,chance_lvl,'one')

%PF
pv = myBinomTest(19,74,chance_lvl,'one')
pv = myBinomTest(10,60,chance_lvl,'one')


pv1 = utest(vec_PC_after, chance_lvl)
pv1 = utest(vec_PC_before, chance_lvl)

pv1 = utest(vec_PC_after, chance_lvl)
pv1 = utest(vec_PC_before, chance_lvl)