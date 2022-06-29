%pie chart for seesion with punits (Y) vs sessions without punits(N)
counts = [89 120-89];
code_names = {'Y', 'N'};
figure;
code_colors         = [ 1 0 0; 0 0.7 0; 0 0 0.7; 0.7 0.7 0.7 ];
ph                  = pie( counts, code_names );
for i= 1:2
set( ph( 2 * i - 1 ), 'FaceColor', code_colors( i, : ), 'EdgeColor', code_colors( i, : ) )
set( ph( 2 * i ), 'String', sprintf( '%s (%d)', code_names{ i }( 1 ), counts( i ) ) )
end
 

%pie chart for the fraction of the punits out of the total units
counts = [286 3541];
code_names = {'P', 'N'};
figure;
code_colors         = [ 1 0 0; 0 0.7 0; 0 0 0.7; 0.7 0.7 0.7 ];
ph                  = pie( counts, code_names );
for i= 1:2
set( ph( 2 * i - 1 ), 'FaceColor', code_colors( i, : ), 'EdgeColor', code_colors( i, : ) )
set( ph( 2 * i ), 'String', sprintf( '%s (%d)', code_names{ i }( 1 ), counts( i ) ) )
end

%show number of punits and nunits per session, sort by number of 
%punits (small to large), add dashed line at mean number of punits

%go over each session in SST structure and find the number of nunits and
%punits

num_units = length(sst.filebase); %number of total units
num_session = length(unique(sst.filebase)); %number of session

session_name = unique(sst.filebase); %sessions name
num_per_sess = zeros(num_session,1); %total units per session
for j = 1:num_session
    for  i=1:num_units
    if ismember(session_name(j),sst.filebase(i))
        num_per_sess(j) =  num_per_sess(j)+1;
    end
    end
end

%now the same for res
pun_per_sess = zeros(num_session,1); %total units per session
pun_units = length(res.filebase);
for j = 1:num_session
    for  i=1:pun_units
    if ismember(session_name(j),res.filebase(i))
        pun_per_sess(j) =  pun_per_sess(j)+1;
    end
    end
end

%now sort by punits
[~,sidx] = sort (pun_per_sess);

pun_per_sess_sort = pun_per_sess(sidx,:);
num_per_sess_sort = num_per_sess(sidx,:);
%plot in bar
figure;
bar(num_per_sess_sort)
hold on
bar(pun_per_sess_sort,'r')
mean_units = mean(num_per_sess_sort);
mean_punits = mean(pun_per_sess_sort);
alines(mean_punits, 'y')
alines(mean_units, 'y')
%plot in fraction

%before sorting
frac_pun =pun_per_sess./ num_per_sess;
[~,sidx2] = sort (frac_pun);
frac_pun_sort = frac_pun(sidx2);
figure;
bar(frac_pun_sort)
mean_frac = mean(frac_pun_sort);
alines(mean_frac, 'y')





%number of punits as function of animal type per session:
%uniqe session with punits:
num_session = length(unique(sst.filebase)); %number of session

sessions_names = unique(res.filebase);
%all animal:
animalname = {'mV99', 'mF84', 'mC41','mK01','mF79','mF93','mP23','mP101','mF105','mO251','mDL5','mF108','mB142'};
%types of all animals
codes_names = {'VIP', 'CaMK', 'FVB-CaMK','CCK','CaMK','CaMK','PV','PV','FVB::c57B','PV','PV','CaMK','FVBxPV'};
sess_whatever = zeros(length(res.filebase),length(animalname));

session_name = unique(sst.filebase); %sessions name
pun_per_sess = zeros(num_session,1); %total units per session
pun_units = length(res.filebase);
for j = 1:num_session
    for  i=1:pun_units
    if ismember(session_name(j),res.filebase(i))
        pun_per_sess(j) =  pun_per_sess(j)+1;
    end
    end
end

for l = 1:length(animalname)
    for f= 1:length( pun_per_sess)
           if ismember(animalname{l},session_name{f})
               
           if l==1
               vip_units(f) =pun_per_sess(f);
           end
                      if ismember(l,[2 5 6 12])
               camK_units(f) =pun_per_sess(f);
                      end
                      if ismember(l,[3])
               fvb_camK_units(f) =pun_per_sess(f);
                      end
                                if ismember(l,[4])
              cck_units(f) =pun_per_sess(f);
                                end
                                                   if ismember(l,[7 8 10 11])
           pv_units(f) =pun_per_sess(f);
                                                   end
                                                   if ismember(l,[9])
              FVB_c57B_units(f) =pun_per_sess(f);
                                                   end
                                                   if ismember(l,[13])
              FVB_pv_units(f) =pun_per_sess(f);
                                                   end
           end
    end
end

%SEM = std(x)/sqrt(length(x));
vip_units_mean = mean(nonzeros(vip_units));
vip_units_SEM = std(nonzeros(vip_units))/sqrt(length(nonzeros(vip_units)));

camK_units_mean = mean(nonzeros(camK_units));
camK_units_SEM = std(nonzeros(camK_units))/sqrt(length(nonzeros(camK_units)));           

fvb_camK_units_mean = mean(nonzeros(fvb_camK_units));
fvb_camK_units_SEM = std(nonzeros(fvb_camK_units))/sqrt(length(nonzeros(fvb_camK_units)));

cck_units_mean = mean(nonzeros(cck_units));
cck_units_SEM = std(nonzeros(cck_units))/sqrt(length(nonzeros(cck_units)));

pv_units_mean = mean(nonzeros(pv_units));
pv_units_SEM = std(nonzeros(pv_units))/sqrt(length(nonzeros(pv_units)));

% FVB_c57B_units_mean = mean(nonzeros(FVB_c57B_units));
% FVB_c57B_units_SEM = std(nonzeros(FVB_c57B_units))/sqrt(length(nonzeros(FVB_c57B_units)));
FVB_c57B_units_mean=[];
FVB_c57B_units_SEM = [];

FVB_pv_units_mean = mean(nonzeros(FVB_pv_units));
FVB_pv_units_SEM = std(nonzeros(FVB_pv_units))/sqrt(length(nonzeros(FVB_pv_units)));


names = {'VIP', 'CaMK', 'FVB-CaMK','CCK','PV','FVB::c57B','FVBxPV'};
counts = [vip_units_mean camK_units_mean fvb_camK_units_mean cck_units_mean pv_units_mean FVB_c57B_units_mean FVB_pv_units_mean];
SEM =  [vip_units_SEM camK_units_SEM fvb_camK_units_SEM cck_units_SEM pv_units_SEM FVB_c57B_units_SEM FVB_pv_units_SEM];
errhigh = counts+(SEM/2);
 errlow  = counts-(SEM/2);
figure;
bar(counts)  
hold on
x= 1:length(counts);
er = errorbar(x,counts,errlow,errhigh); 
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

