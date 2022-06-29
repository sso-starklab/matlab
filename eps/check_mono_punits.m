
load ('D:\_Shirly\Lab\matlab\eran_code\res');
punits_fb = res.filebase;
punits_clu = res.shankclu (:,1:2);
npunits = length(punits_clu);
idx = [];
n_synp = [];

for i=1:npunits
mono = check_mono (punits_fb(i));

r = punits_clu(i) == mono.shankclu;
idx2 = mono.pairsExc (:,1) == r;
n_r = sum (idx2);
if n_r>1
    idx=[idx,i];
    n_synp=[n_synp,n_r];
end
end

punits_pre = length(idx);
for j=1:punits_pre 
    unit = idx(j);
plot_ss (punits_fb(unit), punits_clu(unit));
end