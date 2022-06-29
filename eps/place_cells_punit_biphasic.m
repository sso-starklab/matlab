% Place cells examples:
filename = 'mC41_21';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [2 28]); % nunit
figure, plot_ss (filebase, [1 15]); % punit - used

% ind = find(contains (sst.filebase, 'mC41_21') & isbip ==1)
% for i=1:length(ind)
%     shanclui = sst.shankclu(ind(i),:);
%     figure, plot_ss (filebase, shanclui); % biphasic
% end
figure, plot_ss (filebase, [1 14]); % biphasic

% another possible session for examples
filename = 'mC41_44';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [4 16]); % nunit
figure, plot_ss (filebase, [2 12]); % biphasic

% ind = find(contains (sst.filebase, 'mC41_21') & ispos ==1)
for i=1:length(ind)
    shanclui = sst.shankclu(ind(i),:);
    figure, plot_ss (filebase, shanclui); % punit
end

filename = 'mC41_40';filebase = filebase_lookup(filename,1);
figure, plot_ss (filebase, [2 12]); % biphasic - used

fig1( 1:12 ) = make_PUNIT_figures( 1:12, 'TPR_23-Aug-2021_Z6_extremum_bpi','/media/shirly/Data/_Shirly/Lab/eps/data/', 'sst', sst );