% csd_mat       calculates csd per column (n) for a matrix of m x n 
%
% call          [ csd ] = csd_mat ( mat )
% 
% gets          mat         matrix of m x n (for example channels x samples)
%
% returns       csd         matrix of m x n padded with NaNs with the respective CSDs

% SSo 13-sep-21

function [ csd ] = csd_mat (mat)

s = size (mat);
mat = [ NaN(1,s(2)); mat; NaN(1,s(2))];         % pad with NaNs
csd = [];                                       % intialize 
for c = 1 : s(2)                                % columns
    for i = 2 : (size(mat,1)-1)                 % rows
        j = i-1;
        csd(j,c) = 2 * mat (i,c) - (mat (i+1,c) + mat (i-1,c));
    end
end

end


% mat = sst1.mean{1,1};

% per one column (sample)
% vec = sst1.mean{1,1}(:,16);
% 
% vec = [1; 2; 3; 4; 3; 2; 1];
% vec = [4; 3; 2; 1; 2; 3; 4];
% 
% vec = [ NaN; vec; NaN];
% csd = [];
% for i = 2 : (length(vec)-1)
%    j=i-1;
%    csd(j) = 2 * vec (i) - vec (i+1) - vec (i-1);
% end
% csd = csd';
