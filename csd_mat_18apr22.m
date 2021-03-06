% csd_mat       calculates csd per column (n) for a matrix of m x n 
%
% call          [ csd ] = csd_mat ( mat )
% 
% gets          mat         matrix of m x n (for example channels x samples)
%
% returns       csd         matrix of m x n padded with NaNs with the respective CSDs

% SSo 13-sep-21
% 
% revisions:
% 18-apr-22 (1) allocated space for csd
%           (2) ended function with 'return'
%           (3) added option for odd/even sited csd - issue

function [ csd1, csd2 ] = csd_mat_18apr22 (mat, oddEven)
if isempty(oddEven)
    oddEven = 0;
    csd2 = [];
end
if oddEven
    s = size (mat);
    s = ceil([s(1)/2 s(2)]);
    mat1 = [ NaN(1,s(2)); mat(2:2:end,:); NaN(1,s(2))];     % pad odds mat with NaNs
    mat2 = [ NaN(1,s(2)); mat(1:2:end,:); NaN(1,s(2))];     % pad even mat with NaNs
    csd1 = zeros( size( mat1, 1 ), size( mat1, 2 ) );  % intialize 
    csd2 = zeros( size( mat2, 1 ), size( mat2, 2 ) );  % intialize 
    for c = 1 : s(2)                                % columns
        for i = 2 : (size(mat1,1)-1)                 % rows
            j = i-1;
            csd1(j,c) = 2 * mat1 (i,c) - (mat1 (i+1,c) + mat1 (i-1,c));
            csd2(j,c) = 2 * mat2 (i,c) - (mat2 (i+1,c) + mat2 (i-1,c));
        end
    end
else
    s = size (mat);
    mat = [ NaN(1,s(2)); mat; NaN(1,s(2))];         % pad with NaNs
    csd = zeros( size( mat, 1 ), size( mat, 2 ) );  % intialize 
    for c = 1 : s(2)                                % columns
        for i = 2 : (size(mat,1)-1)                 % rows
            j = i-1;
            csd(j,c) = 2 * mat (i,c) - (mat (i+1,c) + mat (i-1,c));
        end
    end
end
return


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
