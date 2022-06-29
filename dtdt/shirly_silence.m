% example:
% cvs_location = '/home/amir/Downloads/_DTDT.xlsx';
% sheet = 1;
% range1= 'MG12:MJ111'
% [pv,pv2,pv3,pv4] = shirly_silence('cvs_location',cvs_location,'sheet',sheet,'range1',range1);

% Shirly silence
function [pv,pv2,pv3,pv4] = shirly_silence(varargin)

[graphics,cvs_location,sheet,range1,issileven] = ParseArgPairs(...
    { 'graphics','cvs_location','sheet','range1','issileven'}...
    , { 0,[],[],[],1}...
    , varargin{ : } );
if ~isempty(cvs_location) && ~isempty(sheet) && ~isempty(range1)
    try
        [NUM,TXT,mat1] = xlsread(cvs_location,sheet,range1);
        NUM1 = NUM(:,1);
        NUM1(isnan(NUM1)) = 0;
        NUM2 = NUM(:,2);
        NUM2(isnan(NUM2)) = 0;
        NUM3 = NUM(:,3);
        NUM3(isnan(NUM3))=0;
        NUM4 = NUM(:,4);
        NUM4(isnan(NUM4)) = 0;
        mat = [NUM1+NUM3 NUM2+NUM4]; 
    catch
        error('could not open csv file')
    end
else

mat = [1	1		
		1	0
		1	1
1	1		
		1	0
1	1		
		1	1
		1	0
		1	0
		1	0
		1	0
		1	1
1	1		
		1	0
		1	0
		1	0
		1	0
		1	1
1	0		
1	1		
		1	1
		1	0
		1	0
		1	0
		1	1
1	1		
		1	0
		1	0
		1	1
1	1		
1	1		
		1	1
		1	0
		1	1
		1	0
1	1		
1	1		
1	1		
		1	1
		1	0
		1	0
		1	1
1	1		
1	1		
1	0		
		1	1
		1	0
		1	1
		1	1
1	1		
		1	1
		1	1
1	0		
1	1		
		1	0
		1	1
1	1		
		1	1
		1	0
1	0		
		1	1
		1	0
1	1		
1	0		
1	0		
1	0		
1	1		
		1	1
1	0		
1	0		
1	1		
1	1		
		1	1
1	1		
		1	1
		1	0
1	0		
1	1		
1	1		
1	1		
		1	0
1	1		
		1	1];
end

Block_size = 12;
numblock = ceil(length(mat)/Block_size);
if issileven
    block_reg = 1:2:numblock;
    block_sil = 2:2:numblock;
else
    block_reg = 2:2:numblock;
    block_sil = 1:2:numblock;
end

missing_trails = Block_size*numblock-length(mat);
blobk_part = true(Block_size,2);
num_reg = length(block_reg);
num_sil = length(block_sil);
reg_mat=[];
sil_mat = [];
for i = 1:numblock
    if sum(ismember(block_reg,i))
    reg_mat = [reg_mat ; blobk_part];
    sil_mat = [sil_mat ; ~blobk_part];
    else
    reg_mat = [reg_mat ; ~blobk_part];
    sil_mat = [sil_mat ; blobk_part];
    end    
end


  reg_mat = logical(reg_mat(1:(end-missing_trails),:));
  sil_mat = logical(sil_mat(1:(end-missing_trails),:));
    

reg_mat2 = mat(reg_mat(:,1),2);
sil_mat2 = mat(sil_mat(:,1),2);
mat = [length(reg_mat2) sum(reg_mat2);length(sil_mat2) sum(sil_mat2)];
% is there a significant difference between the reg and silence
pv = gtest(mat);

% is the reg part significant:
pv2 = myBinomTest( sum(reg_mat2), length(reg_mat2), 0.5, 'one');

% is the sil part significant:
pv3 = myBinomTest( sum(sil_mat2), length(sil_mat2), 0.5, 'one');

%is the total session significant:
pv4 = myBinomTest( sum(sil_mat2)+sum(reg_mat2), length(sil_mat2)+length(reg_mat2), 0.5, 'one');

% plot by block:
[ bino_ci_exact, bino_ci_norm, bino_se_norm ] = binomial_inlines;
err_reg = bino_se_norm( sum(reg_mat2),length(reg_mat2) )';
err_sil = bino_se_norm( sum(sil_mat2),length(sil_mat2) )';

if graphics
figure;
barwerror([1 2],[sum(reg_mat2)/length(reg_mat2) sum(sil_mat2)/ length(sil_mat2)],[err_reg err_sil]',[0 0 1],0.5)
alines(0.5, 'y')
set(gca,'XTickLabel',{'no-silencing' 'silencing'})
ylabel('Success rate')
location_astrix = sum(reg_mat2)/length(reg_mat2) + err_reg + 0.05;
text(1,location_astrix,'*')
% figname = '/amir1/MATLAB/silencing_simcky';
% save_print(figname)
end
return
