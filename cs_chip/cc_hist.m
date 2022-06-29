foldername = '/media/shirly/Data/_Shirly/Lab/uLED/fig8/correlation_for_review/';
sess = [13 17 21 23 25 26];
seq = 1:5;

mat = NaN (length(seq) ,length(sess));
ntrials = mat;
ses = mat;
seqs = mat;
pval0 = mat;
l=0;
for i = 1:length(sess)
    for j= 1:length(seq)
        try
            L = load ( sprintf( '%smC400_%d_seq%d_cc.mat', foldername, sess(i), seq(j)));
            mat(j,i)          = mean(L.cc); 
            ntrials (j,i)     = length (L.cc); %number of trials
            l=l+1;
            seq_struct.cc{l}  = L.cc;
            ses(j,i)          = sess(i);
            seqs(j,i)         = seq(j);
            pval0(j,i)        = signrank( L.cc, 0 ); 
        catch
        end
    end
end

mat1 = mat(:);
ntrials1 = ntrials(:);
ses1 = ses(:);
seqs1 = seqs(:);
pval01 = pval0(:);
seq_struct.cc = seq_struct.cc';

nanidx = isnan(mat1);
nunits = [3 3 6 6 7 7 7 7 6 4 6 5 5 7 7];
seq_struct.session = ses1(~nanidx);
seq_struct.seqnum = seqs1(~nanidx);
seq_struct.nunits = nunits';
seq_struct.pval   = pval01(~nanidx);
seq_struct.ntrials = ntrials1(~nanidx);



mmat = nanmean(mat1);
semmat = calc_sem(mat1);

figure; histogram(mat1,10)
alines(mmat,'x','LineStyle','--')
text(mmat,2.9,sprintf('%0.2g+-%0.2g',mmat,semmat))
set( gca, 'tickdir', 'out', 'box', 'off' )
xlabel( 'mean correlation coefficient' )
ylabel( 'Number of sessions' )