function [mspk,sspk] = addChanels(mspk, sspk, nchannels)

ch2add           = 8 - nchannels;
nsamples         = size(mspk, 2);
nunits           = size(mspk, 3);
for i = 1:ch2add
    n_ch         = size(mspk, 1);
    for j = 1: nunits 
        z            = sspk(:,:,j); 
        mean_std     = mean(z(:));
        new_std      = repelem(mean_std/5, nsamples);
        new_ch       = (mean_std/5) * randn(1, nsamples);
        [b,~] = butter(3,0.3);
        new_ch = filter(b,1,new_ch);
        mspk(n_ch + 1,:,j) = new_ch;
        sspk(n_ch + 1,:,j) = new_std;
    end
end
end
