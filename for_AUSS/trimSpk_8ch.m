function [spk,ind] = trimSpk_8ch(spk1)

nchannels      =size(spk1,1);
switch nchannels
    
    case 8
    spk          = spk1;
    ind          =(1:8);
    
    case 9
        [extrimum,I]   =max(abs(spk1(:,16)));
        if I<= 5
            ind        =(1:8)';
            spk       = spk1(ind,:);
           
        else
            ind        =(2:9)';
            spk       = spk1(ind,:);
           
        end
        
    case 10
        [extrimum,I]   =max(abs(spk1(:,16)));
        if I<= 5
            ind        =(1:8)';
            spk       = spk1(ind,:);
           
        else
            ind        =(3:10)';
            spk       = spk1(ind,:);
       
        end
    case 11
        [extrimum,I]   =max(abs(spk1(:,16)));
        if I<= 5
            ind        =(1:8)';
            spk       = spk1(ind,:);
         
        elseif I>6
            ind        =(4:11)';
            spk       = spk1(ind,:);
     
        else
            ind        =(3:10)';
            spk       = spk1(ind,:);

        end
%%%%
%writeNPY(xtag,"/home/tali/matlab/AUSS_python/mP31_04/mP31_04.trim.1-1-1.npy")
%%%    
    
end
