% adding probe type to the sst structure
%       High-density 15um    {0}
%       High-density 20um    {1}
%       linear               {2}
%       Dual-sided           {3}
%       else                 {4}

function sst = addprobedensity (sst)

    % cell array of sessions from high-density of 15um spacing probes
    prb_HD15 = {'mC41','mF79', 'mF93', 'mF105', 'mF108', 'mP20' };
                    
    % cell array of sessions from high-density of 20um spacing probes
    prb_HD20 = {'mA234', 'mP23', 'mP101', 'mDL5', 'mB142', 'mS234'};

    % cell array of sessions from linear probes
    prb_LR       = {'mF84', 'mO251','mV99', 'mK01'};
    % cell array of sessions from Dual-sided probes
    prb_DS      = {  'mDS1', 'mDS2'};
  

for i = 1:length (sst.filebase)    
    if contains(sst.filebase{i}, prb_HD15)
        sst.probeDens(i)       = 0;
    elseif  contains(sst.filebase{i}, prb_HD20)
        sst.probeDens(i)       = 1;
    elseif  contains(sst.filebase{i}, prb_LR) 
        sst.probeDens(i)       = 2;
    elseif  contains(sst.filebase{i}, prb_DS) 
        sst.probeDens(i)       = 3;        
    else
        sst.probeDens(i) = 4;
    end
end
sst.probeDens = sst.probeDens';
return