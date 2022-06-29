function rez = set_cutoff(rez)

ops = rez.ops;
dt = 1/1000;

Nk = max(rez.st3(:,2));

% sort by firing rate first
rez.good = zeros(Nk, 1);
for j = 1:Nk
    ix = find(rez.st3(:,2)==j);        
    ss = rez.st3(ix,1)/ops.fs;
    if numel(ss)==0
        continue;
    end
    
    vexp = rez.st3(ix,4);
    
    Th = ops.Th(1);    
    

    fcontamination = 0.1; % acceptable contamination rate

    rez.est_contam_rate(j) = 1;
    while Th>=ops.Th(2)
        st = ss(vexp>Th);
        if isempty(st)
            Th = Th - .5;
            continue;
        end
        [K, Qi, Q00, Q01, rir] = ccg(st, st, 500, dt);
        Q = min(Qi/(max(Q00, Q01)));
        R = min(rir);
        if Q>fcontamination || R>.05                
           break; 
        else
            if Th==ops.Th(1) && Q<.05
                fcontamination = min(.05, max(.01, Q*2));
            end
            rez.good(j) = 1;
            rez.est_contam_rate(j) = Q;
            Th = Th - .5;
        end        
    end
    Th = Th + .5;
    
    rez.Ths(j) = Th;
    rez.st3(ix(vexp<=Th), 2) = 0;
    
    if rem(j,100)==1
%        fprintf('%d \n', j) 
    end
end
% eliminate spikes with id = 0

% we sometimes get NaNs, why?
rez.est_contam_rate(isnan(rez.est_contam_rate)) = 1;

ix = rez.st3(:,2)==0;
rez.st3(ix, :) = [];
if ~isempty(rez.cProj)
    rez.cProj(ix, :) = [];
    rez.cProjPC(ix, :,:) = [];
end
