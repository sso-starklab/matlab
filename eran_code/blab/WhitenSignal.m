%function [y, ARmodel] = WhitenSignal(x, window,CommonAR,ARmodel, ARorder)
%
% whitens the signal 
% if window specified will recompute the model in each window of that size
% (window is in samples,e,g, 300sec*1250 samples
% if CommonAR is set to 1, then will use model from first channel for all
% if ARmodel is specified - use it, not compute fromthe data
% output optionaly the ARmodel for use on the other data to be on the same scale

% revisions
% 04-nov-11 (1) removed addpath (2) replaced Filter0 -> firfilt (zero phase distortion)

function [y, A] = WhitenSignal(x,varargin)

%artype =2; %Signal processing toolbox
artype      = 1; %arfit toolbox
[window,CommonAR, ARmodel,ArOrder] = DefaultArgs(varargin,{[],1,[],1});
ArOrder = ArOrder+1;
Trans = 0;
if size(x,1)<size(x,2)
    x = x';
end
[nT, nCh]  = size(x);
y = zeros(nT,nCh);
if isempty(window)
    seg = [1 nT];
    nwin=1;
else
    nwin = floor(nT/window)+1;
    seg = repmat([1 window],nwin,1)+repmat([0:nwin-1]'*window,1,2);
    if nwin*window>nT
        seg(end,2) =nT;
    end   
end

for j=1:nwin
    if ~isempty(ARmodel) 
        A = ARmodel;
        for i=1:nCh
            y(seg(j,1):seg(j,2),i) = firfilt(x(seg(j,1):seg(j,2),i),A);
        end
    else
        if CommonAR % meaning common model for all channels and segments!!! 
            for i=1:nCh
                if  j==1 && i==1
                    switch artype
                        case 1
                            [ ~, Atmp] = arfit(x(seg(j,1):seg(j,2),i),ArOrder,ArOrder);
                            A = [1 -Atmp];
                        case 2
                            A = arburg(x(seg(j,1):seg(j,2),i),ArOrder);
                    end
                    ARmodel = A;
                end
                y(seg(j,1):seg(j,2),i) = firfilt(x(seg(j,1):seg(j,2),i),A);
            end
        else
            for i=1:nCh
                switch artype
                    case 1
                        [ ~, Atmp] = arfit(x(seg(j,1):seg(j,2),i),ArOrder,ArOrder);
                        A =[1 -Atmp];
                    case 2
                        A = arburg(x(seg(j,1):seg(j,2),i),ArOrder);
                end
                y(seg(j,1):seg(j,2),i) = firfilt(x(seg(j,1):seg(j,2),i), A);
            end
        end
    end
end

if Trans
    y =y';
end

return

