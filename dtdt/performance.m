SucPW = NaN(95,7);
TrlPW = NaN(95,7);
%paste data from DTDT
TrlPW(TrlPW==0)=NaN;
SucPW(SucPW==0)=NaN;

flag335 = ~isnan(SucPW(:,2));
flag79 = ~isnan(SucPW(:,3));
flag271 = ~isnan(SucPW(:,4));
flag109 = ~isnan(SucPW(:,5));
flag352 = ~isnan(SucPW(:,6));
flag353 = ~isnan(SucPW(:,7));


figure, plot(SucPW(flag335,1),SucPW(flag335,2),'LineWidth', 1,'Marker', '.' , 'MarkerSize',10)
hold on,
plot(SucPW(flag79,1),SucPW(flag79,3),'LineWidth', 1, 'Marker', '.' , 'MarkerSize',10)
plot(SucPW(flag271,1),SucPW(flag271,4),'LineWidth', 1, 'Marker', '.' , 'MarkerSize',10)
plot(SucPW(flag109,1),SucPW(flag109,5),'LineWidth', 1,'Marker', '.' , 'MarkerSize',10)
plot(SucPW(flag352,1),SucPW(flag352,6),'LineWidth', 1,'Marker', '.' , 'MarkerSize',10)
plot(SucPW(flag353,1),SucPW(flag353,7),'LineWidth', 1,'Marker', '.' , 'MarkerSize',10)
alines(0.5,'y','lineStyle','--', 'color','k')
alines(0,'x','lineStyle','--', 'color','k')
set( gca, 'tickdir', 'out', 'box', 'off' );
legend ('mA335','m79','m271','m109','m352','m353')

title ('Success rates')

figure, plot(TrlPW(flag335,1),TrlPW(flag335,2),'LineWidth', 1,'Marker', '.' , 'MarkerSize',10)
hold on,
plot(TrlPW(flag79,1),TrlPW(flag79,3),'LineWidth', 1, 'Marker', '.' ,'MarkerSize',10)
plot(TrlPW(flag271,1),TrlPW(flag271,4),'LineWidth', 1, 'Marker', '.' ,'MarkerSize',10)
plot(TrlPW(flag109,1),TrlPW(flag109,5),'LineWidth', 1, 'Marker', '.' ,'MarkerSize',10)
plot(TrlPW(flag352,1),TrlPW(flag352,6),'LineWidth', 1, 'Marker', '.' ,'MarkerSize',10)
plot(TrlPW(flag353,1),TrlPW(flag353,7),'LineWidth', 1, 'Marker', '.' ,'MarkerSize',10)
alines(75,'y','lineStyle','--', 'color','k')
alines(0,'x','lineStyle','--', 'color','k')
set( gca, 'tickdir', 'out', 'box', 'off' );
legend ('mA335','m79','m271','m109','m352','m353')
title ('Number of trials')

flag335 = ~isnan(SucPW2(:,2));
flag79 = ~isnan(SucPW2(:,3));
flag271 = ~isnan(SucPW2(:,4));
flag109 = ~isnan(SucPW2(:,5));
flag352 = ~isnan(SucPW2(:,6));
flag353 = ~isnan(SucPW2(:,7));

num335 = sum(flag335);
num79 = sum(flag79);
num271 = sum(flag271);
num109 = sum(flag109 );
num352 = sum(flag352 );
num353 = sum(flag353 );

SucPWBox = [SucPW2(flag335,2)' SucPW2(flag79,3)' SucPW2(flag271,4)' SucPW2(flag109,5)' SucPW2(flag352,6)' SucPW2(flag353,7)'];
TrlPWBox = [TrlPW(flag335,2)' TrlPW(flag79,3)' TrlPW(flag271,4)' TrlPW(flag109,5)' TrlPW(flag352,6)' TrlPW(flag353,7)'];
grps     = [ones(1,num335) 2*ones(1,num79) 3*ones(1,num271) 4*ones(1,num109) 5*ones(1,num352) 6*ones(1,num353)];

figure,boxplot (SucPWBox, grps,'notch','on','Labels',{'mA335','mP79' ,'mA271','mF109','mA352','mA353'})
alines(0.5,'y','lineStyle','--', 'color','k')
set( gca, 'tickdir', 'out', 'box', 'off' );
title ('Success rates')
for i=1:6
    x = SucPW2(:,i+1);
    pv(i) = signrank( x, 0.5, 'tail', 'right' );
    if pv(i)<0.001
        text(i,0.87,'***')
    elseif pv(i)<0.005
        text(i,0.87,'**')
    elseif pv(i)<0.05
        text(i,0.87,'*')
    end
end
    
figure,boxplot (TrlPWBox, grps,'notch','on','Labels',{'mA335','mP79' ,'mA271','mF109','mA352','mA353'})
set( gca, 'tickdir', 'out', 'box', 'off' );
title ('Number of trials')



m1tr = [60
48
80
168
155
81
72
71
84
92
161
45
71
82
115
107
72
110
96
111
108
142
132
145
66
124
97
144
138];

m2tr = [60
58
68
173
48
116
83
59
72
72
96
68
42
58
95
85
68];

m3tr = [    300
96
168
200
144
156
144
168
132
144
204];

m4tr = [84
84
124
141];

m5tr = [71
84
96
83
108
96
108
72
96
108
93
154
203
180
154
72
204
217
274
76
96];

m6tr = [120
47
96
62
110
120
99
96
72
89
86
96
132
102
108
132
120
152
84
84
100
96
93
64
98
108
91
92
96
88
106
64
80
114];

allmicetr = NaN(34,6);
allmicetr(1:29,1) = m1tr;
allmicetr(1:17,2) = m2tr;
allmicetr(1:11,3) = m3tr;
allmicetr(1:4,4) = m4tr;
allmicetr(1:21,5) = m5tr;
allmicetr(1:34,6) = m6tr;
SucPW_used = SucPW(1:34,2:7);
allmice_succ_tr = round(allmicetr.*SucPW_used);
issig=[];
pvmat=[];
for i = 1:34
    for j = 1:6
        if ~isnan(allmice_succ_tr(i,j))
            pvmat(i,j) = myBinomTest(allmice_succ_tr(i,j),allmicetr(i,j),0.5,'one');
            if pvmat(i,j)<=0.05 && SucPW_used(i,j)>0.5
                issig(i,j) = 1;
            else
                issig(i,j) = 0;
            end
        else
           pvmat(i,j) = NaN;
           issig(i,j) = NaN;
        end
    end
end

for i = 1:6
    %find first succ for each mouse
    sigtemp = issig(:,i);
    aa = find(sigtemp,1,'first');
    bb = sum(~isnan(sigtemp));
    sess_after(i) = bb-aa+1;
    succ_after(i) = sum(sigtemp(aa:bb));
    chance(i) = max(0.05, 1/ sess_after(i));  
    binp = @(obs,tot,p)(1-binocdf(obs-1,tot,p));
    pv_after(i) = binp(succ_after(i), sess_after(i), chance(i));
    [ bino_ci_exact, bino_ci_norm, bino_se_norm ] = binomial_inlines;
    err(i) = bino_se_norm(succ_after(i) ,sess_after(i) );
end

figure,
barwerror([1 2 3 4 5 6],succ_after./sess_after ,err);
alines(chance(1),'y','lineStyle','--', 'color','k') % mA335
alines(chance(2),'y','lineStyle','--', 'color','b') % mP79
alines(chance(3),'y','lineStyle','--', 'color','r') % mA271
alines(chance(4),'y','lineStyle','--', 'color','m') % mF109
alines(chance(5),'y','lineStyle','--', 'color','g') % mA352
alines(chance(6),'y','lineStyle','--', 'color','y') % mA353


