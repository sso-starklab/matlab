%load data 

spk=readspk('F:\dat\mV99_03\mV99_03.spk.1',32,32);
clu=load('F:\dat\mV99_03\mV99_03.clu.1');

% delete numbers of clu

clu(1)
clu(1)=[];
vec=unique(clu);
positive=[];

for i =1:length(vec);
%mean WF
idx=(clu==vec(i));
spk5=spk(:,:,idx);
mWF=mean(spk5,3);


[m j]=max(abs(mWF(:,16)));
check=mWF(j,16);

if check<0
    
[minV Imin]=min((mWF(:,16)));
[maxV Imax]=max(mWF(Imin,17:32));
figure;
plot(mWF(Imin,:));hold on;
plot([16 16+Imax],[minV maxV],'o');
title(['unit',num2str(vec(i))]);

else


%for EPS
[maxV Imax]=max((mWF(:,16)));
[minV Imin]=min(mWF(Imax,17:32));
figure;
plot(mWF(Imax,:));hold on;
plot([16 16+Imin],[maxV minV],'o')
title(['unit',num2str(vec(i))])
positive=[positive, vec(i)];
figure;
plot(mWF');
title(['unit',num2str(vec(i))]);
end

end

idx=(clu==9);
spk9=spk(:,:,idx);
mWF=mean(spk9,3);
[maxV Imax]=max((mWF(:,16)));
[minV Imin]=min(mWF(Imax,17:32));
%param wanted
amp = maxV-minV;
realamp = amp * 2.45 * 10^6 / 192 / 2^16;

%dur=Imax/20; %in ms
dur=Imin/20; %in ms

%for linear probes
filebase='F:\dat\mK01_14\mK01_14';
rezfname='F:\dat\mK01_14\rez2.mat'; 

rez=phy2rez(filebase,rezfname);

filebase='F:\dat\EPS\mV99_11\mV99_11';
rezfname='F:\dat\mK01_14\rez3.mat';

rez=ks2ndm(filebase,rezfname);
