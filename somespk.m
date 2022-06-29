clu=load('F:\dat\mF84_03\mF84_03.clu.1');
clu(1)=[];

vec=unique(clu);
units=[2 4];

idx2 = find( clu == units (1) );
idx4 = find( clu == units (2) );
idxn=[idx2;idx4];

spk2 = readspk( 'F:\dat\mF84_03\mF84_03.spk.1',32,32, idx2);
spk4 = readspk( 'F:\dat\mF84_03\mF84_03.spk.1',32,32, idx4);

mWF=mean(spk4,3);

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