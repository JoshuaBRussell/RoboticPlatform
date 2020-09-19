if(perturb=='D' || perturb=='P')
h = fopen(strcat('DPFOOT','.dat'));
live_data=fread(h);
Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);

% 3 is channel for foot and 6 for platform 
gonio=filtfilt(d1,Input1.data(:,4));
gonio=detrend(gonio,'Linear');
plot(gonio)
[a,b]=findpeaks(gonio);
encoder=filtfilt(d1,Input1.data(:,1));

[c,d]=findpeaks(encoder);
for i=1:length(gonio)-221
c_gonio(i)=gonio(i+221);
end
c_gonio=c_gonio-mean(c_gonio);
DP_foot_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
figure
plot(DP_foot_gonio*c_gonio)
hold on
plot(encoder)

% Platform goniometer
gonio=filtfilt(d1,Input1.data(:,6));
plot(gonio)
[a,b]=findpeaks(gonio);
encoder=filtfilt(d1,Input1.data(:,1));

[c,d]=findpeaks(encoder);
for i=1:length(gonio)-221
c_gonio(i)=gonio(i+221);
end
c_gonio=c_gonio-mean(c_gonio);
DP_plat_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
figure
plot(DP_plat_gonio*c_gonio)
hold on
plot(encoder)

else
h = fopen(strcat('IEFOOT','.dat'));
live_data=fread(h);
Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);

% 3 is channel for foot and 6 for platform 
gonio=filtfilt(d1,Input1.data(:,4));
gonio=detrend(gonio,'Linear');
plot(gonio)
[a,b]=findpeaks(gonio);
encoder=filtfilt(d1,Input1.data(:,2));

[c,d]=findpeaks(encoder);
for i=1:length(gonio)-221
c_gonio(i)=gonio(i+221);
end
c_gonio=c_gonio-mean(c_gonio);
IE_foot_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
figure
plot(IE_foot_gonio*c_gonio)
hold on
plot(encoder)

%% Platform goniometer
gonio=filtfilt(d1,Input1.data(:,6));
plot(gonio)
[a,b]=findpeaks(gonio);
encoder=filtfilt(d1,Input1.data(:,2));

[c,d]=findpeaks(encoder);
for i=1:length(gonio)-221
c_gonio(i)=gonio(i+221);
end
c_gonio=c_gonio-mean(c_gonio);
IE_plat_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
figure
plot(IE_plat_gonio*c_gonio)
hold on
plot(encoder)




end
pause(3);
close all
