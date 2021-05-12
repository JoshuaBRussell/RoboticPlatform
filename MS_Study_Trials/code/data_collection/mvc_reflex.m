
%This assumes only a single foot is on the platform.


TA_EMG_SIG  = 1;
SOL_EMG_SIG = 2;
PL_EMG_SIG  = 3;
GCA_EMG_SIG = 4;

WEIGHT_SIG = 17; 

DP_COP_TORQUE_SIG = 18;
IE_COP_TORQUE_SIG = 8;
DP_TORQUE_SIG = 7;

h = fopen('MVC.DAT')
live_data2=fread(h);
Input3= SimulinkRealTime.utils.getFileScopeData(live_data2);
siz=size(Input3.data)
siz=size(Input3.data)
for i=1:siz(1,1)
if(or(Input3.data(i,2)>10,Input3.data(i,2)<-10))
Input3.data(i,2)=0;
end
end

for i=1:siz(1,1)
if(or(Input3.data(i,1)>10,Input3.data(i,1)<-10))
Input3.data(i,1)=0;
end
end

%----Find EMG DC Bias----%
off_TA=nanmean(tsmovavg(Input3.data(1:20000,TA_EMG_SIG),'s',1000,1));
off_SOL=nanmean(tsmovavg(Input3.data(1:20000,SOL_EMG_SIG),'s',1000,1));
off_PL=nanmean(tsmovavg(Input3.data(1:20000,PL_EMG_SIG),'s',1000,1));
off_GCA=nanmean(tsmovavg(Input3.data(1:20000,GCA_EMG_SIG),'s',1000,1));

assignin('base','off_TA',off_TA);
assignin('base','off_SOL',off_SOL);
assignin('base','off_PL',off_PL);
assignin('base','off_GCA',off_GCA);

%----Find EMG MVC----%
TA_MVC=max((tsmovavg(abs(Input3.data(:,TA_EMG_SIG)-off_TA),'s',1000,1)));
SOL_MVC=max((tsmovavg(abs(Input3.data(:,SOL_EMG_SIG)-off_SOL),'s',1000,1)));
PL_MVC=max((tsmovavg(abs(Input3.data(:,PL_EMG_SIG)-off_PL),'s',1000,1)));
GCA_MVC=max((tsmovavg(abs(Input3.data(:,GCA_EMG_SIG)-off_GCA),'s',1000,1)));

assignin('base','TA_MVC',TA_MVC);
assignin('base','SOL_MVC',SOL_MVC);
assignin('base','PL_MVC',SOL_MVC);
assignin('base','GCA_MVC',GCA_MVC);

for i=2:siz(1,1)
    if(or(Input3.data(i,7)-Input3.data(i-1,7)>100,Input3.data(i,7)-Input3.data(i-1,7)<-100))
        Input3.data(i,7)=Input3.data(i-1,7);
    end
end
for i=2:siz(1,1)
    if(abs(Input3.data(i,8)-Input3.data(i-1,8))>100)
        Input3.data(i,8)=Input3.data(i-1,8);
    end
end
DPbase=nanmean(tsmovavg(Input3.data(1:20000,DP_TORQUE_SIG),'s',1000,1));
IEbase=nanmean(tsmovavg(Input3.data(1:20000,IE_COP_TORQUE_SIG),'s',1000,1));
CoPbase=nanmean(tsmovavg(Input3.data(1:20000,DP_COP_TORQUE_SIG),'s',1000,1));
single_weight = mean(Input3.data(1:20000, WEIGHT_SIG)); %Gets the base weight while a single foot is on the platform. 

assignin('base','DPbase',DPbase);
assignin('base','IEbase',IEbase);
assignin('base','CoPbase',CoPbase);
assignin('base','single_weight', single_weight);


%This uses the torque/(actual weight) on the platform for a single foot (instead of
%just assuming the weight is evenly distributed between the both legs (i.e. torque/(weight/2))).
DP_cop=mean(Input3.data(1:20000, DP_COP_TORQUE_SIG)./Input3.data(1:20000, WEIGHT_SIG)) * 100; %in cm - x100 converts to cm
IE_cop=mean(Input3.data(1:20000, IE_COP_TORQUE_SIG)./Input3.data(1:20000, WEIGHT_SIG)) * 100; %in cm - x100 converts to cm

assignin('base','DP_cop',DP_cop);
assignin('base','IE_cop',IE_cop);




%----EMG Plot----%
figure
subplot(4,1,1)
plot((tsmovavg(abs(Input3.data(:,1)-off_TA),'s',1000,1)));
title("TA");
subplot(4,1,2)
plot((tsmovavg(abs(Input3.data(:,2)-off_SOL),'s',1000,1)));
title("SOL");
subplot(4,1,3)
plot((tsmovavg(abs(Input3.data(:,3)-off_PL),'s',1000,1)));
title("PL");
subplot(4,1,4)
plot((tsmovavg(abs(Input3.data(:,4)-off_GCA),'s',1000,1)));
title("GCA");
fclose(h);