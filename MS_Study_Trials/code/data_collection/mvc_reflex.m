h = fopen(r)
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
off_TA=nanmean(tsmovavg(Input3.data(1:20000,1),'s',1000,1));
 assignin('base','off_TA',off_TA);
off_SOL=nanmean(tsmovavg(Input3.data(1:20000,2),'s',1000,1));
 assignin('base','off_SOL',off_SOL);
off_PL=nanmean(tsmovavg(Input3.data(1:20000,3),'s',1000,1));
 assignin('base','off_PL',off_PL);
 off_GCA=nanmean(tsmovavg(Input3.data(1:20000,4),'s',1000,1));
 assignin('base','off_GCA',off_GCA);
TA_MVC=max((tsmovavg(abs(Input3.data(:,1)-off_TA),'s',1000,1)));
assignin('base','TA_MVC',TA_MVC);
SOL_MVC=max((tsmovavg(abs(Input3.data(:,2)-off_SOL),'s',1000,1)));
assignin('base','SOL_MVC',SOL_MVC);
PL_MVC=max((tsmovavg(abs(Input3.data(:,3)-off_PL),'s',1000,1)));
assignin('base','PL_MVC',SOL_MVC);
GCA_MVC=max((tsmovavg(abs(Input3.data(:,4)-off_GCA),'s',1000,1)));
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
DPbase=nanmean(tsmovavg(Input3.data(1:20000,7),'s',1000,1));
CoPbase=nanmean(tsmovavg(Input3.data(1:20000,18),'s',1000,1));
 assignin('base','DPbase',DPbase);
 assignin('base','CoPbase',CoPbase);
weight=evalin('base','weight'); %Weight is kept in marp
%controller.marp_data not matlab WS
  aweight=weight/2;
  assignin('base','testweight',aweight);
IEbase=nanmean(tsmovavg(Input3.data(1:20000,8),'s',1000,1));
 assignin('base','IEbase',IEbase);
% weight=nanmean(tsmovavg(Input3.data(1:20000,17),'s',1000,1))
%  weight2=nanmean(tsmovavg(Input3.data(60000:70000,17),'s',1000,1));
 
DP_cop=CoPbase*200/weight;
IE_cop=IEbase*200/weight;
%  assignin('base','weight',weight);
assignin('base','DP_cop',DP_cop);
assignin('base','IE_cop',IE_cop);
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