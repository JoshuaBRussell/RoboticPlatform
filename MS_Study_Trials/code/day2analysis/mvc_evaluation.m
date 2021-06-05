h = fopen(strcat(DATA_FOLDER_REL_LOC,'MVC.dat'));
live_data=fread(h);
Input3= SimulinkRealTime.utils.getFileScopeData(live_data);
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
for i=1:siz(1,1)
if(or(Input3.data(i,3)>10,Input3.data(i,3)<-10))
Input3.data(i,3)=0;
end
end
for i=1:siz(1,1)
if(or(Input3.data(i,4)>10,Input3.data(i,4)<-10))
Input3.data(i,4)=0;
end
end
for i=1:siz(1,1)
if(or(Input3.data(i,18)>10,Input3.data(i,18)<-10))
Input3.data(i,18)=0;
end
end
for i=1:siz(1,1)
if(or(Input3.data(i,19)>10,Input3.data(i,19)<-10))
Input3.data(i,19)=0;
end
end

for i=1:siz(1,1)
if(or(Input3.data(i,17)>10,Input3.data(i,17)<-10))
Input3.data(i,17)=0;
end
end

off_TA=nanmean(tsmovavg(Input3.data(1:20000,1),'s',1000,1));
 %assignin('base','off_TA',off_TA);
off_SOL=nanmean(tsmovavg(Input3.data(1:20000,2),'s',1000,1));
%assignin('base','off_SOL',off_SOL);
off_PL=nanmean(tsmovavg(Input3.data(1:20000,3),'s',1000,1));
 %assignin('base','off_TA',off_PL);
off_GCA=nanmean(tsmovavg(Input3.data(1:20000,4),'s',1000,1));
%assignin('base','off_SOL',off_SOL);



weight=nanmean(tsmovavg(Input3.data(1:20000,17),'s',1000,1));
% %assignin('base','off_SOL2',off_SOL2);

mvc_ta=max((tsmovavg(abs(Input3.data(:,1)-off_TA),'s',1000,1)));

mvc_sol=max((tsmovavg(abs(Input3.data(:,2)-off_SOL),'s',1000,1)));

mvc_pl=max((tsmovavg(abs(Input3.data(:,3)-off_PL),'s',1000,1)));

mvc_gca=max((tsmovavg(abs(Input3.data(:,4)-off_GCA),'s',1000,1)));








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
 IEbase=nanmean(tsmovavg(Input3.data(1:20000,8),'s',1000,1));
 Ts=0.0005;
fclose('all')
figure
subplot(4,1,1)
plot((tsmovavg(abs(Input3.data(:,1)-off_TA),'s',1000,1)));
subplot(4,1,2)
plot((tsmovavg(abs(Input3.data(:,2)-off_SOL),'s',1000,1)));
subplot(4,1,3)
plot((tsmovavg(abs(Input3.data(:,3)-off_PL),'s',1000,1)));
subplot(4,1,4)
plot((tsmovavg(abs(Input3.data(:,4)-off_GCA),'s',1000,1)));