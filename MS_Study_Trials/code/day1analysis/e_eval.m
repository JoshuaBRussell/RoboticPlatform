

shift=0;

time=0;
fclose('all')
 d1 = designfilt('lowpassiir','FilterOrder',2,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
 d2 = designfilt('lowpassiir','FilterOrder',2,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
%   d3 = designfilt('bandstopfir','FilterOrder',20,'CutoffFrequency1',10,'CutoffFrequency2',30,'SampleRate',2000);
d3 = designfilt('lowpassfir', 'FilterOrder', 12, 'CutoffFrequency', 3, 'SampleRate', 2000, 'DesignMethod', 'window');

for l=1:1
    
    switch l
        case 1
            h = fopen(strcat(pfile,'.DAT'));
            live_data=fread(h);
            Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
            siz=size(Input1.data);
            for i=1:siz(1,1)
                if(or(Input1.data(i,2)>10,Input1.data(i,2)<-10))
                    Input1.data(i,2)=0;
                end
            end
            
            for i=1:siz(1,1)
                if(or(Input1.data(i,1)>10,Input1.data(i,1)<-10))
                    Input1.data(i,1)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,3)>10,Input1.data(i,3)<-10))
                    Input1.data(i,3)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,4)>10,Input1.data(i,4)<-10))
                    Input1.data(i,4)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,5)>10,Input1.data(i,5)<-10))
                    Input1.data(i,5)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,6)>10,Input1.data(i,6)<-10))
                    Input1.data(i,6)=0;
                end
            end
                    
            end
    emg_data{l}=Input1;
  emg_data{1,l}.data(:,5)=filtfilt(d1,emg_data{1,l}.data(:,5));
  emg_data{1,l}.data(:,16)=filtfilt(d1,emg_data{1,l}.data(:,16));
      emg_data{1,l}.data(:,7)=filtfilt(d2,emg_data{1,l}.data(:,7));
      emg_data{1,l}.data(:,8)=filtfilt(d2,emg_data{1,l}.data(:,8));
%       emg_data{1,l}.data(:,8)=filtfilt(d3,emg_data{1,l}.data(:,8));

end

actual_peaks=find_all_peaks(emg_data,1,1,dir);
sizet=size(actual_peaks);
k=1;

for t=1:sizet(1,1)

    if( actual_peaks(t,2)==signal)
        
        count(k,1)=actual_peaks(t,3);
        count(k,2)=actual_peaks(t,1);
    k=k+1;
    end
end
csize=size(count);





csize=min(size(count),12);
for i=1:csize
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        temp(ran)=((ran-400)/2);
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        dptorque(i,ran)=(emg_data{1,count(i,2)}.data(l+shift,7));
        
        ran=ran+1;
    end
%     dptorque(i,:)=filtfilt(d2,dptorque(i,:));
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        ietorque(i,ran)=(emg_data{1,count(i,2)}.data(l+shift,8));
        
        ran=ran+1;
    end
    ran=1;
  meanpos=mean(emg_data{1,1}.data(:,16));
    for l=count(i,1)-380:count(i,1)+1020
        pos(i,ran)=(emg_data{1,count(i,2)}.data(l+221,16)-meanpos)*IE_plat_gonio*pi/180;
        ran=ran+1;
    end
    
    vel(i,1)=0;
    ran=2;
    for l=count(i,1)-380:count(i,1)+1019
        vel(i,ran)=(pos(i,ran)-pos(i,ran-1))/0.0005;
        ran=ran+1;
    end
    acc(i,1)=0;
    ran=2;
    for l=count(i,1)-380:count(i,1)+1019
        acc(i,ran)=(vel(i,ran)-vel(i,ran-1))/0.0005;
        ran=ran+1;
    end
    offsetdptorque(i)=mean(dptorque(i,1:200));
    offsetietorque(i)=mean(ietorque(i,340:350));
    dptorque(i,:)=dptorque(i,:)-offsetdptorque(i);
    ietorque(i,:)=ietorque(i,:)-offsetietorque(i);

end

posp=mean(pos);
dptorquep=mean(dptorque);
ietorquep=mean(ietorque);
velp=mean(vel);
accp=mean(acc);
fclose('all')

 clear dptorque ietorque actual_peaks count  vel acc velm  temp pos
 %%
 for l=1:3
    
   
             h = fopen(strcat(efile,num2str(l),'.DAT'));
            live_data=fread(h);
            Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
            siz=size(Input1.data);
            for i=1:siz(1,1)
                if(or(Input1.data(i,2)>10,Input1.data(i,2)<-10))
                    Input1.data(i,2)=0;
                end
            end
            
            for i=1:siz(1,1)
                if(or(Input1.data(i,1)>10,Input1.data(i,1)<-10))
                    Input1.data(i,1)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,3)>10,Input1.data(i,3)<-10))
                    Input1.data(i,3)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,4)>10,Input1.data(i,4)<-10))
                    Input1.data(i,4)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,5)>10,Input1.data(i,5)<-10))
                    Input1.data(i,5)=0;
                end
            end
            for i=1:siz(1,1)
                if(or(Input1.data(i,6)>10,Input1.data(i,6)<-10))
                    Input1.data(i,6)=0;
                end
            end
        for i=1:siz(1,1)
                if((Input1.data(i,18)>0))
                    time=Input1.data(i,19);
                    start_point=i;
                    break;
                end
        end
        Input1.data(:,19)=Input1.data(:,19)-time;     
    
        emg_data1{l}=Input1;
        emg_data1{1,l}.data(:,5)=filtfilt(d1,emg_data1{1,l}.data(:,5));
        emg_data1{1,l}.data(:,16)=filtfilt(d1,emg_data1{1,l}.data(:,16));
        emg_data1{1,l}.data(:,14)=filtfilt(d1,emg_data1{1,l}.data(:,14));
      emg_data1{1,l}.data(:,7)=filtfilt(d2,emg_data1{1,l}.data(:,7));
      emg_data1{1,l}.data(:,8)=filtfilt(d2,emg_data1{1,l}.data(:,8));
      emg_data1{1,l}.data(:,9)=filtfilt(d2,emg_data1{1,l}.data(:,9));
      emg_data1{1,l}.data(:,10)=filtfilt(d2,emg_data1{1,l}.data(:,10));
      emg_data1{1,l}.data(:,11)=filtfilt(d2,emg_data1{1,l}.data(:,11));
      emg_data1{1,l}.data(:,12)=filtfilt(d2,emg_data1{1,l}.data(:,12));
      
      emg_data1{1,l}.data(:,17)=filtfilt(d3,emg_data1{1,l}.data(:,17));
    emg_data1{1,l}.data(:,19)=filtfilt(d3,emg_data1{1,l}.data(:,19));
%        emg_data1{1,l}.data(:,8)=filtfilt(d3,emg_data1{1,l}.data(:,8));

end

actual_peaks=find_all_peaks(emg_data1,1,3,dir);
sizet=size(actual_peaks);
k=1;

for t=1:sizet(1,1)

    if( actual_peaks(t,2)==signal)
        
        count(k,1)=actual_peaks(t,3);
        count(k,2)=actual_peaks(t,1);
             time_peaks(k,1)=actual_peaks(t,1);
        time_peaks(k,2)=emg_data1{1,actual_peaks(t,1)}.data(actual_peaks(t,3),19);
   
    k=k+1;
    end
end
% sin{1}=load(strcat('DPC501','.csv'));
% sin{2}=load(strcat('DPC502','.csv'));
% 
% 
% for i=1:2
%     
%     foot{i}=[sin{1,i}(:,2),sin{1,i}(:,3),sin{1,i}(:,4),sin{1,i}(:,5)];
%     shank{i}=[sin{1,i}(:,2),sin{1,i}(:,10),sin{1,i}(:,11),sin{1,i}(:,12)];
% end
% 
% % 
%   time_peaks(:,2)=(floor(time_peaks(:,2)*100)/100)+0.01;
% siz11=size(time_peaks);
%  for i=1:siz11(1,1)
%                 
%     siz=size(foot{1,time_peaks(i,1)})
%             for val=1:siz(1,1)
%                 if(abs(time_peaks(i,2)-foot{1,time_peaks(i,1)}(val,1))<0.005)
%                     value_peaks(i,1)=time_peaks(i,1);
%                     value_peaks(i,2)=val;
%                     break;
%                 end
%             end
%             
% end

csize=size(count);

% abc=figure

csize=min(size(count),30);
for i=1:csize
    
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        temp(ran)=((ran-400)/2);
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        dptorque(i,ran)=(emg_data1{1,count(i,2)}.data(l+shift,7));

           a(i,ran)=(emg_data1{1,count(i,2)}.data(l+shift,8));
        ran=ran+1;
    end
%     dptorque(i,:)=filtfilt(d2,dptorque(i,:));
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        ietorque(i,ran)=(emg_data1{1,count(i,2)}.data(l+shift,8));
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        weight_ie(i,ran)=(emg_data1{1,count(i,2)}.data(l,17));
        
        ran=ran+1;
    end
     mweight1=mean(weight_ie(i,280:380));
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
    
       cop(i,ran)=(a(i,ran))/mweight1;
    ran=ran+1;
    end
    
    ran=1;
    
    meanpos=mean(emg_data1{1,1}.data(:,14));
    for l=count(i,1)-380:count(i,1)+1020
        pos(i,ran)=(emg_data1{1,count(i,2)}.data(l+221,14)-meanpos)*IE_foot_gonio*pi/180;
        ran=ran+1;
    end
    pos(i,:)=pos(i,:)-mean(pos(i,1:380));
    ran=1; 
    meanpos2=mean(emg_data1{1,1}.data(:,15));
   
    for l=count(i,1)-380:count(i,1)+1020
        pos2(i,ran)=(emg_data1{1,count(i,2)}.data(l+221,16)-meanpos2)*IE_plat_gonio*pi/180;
        ran=ran+1;
    end
    pos2(i,:)=pos2(i,:)-mean(pos2(i,1:380));
    ran=1;
  
    pos2(i,:)=pos2(i,:)-mean(pos2(i,1:380));
    vel(i,1)=0;
    ran=2;
    for l=count(i,1)-380:count(i,1)+1019
        vel(i,ran)=(pos(i,ran)-pos(i,ran-1))/0.0005;
        ran=ran+1;
    end
        
    vel2(i,1)=0;
    ran=2;
    for l=count(i,1)-380:count(i,1)+1019
        vel2(i,ran)=(pos2(i,ran)-pos2(i,ran-1))/0.0005;
        ran=ran+1;
    end
    
    acc(i,1)=0;
    ran=2;
    for l=count(i,1)-380:count(i,1)+1019
        acc(i,ran)=(vel(i,ran)-vel(i,ran-1))/0.0005;
        ran=ran+1;
    end
       acc2(i,1)=0;
    ran=2;
    for l=count(i,1)-380:count(i,1)+1019
        acc2(i,ran)=(vel2(i,ran)-vel2(i,ran-1))/0.0005;
        ran=ran+1;
    end 
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        weight1(i,ran)=(emg_data1{1,count(i,2)}.data(l,17));
        
        ran=ran+1;
    end
    mweight1=mean(weight1(i,280:380));
    
     ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        taemg(i,ran)=(emg_data1{1,count(i,2)}.data(l+96,1));
        ran=ran+1;
    end
    
    taemg(i,:)= abs(taemg(i,:)-off_TA)*100/mvc_ta;
    taemg(i,:)=filtfilt(d3,taemg(i,:));
    
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        solemg(i,ran)=(emg_data1{1,count(i,2)}.data(l+96,2));
        ran=ran+1;
    end
    
    solemg(i,:)= abs(solemg(i,:)-off_SOL)*100/mvc_sol;
    solemg(i,:)=filtfilt(d3,solemg(i,:));
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        plemg(i,ran)=(emg_data1{1,count(i,2)}.data(l+96,3));
        ran=ran+1;
    end
    
    plemg(i,:)= abs(plemg(i,:)-off_PL)*100/mvc_pl;
    plemg(i,:)=filtfilt(d3,plemg(i,:));
    ran=1;
    for l=count(i,1)-380:count(i,1)+1020
        gcaemg(i,ran)=(emg_data1{1,count(i,2)}.data(l+96,4));
        ran=ran+1;
    end
    
    gcaemg(i,:)= abs(gcaemg(i,:)-off_GCA)*100/mvc_gca;
    gcaemg(i,:)=filtfilt(d3,gcaemg(i,:));
    ran=1;
    
 
    offsetdptorque(i)=mean(dptorque(i,1:200));
    offsetietorque(i)=mean(ietorque(i,340:360));
    dptorque(i,:)=dptorque(i,:)-offsetdptorque(i);
    ietorque(i,:)=ietorque(i,:)-offsetietorque(i);
    
%         ax1=subplot(3,2,1);
%         hold on
%         
%         plot(temp,pos(i,:),'Color',[0.7 0.7 0.7])
%         axis([-200 500 -300 300])
%         ax2=subplot(3,2,3);
%         hold on
%         
%         plot(temp,vel(i,:),'Color',[0.7 0.7 0.7])
%         axis([-150 250 -inf inf])
%         ax3=subplot(3,2,5);
%         hold on
%         
%         plot(temp,acc(i,:),'Color',[0.7 0.7 0.7])
%         axis([-150 250 -inf inf])
%         
%         ax4=subplot(3,2,2);
%         hold on
%         
%             plot(temp,dptorque(i,:)-dptorquep,'Color',[0.7 0.7 0.7])
%             axis([-200 500 -100 100])
%           indimp(i,:)=regress(dptorque(i,380:550)'-dptorquep(380:550)',[pos(i,380:550)' vel(i,380:550)' ]);
%         
    
end
% for i =1:10
%         if(pos(i,580)*180/pi<2)
%             pos(i,:)=NaN;
%             ietorque(i,:)=NaN;
%             vel(i,:)=NaN;
%             acc(i,:)=NaN;
%         end
% end
%%
posm=nanmean(pos);
pos2m=trimmean(pos2,30);

am=trimmean(a,30);
dptorquem=trimmean(dptorque,30);
ietorquem=nanmean(ietorque);
velm=trimmean(vel,30);
vel2m=trimmean(vel2,30);
accm=trimmean(acc,30);
acc2m=trimmean(acc2,30);
taemgm=trimmean(taemg,30);
solemgm=trimmean(solemg,30);
plemgm=trimmean(plemg,30);
gcaemgm=trimmean(gcaemg,30);
emgbase(1,1)=mean(taemgm(300:400));
emgbase(1,3)=mean(solemgm(300:400));
emgbase(1,2)=mean(plemgm(300:400));
emgbase(1,4)=mean(gcaemgm(300:400));
sweight=mean(weight1(300:400))*50/weight;
emgbase(2,1)=mean(taemgm(400:700));
emgbase(2,3)=mean(solemgm(400:700));
emgbase(2,2)=mean(plemgm(400:700));
emgbase(2,4)=mean(gcaemgm(400:700));
emgbase(3,1)=mean(taemgm(700:1400));
emgbase(3,3)=mean(solemgm(700:1400));
emgbase(3,2)=mean(plemgm(700:1400));
emgbase(3,4)=mean(gcaemgm(400:1400));
copm=mean(cop);
weight_iem=mean(weight_ie);
%%
imp=0;
% load('edynamics.mat')
k=ietorquem-ietorquep;
 C=[posm(380:550)' velm(380:550)' accm(380:550)'];
d=ietorquem(380:550)'-ietorquep(380:550)';
d=d-1*min(k(380:450));


A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1;0 0 1];
B=[0 ;0 ;1000;1000;-0.002;0.004];
imp=lsqlin(C,d,A,B)

vartorque=var(ietorquem(380:550)-ietorquep(380:550)-1*min(k(380:450)));
varimp=var(-(imp(1)*posm(380:550)+imp(2)*velm(380:550)+imp(3)*accm(380:550))+ietorquem(380:550)-ietorquep(380:550));
% plot(goodnessfit(:,1),goodnessfit(:,2));
goodness=100*(1-varimp/vartorque);

 alalal=mean(copm(250:380))*100
 
 
 aaa=ietorquem(350:550)-ietorquep(350:550);
 imp=imp';
imp(4)=goodness;

imp=[imp,emgbase(1,:),alalal,sweight];

fclose('all')
if(plotfig==1)
   fig= figure
ax1=subplot(3,2,1);
    title(strcat('Eversion    ',etitle));

hold on
plot(temp,posm*180/pi,'Color',[0 0 0])
plot(temp,pos2m*180/pi,'Color',[1 0 1])
%  plot(temp,(dptorquem-dptorquep)*pi/180,'k--')
axis([-150 250 -inf inf])
ylabel('angle')
ax2=subplot(3,2,3);
hold on
plot(temp,velm,'Color',[0 0 0])
% plot(temp,velp,'Color',[1 0 1])

axis([-150 250 -inf inf])
ylabel('velocity')
clear pk pl
ax3=subplot(3,2,5);
hold on
plot(temp,weight_iem);
% plot(temp,accm,'Color',[0 0 0])
% [pk,pl]=findpeaks(acc2m);
%  scatter((pl-400)/2,pk,'k');
%  plot(temp,ietorquep,'c');
%  [pk,pl]=findpeaks(ietorquep);
%  scatter((pl-400)/2,pk,'c');
% plot(temp,ietorquep,'Color',[1 0.5 0.5]);
 axis([-150 250 -inf inf])
ylabel('weight')
clear pk pl

ax4=subplot(3,2,2);
hold on

    plot(temp,ietorquem-ietorquep-1*min(k(380:450)),'Color',[0 0 0]);
    hold on
     plot(temp,ietorquep,'c');
     plot(temp,ietorquem,'y');
    
%     plot(temp,accm*0.0530,'r');
    axis([-150 250 -inf inf])
    ylabel('ietorque')
%     imp=regress(ietorquem(350:801)'-ietorquep(350:801)',[ accm(350:801)']);
%  plot(temp,0.3*accm,'r');
%     imp=regress(dptorquem(380:550)'-dptorquep(380:550)'-0.007*accm(380:550)',[posm(380:550)' velm(380:550)' ]);
%      plot(temp,posm*imp(1)+velm*imp(2),'r');

  plot(temp,posm*imp(1)+velm*imp(2)+imp(3)*accm,'r');
ax5=subplot(3,2,4)
hold on
 plot(temp,imp(1)*posm,'r');
  plot(temp,imp(2)*velm,'g');

   plot(temp,imp(1)*posm+imp(2)*velm+imp(3)*accm,'m');
           plot(temp,0.003*accm,'b');
    
   plot(temp,ietorquem-ietorquep-1*min(k(380:450)),'k');
axis([-150 250 -inf inf])
% saveas(abc,strcat('inertia',num2str(6*direction),'_',act),'fig')
% saveas(abc,strcat('inertia',num2str(6*direction),'_',act),'png')


fclose('all')

%       legend(strcat('k=',num2str((imp(1))),',   b=',num2str((imp(2))),',   i=',num2str((imp(3)))),'location','SouthOutside');
 
     legend(strcat('k=',num2str((imp(1))),',   b=',num2str((imp(2))),',   i=',num2str((imp(3)))),'location','SouthOutside');
 
%  axis([-200 500 -50 50])
subplot(3,2,6)
plot(0,0,'r');
plot(0,0,'b');
 legend(strcat('t=',num2str(shift/2),',   goodness=',num2str(goodness)),'location','South');
 

 linkaxes([ax1,ax2,ax3,ax4,ax5],'x');
 saveas(fig,strcat(etitle,'_ie','.jpg'));
end