clear actual_peaks posm velm dptorquem

%----Data File Signal Channels----%
PLAT_DP_ENC_POS_SIG = 5;
FOOT_GON_POS_SIG = 13;
PLAT_GON_POS_SIG = 14;

PERT_TORQUE_SIG = 7;
IE_TORQUE_SIG = 8; 

F1_SIG = 9;
F2_SIG = 10;
F3_SIG = 11;
F4_SIG = 12;

WEIGHT_SIG = 17;
COP_TORQUE_SIG = 19;

NUM_OF_DATA_BLOCKS = 3;

%Chunks intial data from the data files for each trial.
TRIAL_WINDOW_PRE_PERT  = -380;
TRIAL_WINDOW_POST_PERT = 1020;

shift=0;

time=0;
fclose('all')
% d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',15,'DesignMethod','butter','Samplerate',2000);
            d1 = designfilt('lowpassiir','FilterOrder',2,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
%             d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 15, 'SampleRate', 2000, 'DesignMethod', 'window');
%              d2 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 15, 'SampleRate', 2000, 'DesignMethod', 'window');
d2 = designfilt('lowpassiir','FilterOrder',2,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
%   d3 = designfilt('lowpassfir','FilterOrder',12,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d3 = designfilt('lowpassfir', 'FilterOrder', 12, 'CutoffFrequency', 5, 'SampleRate', 2000, 'DesignMethod', 'window');


%----Find Platform Dynamics Using----%
for l=1
    
    switch l
        case 1
            h = fopen(strcat(DATA_FOLDER_REL_LOC,pfile));
            live_data=fread(h);
            Input1= SimulinkRealTime.utils.getFileScopeData(live_data); 
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





csize=min(size(count),10);
for i=1:csize
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        temp(ran)=((ran-400)/2);
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        dptorque(i,ran)=(emg_data{1,count(i,2)}.data(l+shift,7));
        
        ran=ran+1;
    end
    %     dptorque(i,:)=filtfilt(d2,dptorque(i,:));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        ietorque(i,ran)=(emg_data{1,count(i,2)}.data(l,8));
        
        ran=ran+1;
    end
    ran=1;
    meanpos=mean(emg_data{1,1}.data(:,16));
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        pos(i,ran)=(emg_data{1,count(i,2)}.data(l+221,16)-meanpos)*DP_plat_gonio*pi/180;
      
        ran=ran+1;
    end
    
    vel(i,1)=0;
    ran=2;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+1019
        vel(i,ran)=(pos(i,ran)-pos(i,ran-1))/0.0005;
        ran=ran+1;
    end
    acc(i,1)=0;
    ran=2;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+1019
        acc(i,ran)=(vel(i,ran)-vel(i,ran-1))/0.0005;
        ran=ran+1;
    end
    offsetdptorque(i)=mean(dptorque(i,1:200));
    offsetietorque(i)=mean(ietorque(i,1:200));
    dptorque(i,:)=dptorque(i,:)-offsetdptorque(i);
    ietorque(i,:)=ietorque(i,:)-offsetietorque(i);
    pos(i,:)=pos(i,:)-pos(i,380);
end

posp=mean(pos);
dptorquep=mean(dptorque);
ietorquep=mean(ietorque);
velp=mean(vel);
accp=mean(acc);
fclose('all')
% clear dptorque ietorque actual_peaks count  vel acc velm  temp pos
%%



%----Find Impedance from Subject Data----%
for block_index=1:NUM_OF_DATA_BLOCKS
    
    
    h = fopen(strcat(DATA_FOLDER_REL_LOC,efile,num2str(block_index),'.DAT'));
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    siz=size(Input1.data);

    Input1.data(:,19)=Input1.data(:,19)-time;
            
           
    block_data{block_index}=Input1;
    block_data{1,block_index}.data(:,PLAT_DP_ENC_POS_SIG)=filtfilt(d1,block_data{1,block_index}.data(:,PLAT_DP_ENC_POS_SIG));
    
    block_data{1,block_index}.data(:,FOOT_GON_POS_SIG)=filtfilt(d1,block_data{1,block_index}.data(:,FOOT_GON_POS_SIG));
    block_data{1,block_index}.data(:,PLAT_GON_POS_SIG)=filtfilt(d1,block_data{1,block_index}.data(:,PLAT_GON_POS_SIG));
    block_data{1,block_index}.data(:,PERT_TORQUE_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,PERT_TORQUE_SIG));
    block_data{1,block_index}.data(:,IE_TORQUE_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,IE_TORQUE_SIG));
    block_data{1,block_index}.data(:,F1_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F1_SIG));
    block_data{1,block_index}.data(:,F2_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F2_SIG));
    block_data{1,block_index}.data(:,F3_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F3_SIG));
    block_data{1,block_index}.data(:,F4_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F4_SIG));
    block_data{1,block_index}.data(:,WEIGHT_SIG)=filtfilt(d3,block_data{1,block_index}.data(:,WEIGHT_SIG));
    block_data{1,block_index}.data(:,COP_TORQUE_SIG)=filtfilt(d3,block_data{1,block_index}.data(:,COP_TORQUE_SIG));
    
end

actual_peaks=find_all_peaks(block_data,1,3,dir);
sizet=size(actual_peaks);
k=1;

for t=1:sizet(1,1)
    
    if( actual_peaks(t,2)==signal)
        
        count(k,1)=actual_peaks(t,3);
        count(k,2)=actual_peaks(t,1);
        time_peaks(k,1)=actual_peaks(t,1);
        time_peaks(k,2)=block_data{1,actual_peaks(t,1)}.data(actual_peaks(t,3),COP_TORQUE_SIG);
        
        k=k+1;
    end
end


csize=size(count);
%%
% abc=figure

csize=min(size(count),30);
for i=1:csize
    
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        temp(ran)=((ran-400)/2);
        
        ran=ran+1;
    end
     ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        pos2(i,ran)=(block_data{1,count(i,2)}.data(l,PLAT_DP_ENC_POS_SIG));
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        dptorque(i,ran)=(block_data{1,count(i,2)}.data(l+shift,PERT_TORQUE_SIG));
        a(i,ran)=(block_data{1,count(i,2)}.data(l+shift,COP_TORQUE_SIG));
        ran=ran+1;
    end
    a(i,:)=filtfilt(d3,a(i,:));
    %     dptorque(i,:)=filtfilt(d2,dptorque(i,:));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        ietorque(i,ran)=(block_data{1,count(i,2)}.data(l,IE_TORQUE_SIG));
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        weight1(i,ran)=(block_data{1,count(i,2)}.data(l,WEIGHT_SIG));
        
        ran=ran+1;
    end
    mweight1=mean(weight1(i,280:380));
    ran=1;
    
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        
        cop(i,ran)=(a(i,ran)/mweight1);
        ran=ran+1;
    end
    
    ran=1;
    
    meanpos=mean(block_data{1,1}.data(:,FOOT_GON_POS_SIG));
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        pos(i,ran)=(block_data{1,count(i,2)}.data(l+221,FOOT_GON_POS_SIG)-meanpos)*DP_foot_gonio*pi/180;
        ran=ran+1;
    end
    pos(i,:)=pos(i,:)-mean(pos(i,1:380));
    vel(i,1)=0;
    ran=2;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+1019
        vel(i,ran)=(pos(i,ran)-pos(i,ran-1))/0.0005;
        ran=ran+1;
    end
    
    acc(i,1)=0;
    ran=2;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+1019
        acc(i,ran)=(vel(i,ran)-vel(i,ran-1))/0.0005;
        ran=ran+1;
    end
    
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        taemg(i,ran)=(block_data{1,count(i,2)}.data(l+96,1));
        ran=ran+1;
    end
    
    taemg(i,:)= abs(taemg(i,:)-off_TA)*100/mvc_ta;
    taemg(i,:)=filtfilt(d3,taemg(i,:));
    
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        solemg(i,ran)=(block_data{1,count(i,2)}.data(l+96,2));
        ran=ran+1;
    end
    
    solemg(i,:)= abs(solemg(i,:)-off_SOL)*100/mvc_sol;
    solemg(i,:)=filtfilt(d3,solemg(i,:));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        plemg(i,ran)=(block_data{1,count(i,2)}.data(l+96,3));
        ran=ran+1;
    end
    
    plemg(i,:)= abs(plemg(i,:)-off_PL)*100/mvc_pl;
    plemg(i,:)=filtfilt(d3,plemg(i,:));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        gcaemg(i,ran)=(block_data{1,count(i,2)}.data(l+96,4));
        ran=ran+1;
    end
    
    gcaemg(i,:)= abs(gcaemg(i,:)-off_GCA)*100/mvc_gca;
    gcaemg(i,:)=filtfilt(d3,gcaemg(i,:));
    ran=1;
    
    offsetdptorque(i)=mean(dptorque(i,300:350));
    offsetietorque(i)=mean(ietorque(i,1:200));
    dptorque(i,:)=dptorque(i,:)-offsetdptorque(i);
    ietorque(i,:)=ietorque(i,:)-offsetietorque(i);
      
end
%%
posm=nanmean(pos);
pos2m=mean(pos2);
am=mean(a);
dptorquem=nanmean(dptorque);
ietorquem=mean(ietorque);
velm=mean(vel);
accm=mean(acc);
taemgm=mean(taemg);
solemgm=mean(solemg);
plemgm=mean(plemg);
gcaemgm=mean(gcaemg);
sweight=mean(weight1(300:400))*50/weight;
emgbase(1,1)=mean(taemgm(300:400));
emgbase(1,3)=mean(solemgm(300:400));
emgbase(1,2)=mean(plemgm(300:400));
emgbase(1,4)=mean(gcaemgm(300:400));
emgbase(2,1)=mean(taemgm(400:700));
emgbase(2,3)=mean(solemgm(400:700));
emgbase(2,2)=mean(plemgm(400:700));
emgbase(2,4)=mean(gcaemgm(400:700));
emgbase(3,1)=mean(taemgm(700:1400));
emgbase(3,3)=mean(solemgm(700:1400));
emgbase(3,2)=mean(plemgm(700:1400));
emgbase(3,4)=mean(gcaemgm(400:1400));
copm=mean(cop);
%%

imp=0;
C=[posm(380:560)' velm(380:560)' accm(380:560)'];
d=dptorquem(380:560)'-dptorquep(380:560)'-0.007*accm(380:560)';



A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1;0 0 1];
B=[0 ;0 ;1000;1000;-0.005;0.02];
imp=lsqlin(C,d,A,B)

vartorque=var(dptorquem(380:560)-dptorquep(380:560));
varimp=var(-(imp(1)*posm(380:560)+imp(2)*velm(380:560)+imp(3)*accm(380:560))+dptorquem(380:560)-dptorquep(380:560));
% plot(goodnessfit(:,1),goodnessfit(:,2));
goodness=100*(1-varimp/vartorque);
alalal=mean(copm(300:380))*100


aaa=dptorquem(350:550)-dptorquep(350:550);
imp=imp';
imp(4)=goodness;

imp=[imp,emgbase(1,:),alalal,sweight]
 fclose('all')
if(plotfig==1)
   fig= figure
    ax1=subplot(3,2,1);
    if(perturb=='D')
        title(strcat('Dorsiflexion    ',etitle));
    end
    if(perturb=='P')
        title(strcat('Plantarflexion    ',etitle));
    end
    hold on
    plot(temp,posm*180/pi,'Color',[0 0 0])
    plot(temp,posp*180/pi,'Color',[1 0 1])
    plot(temp,pos2m,'Color',[0 1 1])
    
    %  plot(temp,(dptorquem-dptorquep)*pi/180,'k--')
    axis([-150 250 -inf inf])
    ylabel('angle')
    ax2=subplot(3,2,3);
    hold on
    plot(temp,velm,'Color',[0 0 0])
    plot(temp,velp,'Color',[1 0 1])
    axis([-150 250 -inf inf])
    ylabel('velocity')
    clear pk pl
    ax3=subplot(3,2,5);
    hold on
    plot(temp,accm,'Color',[0 0 0])
%     [pk,pl]=findpeaks(accp);
%     scatter((pl-400)/2,pk,'k');
%     plot(temp,dptorquep,'c');
%     [pk,pl]=findpeaks(dptorquep);
%     scatter((pl-400)/2,pk,'c');
%     plot(temp,dptorquep,'Color',[1 0.5 0.5]);
    axis([-150 250 -inf inf])
    ylabel('acceleration')
    clear pk pl
    ax4=subplot(3,2,2);
    hold on
    
    plot(temp,dptorquem-dptorquep,'Color',[0 0 0])
    hold on
    plot(temp,dptorquep,'c');
    plot(temp,dptorquem,'y');
    
    %     plot(temp,accm*0.0530,'r');
    axis([-150 250 -inf inf])
    ylabel('dptorque')
    %     imp=regress(ietorquem(350:801)'-ietorquep(350:801)',[ accm(350:801)']);
    %  plot(temp,0.3*accm,'r');
    %     imp=regress(dptorquem(380:560)'-dptorquep(380:560)'-0.007*accm(380:560)',[posm(380:560)' velm(380:560)' ]);
    %      plot(temp,posm*imp(1)+velm*imp(2),'r');
    plot(temp,posm*imp(1)+velm*imp(2)+imp(3)*accm,'r');
    ax5=subplot(3,2,4)
    hold on
    plot(temp,imp(1)*posm,'r');
    plot(temp,imp(2)*velm,'g');
    
    plot(temp,imp(1)*posm+imp(2)*velm+imp(3)*accm,'m');
    plot(temp,imp(3)*accm,'b');
    
    plot(temp,dptorquem-dptorquep,'k');
    axis([-150 250 -inf inf])
    % saveas(abc,strcat('inertia',num2str(6*direction),'_',act),'fig')
    % saveas(abc,strcat('inertia',num2str(6*direction),'_',act),'png')
    
    
   
    
    %       legend(strcat('k=',num2str((imp(1))),',   b=',num2str((imp(2))),',   i=',num2str((imp(3)))),'location','SouthOutside');
    
    legend(strcat('k=',num2str((imp(1))),',   b=',num2str((imp(2))),',   i=',num2str((imp(3)))),'location','SouthOutside');
    
    %  axis([-200 500 -50 50])
    subplot(3,2,6)
    
    plot(0,0,'r');
    plot(0,0,'b');
    legend(strcat('t=',num2str(shift/2),',   goodness=',num2str(goodness)),'location','South');
    
    
    linkaxes([ax1,ax2,ax3,ax4,ax5],'x')
    saveas(fig,strcat(etitle,'_dp2','.jpg'));
end