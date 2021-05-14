
%----Data File Signal Channels----%
PLAT_DP_ENC_POS_SIG = 5;
FOOT_GON_POS_SIG = 13;
PLAT_GON_POS_SIG = 14;

DP_TORQUE_SIG = 7; 
PERT_TORQUE_SIG = 8;

F1_SIG = 9;
F2_SIG = 10;
F3_SIG = 11;
F4_SIG = 12;

WEIGHT_SIG = 17;
COP_TORQUE_SIG = 19;

SAMPLE_RATE_HZ = 2000;
SAMPLE_PERIOD = 1/SAMPLE_RATE_HZ;

NUM_OF_DATA_BLOCKS = 3;

%Chunks intial data from the data files for each trial.
TRIAL_WINDOW_PRE_PERT  = -380;
TRIAL_WINDOW_POST_PERT = 1020;

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
            h = fopen(strcat(DATA_FOLDER_REL_LOC, pfile,'.DAT'));
            live_data=fread(h);
            Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
            siz=size(Input1.data);
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
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        dptorque(i,ran)=(emg_data{1,count(i,2)}.data(l+shift,7));
        
        ran=ran+1;
    end
%     dptorque(i,:)=filtfilt(d2,dptorque(i,:));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        ietorque(i,ran)=(emg_data{1,count(i,2)}.data(l+shift,8));
        
        ran=ran+1;
    end
    ran=1;
  meanpos=mean(emg_data{1,1}.data(:,16));
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        pos(i,ran)=(emg_data{1,count(i,2)}.data(l+221,16)-meanpos)*IE_plat_gonio*pi/180;
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
 for block_index=1:NUM_OF_DATA_BLOCKS
    

        h = fopen(strcat(DATA_FOLDER_REL_LOC,efile,num2str(block_index),'.DAT'));
        live_data=fread(h);
        Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
        siz=size(Input1.data);
            
        for i=1:siz(1,1)
                if((Input1.data(i,18)>0))
                    time=Input1.data(i,WEIGHT_SIG);
                    start_point=i;
                    break;
                end
        end
        Input1.data(:,WEIGHT_SIG)=Input1.data(:,WEIGHT_SIG)-time;     
    
        block_data{block_index}=Input1;
        block_data{1,block_index}.data(:,5)=filtfilt(d1,block_data{1,block_index}.data(:,5));
        block_data{1,block_index}.data(:,FOOT_GON_POS_SIG)=filtfilt(d1,block_data{1,block_index}.data(:,FOOT_GON_POS_SIG));
        block_data{1,block_index}.data(:,PLAT_GON_POS_SIG)=filtfilt(d1,block_data{1,block_index}.data(:,PLAT_GON_POS_SIG));
        block_data{1,block_index}.data(:,DP_TORQUE_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,DP_TORQUE_SIG));
        block_data{1,block_index}.data(:,PERT_TORQUE_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,PERT_TORQUE_SIG));
        block_data{1,block_index}.data(:,F1_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F1_SIG));
        block_data{1,block_index}.data(:,F2_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F2_SIG));
        block_data{1,block_index}.data(:,F3_SIG)=filtfilt(d2,block_data{1,block_index}.data(:,F3_SIG));
        
        temp = block_data{1,block_index}.data(:,F4_SIG);
        temp(isnan(temp))=0
        block_data{1,block_index}.data(:,F4_SIG)=filtfilt(d2,temp);
      
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
        time_peaks(k,2)=block_data{1,actual_peaks(t,1)}.data(actual_peaks(t,3),WEIGHT_SIG);
   
    k=k+1;
    end
end



%Remove any "Dud" perturbations from being indexed
temp_count = [];
for block_index = 1:NUM_OF_DATA_BLOCKS
   block_start_index = min(find(count(:, 2) == block_index)); 
   true_block_indices = count(block_start_index:block_start_index + 10 - 1, :);
   temp_count = [temp_count; true_block_indices];
end

count = temp_count;

csize=size(count);
csize=min(size(count),30);

for i=1:csize
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        %dptorque(i,ran)=(block_data{1,count(i,2)}.data(l+shift,DP_TORQUE_SIG));

        a(i,ran)=(block_data{1,count(i,2)}.data(l+shift,COP_TORQUE_SIG));
        ran=ran+1;
    end
%     dptorque(i,:)=filtfilt(d2,dptorque(i,:));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        ietorque(i,ran)=(block_data{1,count(i,2)}.data(l+shift,PERT_TORQUE_SIG));
        
        ran=ran+1;
    end
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        weight_ie(i,ran)=(block_data{1,count(i,2)}.data(l,WEIGHT_SIG));
        
        ran=ran+1;
    end
     mweight1=mean(weight_ie(i,280:380));
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
    
       cop(i,ran)=(a(i,ran))/mweight1;
    ran=ran+1;
    end
    
    ran=1;
    
    meanpos=mean(block_data{1,1}.data(:,FOOT_GON_POS_SIG));
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        foot_pos(i,ran)=(block_data{1,count(i,2)}.data(l+221,FOOT_GON_POS_SIG)-meanpos)*IE_foot_gonio*pi/180;
        ran=ran+1;
    end
    foot_pos(i,:)=foot_pos(i,:)-mean(foot_pos(i,1:380));
    
    ran=1; 
    mean_plat_pos=mean(block_data{1,1}.data(:,PLAT_GON_POS_SIG));
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        plat_pos(i,ran)=(block_data{1,count(i,2)}.data(l+221,PLAT_GON_POS_SIG)-mean_plat_pos)*IE_plat_gonio*pi/180;
        ran=ran+1;
    end
    plat_pos(i,:)=plat_pos(i,:)-mean(plat_pos(i,1:380));
    ran=1;
    
    ran=1;
    for l=count(i,1)+TRIAL_WINDOW_PRE_PERT:count(i,1)+TRIAL_WINDOW_POST_PERT
        weight1(i,ran)=(block_data{1,count(i,2)}.data(l,WEIGHT_SIG));
        
        ran=ran+1;
    end
    mweight1=mean(weight1(i,280:380));
    
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
    
 
%     offsetdptorque(i)=mean(dptorque(i,1:200));
    offsetietorque(i)=mean(ietorque(i,340:360));
%     dptorque(i,:)=dptorque(i,:)-offsetdptorque(i);
    ietorque(i,:)=ietorque(i,:)-offsetietorque(i);
    

end

[vel, acc] = get_derivatives(foot_pos, SAMPLE_PERIOD);
[vel2, acc2] = get_derivatives(plat_pos, SAMPLE_PERIOD);

posm=nanmean(foot_pos);
plat_posm=trimmean(plat_pos,30);
am=trimmean(a,30);
% dptorquem=trimmean(dptorque,30);
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
    time_axis_data = TRIAL_WINDOW_PRE_PERT/2:0.5:TRIAL_WINDOW_POST_PERT/2; %Each sample is 0.5 ms
    fig= figure
    ax1=subplot(3,2,1);
    title(strcat('Eversion    ',etitle));

    hold on
    plot(time_axis_data,posm*180/pi,'Color',[0 0 0])
    plot(time_axis_data,plat_posm*180/pi,'Color',[1 0 1])
    %  plot(temp,(dptorquem-dptorquep)*pi/180,'k--')
    axis([-150 250 -inf inf])
    ylabel('angle')
    ax2=subplot(3,2,3);
    hold on
    plot(time_axis_data,velm,'Color',[0 0 0])
    % plot(temp,velp,'Color',[1 0 1])

    axis([-150 250 -inf inf])
    ylabel('velocity')
    clear pk pl
    ax3=subplot(3,2,5);
    hold on
    plot(time_axis_data,weight_iem);
    % plot(time_axis_data,accm,'Color',[0 0 0])
    % [pk,pl]=findpeaks(acc2m);
    %  scatter((pl-400)/2,pk,'k');
    %  plot(time_axis_data,ietorquep,'c');
    %  [pk,pl]=findpeaks(ietorquep);
    %  scatter((pl-400)/2,pk,'c');
    % plot(time_axis_data,ietorquep,'Color',[1 0.5 0.5]);
     axis([-150 250 -inf inf])
    ylabel('weight')
    clear pk pl

    ax4=subplot(3,2,2);
    hold on

    plot(time_axis_data,ietorquem-ietorquep-1*min(k(380:450)),'Color',[0 0 0]);
    hold on
    plot(time_axis_data,ietorquep,'c');
    plot(time_axis_data,ietorquem,'y');
    
%     plot(time_axis_data,accm*0.0530,'r');
    axis([-150 250 -inf inf])
    ylabel('ietorque')
%     imp=regress(ietorquem(350:801)'-ietorquep(350:801)',[ accm(350:801)']);
%  plot(time_axis_data,0.3*accm,'r');
%     imp=regress(dptorquem(380:550)'-dptorquep(380:550)'-0.007*accm(380:550)',[posm(380:550)' velm(380:550)' ]);
%      plot(time_axis_data,posm*imp(1)+velm*imp(2),'r');

    plot(time_axis_data,posm*imp(1)+velm*imp(2)+imp(3)*accm,'r');
    ax5=subplot(3,2,4)
    hold on
    plot(time_axis_data,imp(1)*posm,'r');
    plot(time_axis_data,imp(2)*velm,'g');

    plot(time_axis_data,imp(1)*posm+imp(2)*velm+imp(3)*accm,'m');
    plot(time_axis_data,0.003*accm,'b');

    plot(time_axis_data,ietorquem-ietorquep-1*min(k(380:450)),'k');
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