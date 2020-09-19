close all
clear all
mvc_evaluation;
% Insert subject initial and name. Make sure it matches the format for
% naming
sub_initial='C';
sub_name='CIE';
repetitions=5;
shift=0;
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
t=gonio_values_func('E');
IE_foot_gonio=t(1);
IE_plat_gonio=t(2);
close all;
%%
i=0;
p1=1;
p2=1;
p3=1;
for trials=1:repetitions
    transfer_p1=p1;
    transfer_p2=p2;
    transfer_p3=p3;
    
    if(trials<10)
        h = fopen(strcat(sub_initial,'IE0',num2str(trials),'.dat'));
    else
        h = fopen(strcat(sub_initial,'E',num2str(trials),'.dat'));
    end
    
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    siz=size(Input1.data);
    Img_flag=Input1.data(:,6);
    [x,img_st]=findpeaks(diff(Img_flag));
    %             img_st=round(img_st(1)/20);
    %             if(trials<10)
    %             Img=csvread(strcat(sub_name,'_00',num2str(trials),'.csv'));
    %             else
    %             Img=csvread(strcat(sub_name,'_0',num2str(trials),'.csv'));
    %             end
    %             %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','VNEamplerate',2000);
    %d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 20, 'VNEampleRate', 2000, 'DesignMethod', 'NEindow');
    pert_torque=filtfilt(d1,Input1.data(:,8));
    ta=Input1.data(:,1);
    ta=abs(ta-off_TA)*100/mvc_ta;
    sol=Input1.data(:,2);
    sol=abs(sol-off_SOL)*100/mvc_sol;
    pl=Input1.data(:,3);
    pl=abs(pl-off_PL)*100/mvc_pl;
    gca=Input1.data(:,4);
    gca=abs(gca-off_GCA)*100/mvc_gca;
    w1=filtfilt(d1,Input1.data(:,18));
    flag=Input1.data(:,17);
    foot_pos_data=filtfilt(d1,Input1.data(:,13));
    foot_pos_data=((foot_pos_data-mean(foot_pos_data))*IE_foot_gonio);
    plat_pos_data=filtfilt(d1,Input1.data(:,16));
    plat_pos_data=((plat_pos_data-mean(plat_pos_data))*IE_plat_gonio);
    [test,peaks]=findpeaks(Input1.data(:,17));
 
    i=0;
    for i=1:length(peaks)
        time=[-200:0.5:900];
        if test(i)==1
            
            weight1(p1,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p1_plat_torque(p1,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p1_plat_pos(p1,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p1_foot_pos(p1,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            %                     img1_pos(p1)=getmin(peaks(i),img_st(1),Img);
            p1=p1+1;
            
        end
        if test(i)==2
            weight2(p2,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p2_plat_torque(p2,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p2_plat_pos(p2,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p2_foot_pos(p2,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            %                     img2_pos(p2)=getmin(peaks(i),img_st(1),Img);
            p2=p2+1;
        end
      
        
    end
    i=0;
  
    
    fclose all;
    %             clear live_data Input1 peaks Img;
    
    
    
    
    
    
end

p1_er=0;
p2_er=0;


%%
%
% for i=1:40
%     med=median(img1_pos);
%     [ra idx]=min(abs(img1_pos-med))
%     if(img1_pos(i)< (mean(img1_pos)-std(img1_pos,'omitNaN'))||img1_pos(i)> (mean(img1_pos)+std(img1_pos,'omitNaN')))
%     p1_plat_torque(i,:)=NaN;
%     p1_foot_pos(i,:)=NaN;
%     p1_er=p1_er+1;
%     end
% end
%
% for i=1:40
%     med=median(img2_pos);
%     [ra idx]=min(abs(img2_pos-med))
%     if(img2_pos(i)< (mean(img2_pos)-std(img2_pos))||img2_pos(i)> (mean(img2_pos)+std(img2_pos)))
%     p2_plat_torque(i,:)=NaN;
%     p2_foot_pos(i,:)=NaN;
%     p2_er=p2_er+1;
%     end
% end
% for i=1:40
%     med=median(img3_pos);
%     [ra idx]=min(abs(img3_pos-med))
%     if(img3_pos(i)< (mean(img3_pos)-std(img3_pos))||img3_pos(i)> (mean(img3_pos)+std(img3_pos)))
%     p3_plat_torque(i,:)=NaN;
%     p3_foot_pos(i,:)=NaN;
%     p3_er=p3_er+1;
%     end
% end





%%
weight1m=nanmean(weight1);
weight2m=nanmean(weight2);



p1_plat_torquem=trimmean(p1_plat_torque,30);
p1_plat_posm=trimmean(p1_plat_pos,30);
p1_foot_posm=trimmean(p1_foot_pos,30);
p2_plat_torquem=trimmean(p2_plat_torque,30);
p2_plat_posm=trimmean(p2_plat_pos,30);
p2_foot_posm=trimmean(p2_foot_pos,30);

% 
% for i=1:40
%     
%     ta_emg(i,:)=filtfilt(d3,ta_emg(i,:));
%     pl_emg(i,:)=filtfilt(d3,pl_emg(i,:));
%     sol_emg(i,:)=filtfilt(d3,sol_emg(i,:));
%     gca_emg(i,:)=filtfilt(d3,gca_emg(i,:));
%     ta_emgh(i,:)=filtfilt(d3,ta_emgh(i,:));
%     pl_emgh(i,:)=filtfilt(d3,pl_emgh(i,:));
%     sol_emgh(i,:)=filtfilt(d3,sol_emgh(i,:));
%     gca_emgh(i,:)=filtfilt(d3,gca_emgh(i,:));
%     
% end
% ta_emgm=trimmean(ta_emg,30);
% sol_emgm=trimmean(sol_emg,30);
% pl_emgm=trimmean(pl_emg,30);
% gca_emgm=trimmean(gca_emg,30);
% ta_emghm=trimmean(ta_emg,30);
% sol_emghm=trimmean(sol_emgh,30);
% pl_emghm=trimmean(pl_emgh,30);
% gca_emghm=trimmean(gca_emgh,30);



    rigid_foot=p2_foot_pos;
    comp_foot=p1_foot_pos;
    
    rigid_plat=p2_plat_pos;
    comp_plat=p1_plat_pos;
   
    save('mshaptic.mat','rigid_foot','comp_foot','rigid_plat','comp_plat');
    

 

%%
 variance_analysis