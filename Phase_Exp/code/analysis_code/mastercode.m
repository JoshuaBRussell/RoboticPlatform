close all
clear all
mvc_evaluation;
% Insert subject initial and name. Make sure it matches the format for naming
sub_initial='C';
sub_name='Clayton';
% insert lower limit of inertia of foot in the fit
lim=0.005;
% change these flags to 1 for figures (normal fit and constrained fit)
plot_figs=1;
plot_figs_constrained=1;
% Change to 1 to get ank;le position histogram
plot_hist=1;
%  Change to get torque plot comparison figure
plot_torque=1;
shift=0;
%add the number of loops you want to run bootstrapping here
loops=100;
% Add the trials you want to exclude in here
exclude=[0];%exclude=[1,2];
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
%% Section to calculate goniometer gains
t=gonio_values_func;
DP_foot_gonio=t(1);
DP_plat_gonio=t(2);
%close all    %uncomment to close all figures
%%
i=0;
p1=1;
p2=1;
p3=1;
p0=1;
p4=1;
p5=1;
p6=1;
p7=1;

for trials=1:16
    
  if(ismember(trials,exclude)==0)  
    if(trials<10)
        h = fopen(strcat(sub_initial,'W0',num2str(trials),'.dat'));
    else
        h = fopen(strcat(sub_initial,'W',num2str(trials),'.dat'));
    end
    
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    siz=size(Input1.data);
    Img_flag=Input1.data(:,6);
    [x,img_st]=findpeaks(diff(Img_flag));
    img_st=round(img_st(1)/20);
    if(trials<10)
        Img=csvread(strcat(sub_name,'_00',num2str(trials),'.csv'));
    else
        Img=csvread(strcat(sub_name,'_0',num2str(trials),'.csv'));
    end
    
    %             %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','CWamplerate',2000);
    %d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 20, 'CWampleRate', 2000, 'DesignMethod', 'window');
    pert_torque=filtfilt(d1,Input1.data(:,7));
    ta=Input1.data(:,1);
    ta=abs(ta-off_TA)*100/mvc_ta;
    sol=Input1.data(:,2);
    sol=abs(sol-off_SOL)*100/mvc_sol;
    pl=Input1.data(:,3);
    pl=abs(pl-off_PL)*100/mvc_pl;
    gca=Input1.data(:,4);
    gca=abs(gca-off_GCA)*100/mvc_gca;
    w1=filtfilt(d1,Input1.data(:,18));
    cop=filtfilt(d1,Input1.data(:,19));
    flag=Input1.data(:,17);
    foot_pos_data=filtfilt(d1,Input1.data(:,14));
    foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
    plat_pos_data=filtfilt(d1,Input1.data(:,16));
    plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
    [test,peaks]=findpeaks(Input1.data(:,17));
    for i=1:length(peaks)
        
        time=[-200:0.5:900];
        if test(i)==1
            
            weight1(p1,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p1_plat_torque(p1,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p1_plat_pos(p1,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p1_foot_pos(p1,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img1_pos(p1)=getmin(peaks(i),img_st,Img);
            cop1(p1,:)=cop(peaks(i)-400:peaks(i)+1800);
            p1=p1+1;
            
        end
        if test(i)==2
            weight2(p2,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p2_plat_torque(p2,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p2_plat_pos(p2,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p2_foot_pos(p2,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img2_pos(p2)=getmin(peaks(i),img_st,Img);
            cop2(p2,:)=cop(peaks(i)-400:peaks(i)+1800);
            p2=p2+1;
        end
        if test(i)==3
            weight3(p3,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p3_plat_torque(p3,:)=1*(pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50));
            p3_plat_pos(p3,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p3_foot_pos(p3,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img3_pos(p3)=getmin(peaks(i),img_st,Img);
            cop3(p3,:)=cop(peaks(i)-400:peaks(i)+1800);
            p3=p3+1;
        end
        if test(i)==4
            ta_emg(p0,:)=ta(peaks(i)-400:peaks(i)+1800);
            sol_emg(p0,:)=sol(peaks(i)-400:peaks(i)+1800);
            pl_emg(p0,:)=pl(peaks(i)-400:peaks(i)+1800);
            gca_emg(p0,:)=gca(peaks(i)-400:peaks(i)+1800);
            weight4(p0,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p0_plat_torque(p0,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p0_plat_pos(p0,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p0_foot_pos(p0,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img0_pos(p0)=getmin(peaks(i),img_st,Img);
            cop4(p4,:)=cop(peaks(i)-400:peaks(i)+1800);
            p0=p0+1;
        end
        if test(i)==5
            weight3(p4,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p4_plat_torque(p4,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p4_plat_pos(p4,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p4_foot_pos(p4,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img4_pos(p4)=getmin(peaks(i),img_st,Img);
            p4=p4+1;
        end
        if test(i)==6
            weight3(p5,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p5_plat_torque(p5,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p5_plat_pos(p5,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p5_foot_pos(p5,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img5_pos(p5)=getmin(peaks(i),img_st,Img);
            p5=p5+1;
        end
        if test(i)==7
            weight3(p6,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p6_plat_torque(p6,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p6_plat_pos(p6,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p6_foot_pos(p6,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img6_pos(p6)=getmin(peaks(i),img_st,Img);
            p6=p6+1;
        end
        if test(i)==8
            ta_emgh(p7,:)=ta(peaks(i)-400:peaks(i)+1800);
            sol_emgh(p7,:)=sol(peaks(i)-400:peaks(i)+1800);
            pl_emgh(p7,:)=pl(peaks(i)-400:peaks(i)+1800);
            gca_emgh(p7,:)=gca(peaks(i)-400:peaks(i)+1800);
            weight7(p7,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p7_plat_torque(p7,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p7_plat_pos(p7,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p7_foot_pos(p7,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img7_pos(p7)=getmin(peaks(i),img_st,Img);
            cop7(p7,:)=cop(peaks(i)-400:peaks(i)+1800);
            p7=p7+1;
        end
        
        
    end
    
    fclose all;
    clear live_data Input1 peaks Img;
    
  end
    
end
p0_er=0;
p1_er=0;
p2_er=0;
p3_er=0;
p4_er=0;
p5_er=0;
p6_er=0;
p7_er=0;

for i=1:p1-1
if(img1_pos(i)>500)
img1_pos(i)=NaN;
end
end
for i=1:p2-1
if(img2_pos(i)>500)
img2_pos(i)=NaN;
end
end
for i=1:p3-1
if(img3_pos(i)>500)
img3_pos(i)=NaN;
end
end
for i=1:p0-1
if(img0_pos(i)>500)
img0_pos(i)=NaN;
end
end
for i=1:p4-1
if(img4_pos(i)>500)
img4_pos(i)=NaN;
end
end
for i=1:p5-1
if(img5_pos(i)>500)
img5_pos(i)=NaN;
end
end
for i=1:p6-1
if(img6_pos(i)>500)
img6_pos(i)=NaN;
end
end
for i=1:p7-1
if(img7_pos(i)>500)
img7_pos(i)=NaN;
end
end

if plot_hist==1
im=[img0_pos img1_pos img2_pos img3_pos img4_pos img5_pos img6_pos img7_pos];
end
%%

   

for i=1:p1-1

    med=median(img1_pos);
    [ra idx]=min(abs(img1_pos-med))
    if(img1_pos(i)< (mean(img1_pos)-std(img1_pos,'omitNaN'))||img1_pos(i)> (mean(img1_pos)+std(img1_pos,'omitNaN')))
        p1_plat_torque(i,:)=NaN;
        p1_foot_pos(i,:)=NaN;
        p1_er=p1_er+1;
    end
end

for i=1:p2-1
    med=median(img2_pos);
    [ra idx]=min(abs(img2_pos-med))
    if(img2_pos(i)< (mean(img2_pos)-std(img2_pos))||img2_pos(i)> (mean(img2_pos)+std(img2_pos)))
        p2_plat_torque(i,:)=NaN;
        p2_foot_pos(i,:)=NaN;
        p2_er=p2_er+1;
    end
end
for i=1:p3-1
    med=median(img3_pos);
    [ra idx]=min(abs(img3_pos-med))
    if(img3_pos(i)< (mean(img3_pos)-std(img3_pos))||img3_pos(i)> (mean(img3_pos)+std(img3_pos)))
        p3_plat_torque(i,:)=NaN;
        p3_foot_pos(i,:)=NaN;
        p3_er=p3_er+1;
    end
end
for i=1:p4-1
    med=median(img4_pos);
    [ra idx]=min(abs(img4_pos-med))
    if(img4_pos(i)< (mean(img4_pos)-std(img4_pos))||img4_pos(i)> (mean(img4_pos)+std(img4_pos)))
        p4_plat_torque(i,:)=NaN;
        p4_foot_pos(i,:)=NaN;
        p4_er=p4_er+1;
    end
end
for i=1:p5-1
    med=median(img5_pos);
    [ra idx]=min(abs(img5_pos-med))
    
    if(img5_pos(i)< (mean(img5_pos)-std(img5_pos))||img5_pos(i)> (mean(img5_pos)+std(img5_pos)))
        p5_plat_torque(i,:)=NaN;
        p5_foot_pos(i,:)=NaN;
        p5_er=p5_er+1;
    end
end
for i=1:p6-1
    med=median(img6_pos);
    [ra idx]=min(abs(img6_pos-med))
    if(img6_pos(i)< (mean(img6_pos)-std(img6_pos))||img6_pos(i)> (mean(img6_pos)+std(img6_pos)))
        p6_plat_torque(i,:)=NaN;
        p6_foot_pos(i,:)=NaN;
        p6_er=p6_er+1;
    end
end
for i=1:p7-1
    med=median(img7_pos);
    [ra idx]=min(abs(img7_pos-med))
    if(img7_pos(i)< (mean(img7_pos)-std(img7_pos))||img7_pos(i)> (mean(img7_pos)+std(img7_pos)))
        p7_plat_torque(i,:)=NaN;
        p7_foot_pos(i,:)=NaN;
        p7_er=p7_er+1;
    end
end
for i=1:p0-1
    med=median(img0_pos);
    [ra idx]=min(abs(img0_pos-med))
    if(img0_pos(i)< (mean(img0_pos)-std(img0_pos))||img0_pos(i)> (mean(img0_pos)+std(img0_pos)))
        p0_plat_torque(i,:)=NaN;
        p0_foot_pos(i,:)=NaN;
        p0_er=p0_er+1;
    end
end




%%
weight1m=nanmean(weight1);
weight2m=nanmean(weight2);
weight3m=nanmean(weight3);
weight4m=trimmean(weight4,30);
weight7m=trimmean(weight7,30);
cop4m=trimmean(cop4,30);
cop7m=trimmean(cop7,30);
p0_plat_torquem=trimmean(p0_plat_torque,30);
p0_plat_posm=trimmean(p0_plat_pos,30);
p0_foot_posm=nanmean(p0_foot_pos);
p1_plat_torquem=trimmean(p1_plat_torque,30);
p1_plat_posm=trimmean(p1_plat_pos,30);
p1_foot_posm=trimmean(p1_foot_pos,30);
p2_plat_torquem=trimmean(p2_plat_torque,30);
p2_plat_posm=trimmean(p2_plat_pos,30);
p2_foot_posm=trimmean(p2_foot_pos,30);
p3_plat_torquem=trimmean(p3_plat_torque,30);
p3_plat_posm=trimmean(p3_plat_pos,30);
p3_foot_posm=trimmean(p3_foot_pos,30);
p4_plat_torquem=trimmean(p4_plat_torque,30);
p4_plat_posm=trimmean(p4_plat_pos,30);
p4_foot_posm=trimmean(p4_foot_pos,30);
p5_plat_torquem=trimmean(p5_plat_torque,30);
p5_plat_posm=trimmean(p5_plat_pos,30);
p5_foot_posm=trimmean(p5_foot_pos,30);
p6_plat_torquem=trimmean(p6_plat_torque,30);
p6_plat_posm=trimmean(p6_plat_pos,30);
p6_foot_posm=trimmean(p6_foot_pos,30);
p7_plat_torquem=trimmean(p7_plat_torque,30);
p7_plat_posm=nanmean(p7_plat_pos);
p7_foot_posm=nanmean(p7_foot_pos);

for i=1:p0-1
    
    ta_emg(i,:)=filtfilt(d3,ta_emg(i,:));
    pl_emg(i,:)=filtfilt(d3,pl_emg(i,:));
    sol_emg(i,:)=filtfilt(d3,sol_emg(i,:));
    gca_emg(i,:)=filtfilt(d3,gca_emg(i,:));
    
end

for i=1:p7-1
    
   
    ta_emgh(i,:)=filtfilt(d3,ta_emgh(i,:));
    pl_emgh(i,:)=filtfilt(d3,pl_emgh(i,:));
    sol_emgh(i,:)=filtfilt(d3,sol_emgh(i,:));
    gca_emgh(i,:)=filtfilt(d3,gca_emgh(i,:));
    
    
end

ta_emgm=trimmean(ta_emg,30);
sol_emgm=trimmean(sol_emg,30);
pl_emgm=trimmean(pl_emg,30);
gca_emgm=trimmean(gca_emg,30);
ta_emghm=trimmean(ta_emgh,30);
sol_emghm=trimmean(sol_emgh,30);
pl_emghm=trimmean(pl_emgh,30);
gca_emghm=trimmean(gca_emgh,30);

%%
excluded=[p1,p2,p3,p4,p5,p6];
analysis_value=min(excluded);
for i=1:analysis_value-1
    
    diff_p1_plat_pos(i,:)=p1_plat_pos(i,:);
    diff_p2_plat_pos(i,:)=p2_plat_pos(i,:);
    diff_p3_plat_pos(i,:)=p3_plat_pos(i,:);
    diff_p1_plat_torque(i,:)=p1_plat_torque(i,:)-p0_plat_torquem;
    diff_p2_plat_torque(i,:)=p2_plat_torque(i,:)-p0_plat_torquem;
    diff_p3_plat_torque(i,:)=p3_plat_torque(i,:)-p0_plat_torquem;
    diff_p1_foot_pos(i,:)=p1_foot_pos(i,:)-p0_foot_posm;
    diff_p2_foot_pos(i,:)=p2_foot_pos(i,:)-p0_foot_posm;
    diff_p3_foot_pos(i,:)=p3_foot_pos(i,:)-p0_foot_posm;
    
    diff_p1_foot_pos(i,:)=diff_p1_foot_pos(i,:)-diff_p1_foot_pos(i,700);
    diff_p2_foot_pos(i,:)=diff_p2_foot_pos(i,:)-diff_p2_foot_pos(i,1000);
    diff_p3_foot_pos(i,:)=diff_p3_foot_pos(i,:)-diff_p3_foot_pos(i,1280);
    p1_plat_vel(i,1)=0;
    for l=2:length(p1_plat_posm)
        p1_plat_vel(i,l)=(p1_plat_pos(i,l)-p1_plat_pos(i,l-1))/0.0005;
    end
    p2_plat_vel(i,1)=0;
    for l=2:length(p2_plat_posm)
        p2_plat_vel(i,l)=(p2_plat_pos(i,l)-p2_plat_pos(i,l-1))/0.0005;
    end
    p3_plat_vel(i,1)=0;
    for l=2:length(p2_plat_posm)
        p3_plat_vel(i,l)=(p3_plat_pos(i,l)-p3_plat_pos(i,l-1))/0.0005;
    end
    p1_foot_vel(i,1)=0;
    for l=2:length(p1_foot_posm)
        p1_foot_vel(i,l)=(p1_foot_pos(i,l)-p1_foot_pos(i,l-1))/0.0005;
    end
    p2_foot_vel(i,1)=0;
    for l=2:length(p2_foot_posm)
        p2_foot_vel(i,l)=(p2_foot_pos(i,l)-p2_foot_pos(i,l-1))/0.0005;
    end
    p3_foot_vel(i,1)=0;
    for l=2:length(p3_foot_posm)
        p3_foot_vel(i,l)=(p3_foot_pos(i,l)-p3_foot_pos(i,l-1))/0.0005;
    end
    
    diff_p1_plat_vel(i,1)=0;
    for l=2:length(p1_plat_posm)
        diff_p1_plat_vel(i,l)=(diff_p1_plat_pos(i,l)-diff_p1_plat_pos(i,l-1))/0.0005;
    end
    diff_p2_plat_vel(i,1)=0;
    for l=2:length(p2_plat_posm)
        diff_p2_plat_vel(i,l)=(diff_p2_plat_pos(i,l)-diff_p2_plat_pos(i,l-1))/0.0005;
    end
    diff_p3_plat_vel(i,1)=0;
    for l=2:length(p3_plat_posm)
        diff_p3_plat_vel(i,l)=(diff_p3_plat_pos(i,l)-diff_p3_plat_pos(i,l-1))/0.0005;
    end
    
    %       diff_p1_foot_pos(i,:)=filtfilt(d1,diff_p1_foot_pos(i,:));
    diff_p1_foot_vel(i,1)=0;
    for l=2:length(p1_foot_posm)
        diff_p1_foot_vel(i,l)=(diff_p1_foot_pos(i,l)-diff_p1_foot_pos(i,l-1))/0.0005;
    end
    %    diff_p1_foot_vel(i,:)=filtfilt(d1,diff_p1_foot_vel(i,:));
    diff_p2_foot_vel(i,1)=0;
    for l=2:length(p2_foot_posm)
        diff_p2_foot_vel(i,l)=(diff_p2_foot_pos(i,l)-diff_p2_foot_pos(i,l-1))/0.0005;
    end
    
    diff_p3_foot_vel(i,1)=0;
    for l=2:length(p3_foot_posm)
        diff_p3_foot_vel(i,l)=(diff_p3_foot_pos(i,l)-diff_p3_foot_pos(i,l-1))/0.0005;
    end
    
    p1_plat_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        p1_plat_acc(i,l)=(p1_plat_vel(i,l)-p1_plat_vel(i,l-1))/0.0005;
    end
    p2_plat_acc(i,1)=0;
    for l=2:length(p2_plat_posm)
        p2_plat_acc(i,l)=(p2_plat_vel(i,l)-p2_plat_vel(i,l-1))/0.0005;
    end
    p3_plat_acc(i,1)=0;
    for l=2:length(p3_plat_posm)
        p3_plat_acc(i,l)=(p3_plat_vel(i,l)-p3_plat_vel(i,l-1))/0.0005;
    end
    p1_foot_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        p1_foot_acc(i,l)=(p1_foot_vel(i,l)-p1_foot_vel(i,l-1))/0.0005;
    end
    p2_foot_acc(i,1)=0;
    for l=2:length(p2_plat_posm)
        p2_foot_acc(i,l)=(p2_foot_vel(i,l)-p2_foot_vel(i,l-1))/0.0005;
    end
    p3_foot_acc(i,1)=0;
    for l=2:length(p3_plat_posm)
        p3_foot_acc(i,l)=(p3_foot_vel(i,l)-p3_foot_vel(i,l-1))/0.0005;
    end
    
    diff_p1_plat_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        diff_p1_plat_acc(i,l)=(diff_p1_plat_vel(i,l)-diff_p1_plat_vel(i,l-1))/0.0005;
    end
    diff_p2_plat_acc(i,1)=0;
    for l=2:length(p2_plat_posm)
        diff_p2_plat_acc(i,l)=(diff_p2_plat_vel(i,l)-diff_p2_plat_vel(i,l-1))/0.0005;
    end
    diff_p3_plat_acc(i,1)=0;
    for l=2:length(p3_plat_posm)
        diff_p3_plat_acc(i,l)=(diff_p3_plat_vel(i,l)-diff_p3_plat_vel(i,l-1))/0.0005;
    end
    diff_p1_foot_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        diff_p1_foot_acc(i,l)=(diff_p1_foot_vel(i,l)-diff_p1_foot_vel(i,l-1))/0.0005;
    end
    diff_p2_foot_acc(i,1)=0;
    for l=2:length(p2_plat_posm)
        diff_p2_foot_acc(i,l)=(diff_p2_foot_vel(i,l)-diff_p2_foot_vel(i,l-1))/0.0005;
    end
    diff_p3_foot_acc(i,1)=0;
    for l=2:length(p3_plat_posm)
        diff_p3_foot_acc(i,l)=(diff_p3_foot_vel(i,l)-diff_p3_foot_vel(i,l-1))/0.0005;
    end
    diff_p4_plat_pos(i,:)=p4_plat_pos(i,:);
    diff_p5_plat_pos(i,:)=p5_plat_pos(i,:);
    diff_p6_plat_pos(i,:)=p6_plat_pos(i,:);
    diff_p4_plat_torque(i,:)=p4_plat_torque(i,:)-p7_plat_torquem;
    diff_p5_plat_torque(i,:)=p5_plat_torque(i,:)-p7_plat_torquem;
    diff_p6_plat_torque(i,:)=p6_plat_torque(i,:)-p7_plat_torquem;
    diff_p4_foot_pos(i,:)=p4_foot_pos(i,:)-p7_foot_posm;
    diff_p5_foot_pos(i,:)=p5_foot_pos(i,:)-p7_foot_posm;
    diff_p6_foot_pos(i,:)=p6_foot_pos(i,:)-p7_foot_posm;
    
    diff_p4_foot_pos(i,:)=diff_p4_foot_pos(i,:)-diff_p4_foot_pos(i,720);
    diff_p5_foot_pos(i,:)=diff_p5_foot_pos(i,:)-diff_p5_foot_pos(i,1020);
    diff_p6_foot_pos(i,:)=diff_p6_foot_pos(i,:)-diff_p6_foot_pos(i,1300);
    p4_plat_vel(i,1)=0;
    for l=2:length(p4_plat_posm)
        p4_plat_vel(i,l)=(p4_plat_pos(i,l)-p4_plat_pos(i,l-1))/0.0005;
    end
    p5_plat_vel(i,1)=0;
    for l=2:length(p5_plat_posm)
        p5_plat_vel(i,l)=(p5_plat_pos(i,l)-p5_plat_pos(i,l-1))/0.0005;
    end
    p6_plat_vel(i,1)=0;
    for l=2:length(p5_plat_posm)
        p6_plat_vel(i,l)=(p6_plat_pos(i,l)-p6_plat_pos(i,l-1))/0.0005;
    end
    p4_foot_vel(i,1)=0;
    for l=2:length(p4_foot_posm)
        p4_foot_vel(i,l)=(p4_foot_pos(i,l)-p4_foot_pos(i,l-1))/0.0005;
    end
    p5_foot_vel(i,1)=0;
    for l=2:length(p5_foot_posm)
        p5_foot_vel(i,l)=(p5_foot_pos(i,l)-p5_foot_pos(i,l-1))/0.0005;
    end
    p6_foot_vel(i,1)=0;
    for l=2:length(p6_foot_posm)
        p6_foot_vel(i,l)=(p6_foot_pos(i,l)-p6_foot_pos(i,l-1))/0.0005;
    end
    
    diff_p4_plat_vel(i,1)=0;
    for l=2:length(p4_plat_posm)
        diff_p4_plat_vel(i,l)=(diff_p4_plat_pos(i,l)-diff_p4_plat_pos(i,l-1))/0.0005;
    end
    diff_p5_plat_vel(i,1)=0;
    for l=2:length(p5_plat_posm)
        diff_p5_plat_vel(i,l)=(diff_p5_plat_pos(i,l)-diff_p5_plat_pos(i,l-1))/0.0005;
    end
    diff_p6_plat_vel(i,1)=0;
    for l=2:length(p6_plat_posm)
        diff_p6_plat_vel(i,l)=(diff_p6_plat_pos(i,l)-diff_p6_plat_pos(i,l-1))/0.0005;
    end
    
    
    diff_p4_foot_vel(i,1)=0;
    for l=2:length(p4_foot_posm)
        diff_p4_foot_vel(i,l)=(diff_p4_foot_pos(i,l)-diff_p4_foot_pos(i,l-1))/0.0005;
    end
    diff_p5_foot_vel(i,1)=0;
    for l=2:length(p5_foot_posm)
        diff_p5_foot_vel(i,l)=(diff_p5_foot_pos(i,l)-diff_p5_foot_pos(i,l-1))/0.0005;
    end
    
    diff_p6_foot_vel(i,1)=0;
    for l=2:length(p6_foot_posm)
        diff_p6_foot_vel(i,l)=(diff_p6_foot_pos(i,l)-diff_p6_foot_pos(i,l-1))/0.0005;
    end
    
    p4_plat_acc(i,1)=0;
    for l=2:length(p4_plat_posm)
        p4_plat_acc(i,l)=(p4_plat_vel(i,l)-p4_plat_vel(i,l-1))/0.0005;
    end
    p5_plat_acc(i,1)=0;
    for l=2:length(p5_plat_posm)
        p5_plat_acc(i,l)=(p5_plat_vel(i,l)-p5_plat_vel(i,l-1))/0.0005;
    end
    p6_plat_acc(i,1)=0;
    for l=2:length(p6_plat_posm)
        p6_plat_acc(i,l)=(p6_plat_vel(i,l)-p6_plat_vel(i,l-1))/0.0005;
    end
    p4_foot_acc(i,1)=0;
    for l=2:length(p4_plat_posm)
        p4_foot_acc(i,l)=(p4_foot_vel(i,l)-p4_foot_vel(i,l-1))/0.0005;
    end
    p5_foot_acc(i,1)=0;
    for l=2:length(p5_plat_posm)
        p5_foot_acc(i,l)=(p5_foot_vel(i,l)-p5_foot_vel(i,l-1))/0.0005;
    end
    p6_foot_acc(i,1)=0;
    for l=2:length(p6_plat_posm)
        p6_foot_acc(i,l)=(p6_foot_vel(i,l)-p6_foot_vel(i,l-1))/0.0005;
    end
    
    diff_p4_plat_acc(i,1)=0;
    for l=2:length(p4_plat_posm)
        diff_p4_plat_acc(i,l)=(diff_p4_plat_vel(i,l)-diff_p4_plat_vel(i,l-1))/0.0005;
    end
    diff_p5_plat_acc(i,1)=0;
    for l=2:length(p5_plat_posm)
        diff_p5_plat_acc(i,l)=(diff_p5_plat_vel(i,l)-diff_p5_plat_vel(i,l-1))/0.0005;
    end
    diff_p6_plat_acc(i,1)=0;
    for l=2:length(p6_plat_posm)
        diff_p6_plat_acc(i,l)=(diff_p6_plat_vel(i,l)-diff_p6_plat_vel(i,l-1))/0.0005;
    end
    diff_p4_foot_acc(i,1)=0;
    for l=2:length(p4_plat_posm)
        diff_p4_foot_acc(i,l)=(diff_p4_foot_vel(i,l)-diff_p4_foot_vel(i,l-1))/0.0005;
    end
    diff_p5_foot_acc(i,1)=0;
    for l=2:length(p5_plat_posm)
        diff_p5_foot_acc(i,l)=(diff_p5_foot_vel(i,l)-diff_p5_foot_vel(i,l-1))/0.0005;
    end
    diff_p6_foot_acc(i,1)=0;
    for l=2:length(p6_plat_posm)
        diff_p6_foot_acc(i,l)=(diff_p6_foot_vel(i,l)-diff_p6_foot_vel(i,l-1))/0.0005;
    end
    
end

diff_p1_plat_torquem=trimmean(diff_p1_plat_torque,30);
diff_p1_plat_posm=trimmean(diff_p1_plat_pos,30);
diff_p1_foot_posm=trimmean(diff_p1_foot_pos,30);
diff_p2_plat_torquem=trimmean(diff_p2_plat_torque,30);
diff_p2_plat_posm=trimmean(diff_p2_plat_pos,30);
diff_p2_foot_posm=trimmean(diff_p2_foot_pos,30);
diff_p3_plat_torquem=trimmean(diff_p3_plat_torque,30);
diff_p3_plat_posm=trimmean(diff_p3_plat_pos,30);
diff_p3_foot_posm=trimmean(diff_p3_foot_pos,30);

diff_p1_plat_accm=trimmean(diff_p1_plat_acc,30);
diff_p1_foot_accm=trimmean(diff_p1_foot_acc,30);
diff_p2_plat_accm=trimmean(diff_p2_plat_acc,30);
diff_p2_foot_accm=trimmean(diff_p2_foot_acc,30);
diff_p3_plat_accm=trimmean(diff_p3_plat_acc,30);
diff_p3_foot_accm=trimmean(diff_p3_foot_acc,30);
diff_p4_plat_torquem=trimmean(diff_p4_plat_torque,30);
diff_p4_plat_posm=trimmean(diff_p4_plat_pos,30);
diff_p4_foot_posm=trimmean(diff_p4_foot_pos,30);
diff_p5_plat_torquem=trimmean(diff_p5_plat_torque,30);
diff_p5_plat_posm=trimmean(diff_p5_plat_pos,30);
diff_p5_foot_posm=trimmean(diff_p5_foot_pos,30);
diff_p6_plat_torquem=trimmean(diff_p6_plat_torque,30);
diff_p6_plat_posm=trimmean(diff_p6_plat_pos,30);
diff_p6_foot_posm=trimmean(diff_p6_foot_pos,30);

diff_p4_plat_accm=trimmean(diff_p4_plat_acc,30);
diff_p4_foot_accm=trimmean(diff_p4_foot_acc,30);
diff_p5_plat_accm=trimmean(diff_p5_plat_acc,30);
diff_p5_foot_accm=trimmean(diff_p5_foot_acc,30);
diff_p6_plat_accm=trimmean(diff_p6_plat_acc,30);
diff_p6_foot_accm=trimmean(diff_p6_foot_acc,30);

%%
exc_rigid=[p1,p2,p3,p4,p5,p6];
analysis_value=min(exc_rigid);
for i=1:analysis_value-1
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torque(i,:)-0.1945*diff_p1_plat_accm;%+0.02*diff_p1_foot_accm;%+1*diff_p1_foot_velm;
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torqueimp(i,:)-diff_p1_plat_torqueimp(i,700);
    
    diff_p1_foot_vel(i,:)=diff_p1_foot_vel(i,:);
    diff_p1_foot_acc(i,:)=diff_p1_foot_acc(i,:);
    
    
    
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torque(i,:)-0.1945*diff_p2_plat_accm;%+0.013*diff_p2_foot_accm;
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torqueimp(i,:)-diff_p2_plat_torqueimp(i,1000);
    
    diff_p2_foot_vel(i,:)=diff_p2_foot_vel(i,:);
    diff_p2_foot_acc(i,:)=diff_p2_foot_acc(i,:);
    
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torque(i,:)-0.1945*diff_p3_plat_accm;%-0.09*diff_p3_foot_accm;
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torqueimp(i,:)-diff_p3_plat_torqueimp(i,1280);
    
    diff_p3_foot_vel(i,:)=diff_p3_foot_vel(i,:);
    diff_p3_foot_acc(i,:)=diff_p3_foot_acc(i,:);
end





% for i=1:40
%     if(min(diff_p1_plat_torqueimp(i,700:740))<-2.5)
%         diff_p1_plat_torqueimp(i,:)=NaN;
%         diff_p1_foot_pos(i,:)=NaN;
%         p1_er=p1_er+1;
%     end
% end
% for i=1:40
%     if(min(diff_p2_plat_torqueimp(i,1000:1040))<-1.5)
%         diff_p2_plat_torqueimp(i,:)=NaN;
%         diff_p2_foot_pos(i,:)=NaN;
%         p2_er=p2_er+1;
%     end
% end
% 
% for i=1:40
%     if(min(diff_p3_plat_torqueimp(i,1280:1320))<-0.8)
%         diff_p3_plat_torqueimp(i,:)=NaN;
%         diff_p3_foot_pos(i,:)=NaN;
%         p3_er=p3_er+1;
%     end
% end
% 
% for i=1:40
%     if(min(diff_p1_foot_pos(i,700:740))<-0.005)
%         diff_p1_foot_pos(i,:)=NaN;
%         diff_p1_foot_pos(i,:)=NaN;
%         p1_er=p1_er+1;
%     end
% end
% for i=1:40
%     if(min(diff_p2_foot_pos(i,1000:1040))<-0.003)
%         diff_p2_foot_pos(i,:)=NaN;
%         diff_p2_foot_pos(i,:)=NaN;
%         p2_er=p2_er+1;
%     end
% end
% 
% for i=1:40
%     if(min(diff_p3_foot_pos(i,1280:1320))<-0.003)
%         diff_p3_foot_pos(i,:)=NaN;
%         diff_p3_foot_pos(i,:)=NaN;
%         p3_er=p3_er+1;
%     end
% end


diff_p1_foot_accm=diff_p1_foot_accm;
diff_p2_foot_accm=trimmean(diff_p2_foot_acc,30);
diff_p3_foot_accm=trimmean(diff_p3_foot_acc,30);
diff_p1_plat_torqueimpm=trimmean(diff_p1_plat_torqueimp,30);
diff_p1_plat_torqueimpm=trimmean(diff_p1_plat_torqueimp,30);
diff_p2_plat_torqueimpm=trimmean(diff_p2_plat_torqueimp,30);
diff_p2_plat_torqueimpm=trimmean(diff_p2_plat_torqueimp,30);
diff_p3_plat_torqueimpm=trimmean(diff_p3_plat_torqueimp,30);
diff_p3_plat_torqueimpm=trimmean(diff_p3_plat_torqueimp,30);
diff_p1_foot_posm=trimmean(diff_p1_foot_pos,30);
diff_p2_plat_posm=trimmean(diff_p2_plat_pos,30);
diff_p2_foot_posm=trimmean(diff_p2_foot_pos,30);
diff_p3_plat_posm=trimmean(diff_p3_plat_pos,30);
diff_p3_foot_posm=trimmean(diff_p3_foot_pos,30);
diff_p1_plat_velm=trimmean(diff_p1_plat_vel,30);
diff_p1_foot_velm=trimmean(diff_p1_foot_vel,30);
diff_p2_plat_velm=trimmean(diff_p2_plat_vel,30);
diff_p2_foot_velm=trimmean(diff_p2_foot_vel,30);
diff_p3_plat_velm=trimmean(diff_p3_plat_vel,30);
diff_p3_foot_velm=trimmean(diff_p3_foot_vel,30);
diff_p2_plat_torqueimpm=trimmean(diff_p2_plat_torqueimp,30);
diff_p3_plat_torqueimpm=trimmean(diff_p3_plat_torqueimp,30);

for i=1:analysis_value-1
    C1=[diff_p1_foot_pos(i,700:900)' diff_p1_foot_vel(i,700:900)' diff_p1_foot_acc(i,700:900)'];
     d1=diff_p1_plat_torqueimp(i,700:900)';
     A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B1=[0 ;0 ;1000;1000;-1*lim;0.07];
     if(isnan(d1(1))==0)
     p1imp(i,:)=lsqlin(C1,d1,A1,B1);
     
     else
      p1imp(i,:)=[NaN NaN NaN];   
     end
     vartor3=var(diff_p1_plat_torqueimp(i,700:900));
varimp3=var(diff_p1_plat_torqueimp(i,700:900)-(diff_p1_foot_pos(i,700:900)*p1imp(i,1)+diff_p1_foot_vel(i,700:900)*p1imp(i,2)+diff_p1_foot_acc(i,700:900)*p1imp(i,3)));
p1impg(i)=100*(1-(varimp3/vartor3));
end

p1impm=regress(diff_p1_plat_torqueimpm(700:900)',[diff_p1_foot_posm(700:900)' diff_p1_foot_velm(700:900)' diff_p1_foot_accm(700:900)' ]);
    C=[diff_p1_foot_posm(700:900)' diff_p1_foot_velm(700:900)' diff_p1_foot_accm(700:900)'];
    d=diff_p1_plat_torqueimpm(700:900)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B=[0 ;0 ;1000;1000;-1*lim;0.07];
    p1impm2=lsqlin(C,d,A,B);
    
  
for i=1:analysis_value-1
    
 C1=[diff_p2_foot_pos(i,1000:1200)' diff_p2_foot_vel(i,1000:1200)' diff_p2_foot_acc(i,1000:1200)'];
     d1=diff_p2_plat_torqueimp(i,1000:1200)';
     A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B1=[0 ;0 ;1000;1000;-1*lim;0.07];
     
   if(isnan(d1(1))==0)
     p2imp(i,:)=lsqlin(C1,d1,A1,B1);
     
     else
      p2imp(i,:)=[NaN NaN NaN];   
     end
    
     vartor3=var(diff_p2_plat_torqueimp(i,1000:1200));
varimp3=var(diff_p2_plat_torqueimp(i,1000:1200)-(diff_p2_foot_pos(i,1000:1200)*p2imp(i,1)+diff_p2_foot_vel(i,1000:1200)*p2imp(i,2)+diff_p2_foot_acc(i,1000:1200)*p2imp(i,3)));
p2impg(i)=100*(1-(varimp3/vartor3));
    
end
 p2impm=regress(diff_p2_plat_torqueimpm(1000:1200)',[diff_p2_foot_posm(1000:1200)' diff_p2_foot_velm(1000:1200)' diff_p2_foot_accm(1000:1200)']);
    C=[diff_p2_foot_posm(1000:1200)' diff_p2_foot_velm(1000:1200)' diff_p2_foot_accm(1000:1200)'];
    d=diff_p2_plat_torqueimpm(1000:1200)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B=[0 ;0 ;1000;1000;-1*lim;0.07];
    p2impm2=lsqlin(C,d,A,B);

for i=1:analysis_value-1
    
   C1=[diff_p3_foot_pos(i,1280:1480)' diff_p3_foot_vel(i,1280:1480)' diff_p3_foot_acc(i,1280:1480)'];
     d1=diff_p3_plat_torqueimp(i,1280:1480)';
     A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B1=[0 ;0 ;1000;1000;-1*lim;0.07];
     
     if(isnan(d1(1))==0)
     p3imp(i,:)=lsqlin(C1,d1,A1,B1);
     
     else
      p3imp(i,:)=[NaN NaN NaN];   
     end
    
     vartor3=var(diff_p3_plat_torqueimp(i,1280:1480));
varimp3=var(diff_p3_plat_torqueimp(i,1280:1480)-(diff_p3_foot_pos(i,1280:1480)*p3imp(i,1)+diff_p3_foot_vel(i,1280:1480)*p3imp(i,2)+diff_p3_foot_acc(i,1280:1480)*p3imp(i,3)));
p3impg(i)=100*(1-(varimp3/vartor3));
end

p3impm=regress(diff_p3_plat_torqueimpm(1280:1480)',[diff_p3_foot_posm(1280:1480)' diff_p3_foot_velm(1280:1480)' diff_p3_foot_accm(1280:1480)']);
    
    C=[diff_p3_foot_posm(1280:1480)' diff_p3_foot_velm(1280:1480)' diff_p3_foot_accm(1280:1480)'];
    d=diff_p3_plat_torqueimpm(1280:1480)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B=[0 ;0 ;1000;1000;-1*lim;0.07];
    p3impm2=lsqlin(C,d,A,B);

exc_haptic=[p4,p5,p6];
% analysis_value=min(exc_haptic);

for i=1:analysis_value-1
    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torque(i,:)-0.1945*diff_p4_plat_accm;%+0.02*diff_p4_foot_accm;%+1*diff_p4_foot_velm;
    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torqueimp(i,:)-diff_p4_plat_torqueimp(i,720);
    
    diff_p4_foot_vel(i,:)=diff_p4_foot_vel(i,:);
    diff_p4_foot_acc(i,:)=diff_p4_foot_acc(i,:);
    
    
    
    diff_p5_plat_torqueimp(i,:)=diff_p5_plat_torque(i,:)-0.1945*diff_p5_plat_accm;%+0.013*diff_p5_foot_accm;
    diff_p5_plat_torqueimp(i,:)=diff_p5_plat_torqueimp(i,:)-diff_p5_plat_torqueimp(i,1020);
    
    diff_p5_foot_vel(i,:)=diff_p5_foot_vel(i,:);
    diff_p5_foot_acc(i,:)=diff_p5_foot_acc(i,:);
    
    diff_p6_plat_torqueimp(i,:)=diff_p6_plat_torque(i,:)-0.1945*diff_p6_plat_accm;%-0.09*diff_p6_foot_accm;
    diff_p6_plat_torqueimp(i,:)=diff_p6_plat_torqueimp(i,:)-diff_p6_plat_torqueimp(i,1300);
    
    diff_p6_foot_vel(i,:)=diff_p6_foot_vel(i,:);
    diff_p6_foot_acc(i,:)=diff_p6_foot_acc(i,:);
end
diff_p5_plat_torqueimp(34,:)=diff_p5_plat_torqueimp(33,:);

% for i=1:40
%     if(min(diff_p4_plat_torqueimp(i,720:760))<-1.2)
%         diff_p4_plat_torqueimp(i,:)=NaN;
%         diff_p4_foot_pos(i,:)=NaN;
%         p4_er=p4_er+1;
%     end
% end
% for i=1:40
%     if(min(diff_p5_plat_torqueimp(i,1020:1060))<-1.5)
%         diff_p5_plat_torqueimp(i,:)=NaN;
%         diff_p5_foot_pos(i,:)=NaN;
%         p5_er=p5_er+1;
%     end
% end
% 
% for i=1:40
%     if(min(diff_p4_foot_pos(i,720:760))<-0.005)
%         diff_p4_foot_pos(i,:)=NaN;
%         diff_p4_foot_pos(i,:)=NaN;
%         p4_er=p4_er+1;
%     end
% end
% for i=1:40
%     if(min(diff_p5_foot_pos(i,1020:1060))<-0.003)
%         diff_p5_foot_pos(i,:)=NaN;
%         diff_p5_foot_pos(i,:)=NaN;
%         p5_er=p5_er+1;
%     end
% end
% 
% for i=1:40
%     if(min(diff_p6_foot_pos(i,1300:1340))<-0.003)
%         diff_p6_foot_pos(i,:)=NaN;
%         diff_p6_foot_pos(i,:)=NaN;
%         p6_er=p6_er+1;
%     end
% end

diff_p4_foot_accm=diff_p4_foot_accm;
diff_p5_foot_accm=trimmean(diff_p5_foot_acc,30);
diff_p6_foot_accm=trimmean(diff_p6_foot_acc,30);
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);
diff_p5_plat_torqueimpm=trimmean(diff_p5_plat_torqueimp,30);
diff_p5_plat_torqueimpm=trimmean(diff_p5_plat_torqueimp,30);
diff_p6_plat_torqueimpm=trimmean(diff_p6_plat_torqueimp,30);
diff_p6_plat_torqueimpm=trimmean(diff_p6_plat_torqueimp,30);
diff_p4_foot_posm=trimmean(diff_p4_foot_pos,30);
diff_p5_plat_posm=trimmean(diff_p5_plat_pos,30);
diff_p5_foot_posm=trimmean(diff_p5_foot_pos,30);
diff_p6_plat_posm=trimmean(diff_p6_plat_pos,30);
diff_p6_foot_posm=trimmean(diff_p6_foot_pos,30);
diff_p4_plat_velm=trimmean(diff_p4_plat_vel,30);
diff_p4_foot_velm=trimmean(diff_p4_foot_vel,30);
diff_p5_plat_velm=trimmean(diff_p5_plat_vel,30);
diff_p5_foot_velm=trimmean(diff_p5_foot_vel,30);
diff_p6_plat_velm=trimmean(diff_p6_plat_vel,30);
diff_p6_foot_velm=trimmean(diff_p6_foot_vel,30);
diff_p5_plat_torqueimpm=trimmean(diff_p5_plat_torqueimp,30);
diff_p6_plat_torqueimpm=trimmean(diff_p6_plat_torqueimp,30);




for i=1:analysis_value-1
    
%     p4imp(i,:)=regress(diff_p4_plat_torqueimp(i,720:920)',[diff_p4_foot_pos(i,720:920)' diff_p4_foot_vel(i,720:920)' diff_p4_foot_acc(i,720:920)']);
   C1=[diff_p4_foot_pos(i,720:920)' diff_p4_foot_vel(i,720:920)' diff_p4_foot_acc(i,720:920)'];
     d1=diff_p4_plat_torqueimp(i,720:920)';
     A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B1=[0 ;0 ;1000;1000;-1*lim;0.07];
     
      if(isnan(d1(1))==0)
     p4imp(i,:)=lsqlin(C1,d1,A1,B1);
     
     else
      p4imp(i,:)=[NaN NaN NaN];   
     end
    
     vartor3=var(diff_p4_plat_torqueimp(i,720:920));
varimp4=var(diff_p4_plat_torqueimp(i,720:920)-(diff_p4_foot_pos(i,720:920)*p4imp(i,1)+diff_p4_foot_vel(i,720:920)*p4imp(i,2)+diff_p4_foot_acc(i,720:920)*p4imp(i,3)));
p4impg(i)=100*(1-(varimp4/vartor3));
    
end
 p4impm=regress(diff_p4_plat_torqueimpm(720:920)',[diff_p4_foot_posm(720:920)' diff_p4_foot_velm(720:920)' diff_p4_foot_accm(720:920)']);
    C=[diff_p4_foot_posm(720:920)' diff_p4_foot_velm(720:920)' diff_p4_foot_accm(720:920)'];
    d=diff_p4_plat_torqueimpm(720:920)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B=[0 ;0 ;1000;1000;-1*lim;0.07];
    p4impm2=lsqlin(C,d,A,B);


for i=1:analysis_value-1
    
   C1=[diff_p5_foot_pos(i,1020:1220)' diff_p5_foot_vel(i,1020:1220)' diff_p5_foot_acc(i,1020:1220)'];
     d1=diff_p5_plat_torqueimp(i,1020:1220)';
     A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B1=[0 ;0 ;1000;1000;-1*lim;0.07];
     
     if(isnan(d1(1))==0)
     p5imp(i,:)=lsqlin(C1,d1,A1,B1);
     
     else
      p5imp(i,:)=[NaN NaN NaN];   
     end
    
     vartor3=var(diff_p5_plat_torqueimp(i,1020:1220));
varimp5=var(diff_p5_plat_torqueimp(i,1020:1220)-(diff_p5_foot_pos(i,1020:1220)*p5imp(i,1)+diff_p5_foot_vel(i,1020:1220)*p5imp(i,2)+diff_p5_foot_acc(i,1020:1220)*p5imp(i,3)));
p5impg(i)=100*(1-(varimp5/vartor3));
    
    
end
    p5impm=regress(diff_p5_plat_torqueimpm(1020:1220)',[diff_p5_foot_posm(1020:1220)' diff_p5_foot_velm(1020:1220)' diff_p5_foot_accm(1020:1220)']);
    C=[diff_p5_foot_posm(1020:1220)' diff_p5_foot_velm(1020:1220)' diff_p5_foot_accm(1020:1220)'];
    d=diff_p5_plat_torqueimpm(1020:1220)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B=[0 ;0 ;1000;1000;-1*lim;0.07];
    p5impm2=lsqlin(C,d,A,B);


for i=1:analysis_value-1
C1=[diff_p6_foot_pos(i,1300:1500)' diff_p6_foot_vel(i,1300:1500)' diff_p6_foot_acc(i,1300:1500)'];
     d1=diff_p6_plat_torqueimp(i,1300:1500)';
     A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
      B1=[0 ;0 ;1000;1000;-1*lim;0.07];
     
     if(isnan(d1(1))==0)
     p6imp(i,:)=lsqlin(C1,d1,A1,B1);
     
     else
      p6imp(i,:)=[NaN NaN NaN];   
     end
    
     vartor3=var(diff_p6_plat_torqueimp(i,1300:1500));
varimp6=var(diff_p6_plat_torqueimp(i,1300:1500)-(diff_p6_foot_pos(i,1300:1500)*p6imp(i,1)+diff_p6_foot_vel(i,1300:1500)*p6imp(i,2)+diff_p6_foot_acc(i,1300:1500)*p6imp(i,3)));
p6impg(i)=100*(1-(varimp6/vartor3));
    
end
p6impm=regress(diff_p6_plat_torqueimpm(1300:1500)',[diff_p6_foot_posm(1300:1500)' diff_p6_foot_velm(1300:1500)' diff_p6_foot_accm(1300:1500)']);
    
    C=[diff_p6_foot_posm(1300:1500)' diff_p6_foot_velm(1300:1500)' diff_p6_foot_accm(1300:1500)'];
    d=diff_p6_plat_torqueimpm(1300:1500)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    p6impm2=lsqlin(C,d,A,B);


%%
impedance=[p1impm';p2impm';p3impm';p4impm';p5impm';p6impm'];
impedance_constrained=[p1impm2';p2impm2';p3impm2';p4impm2';p5impm2';p6impm2'];
goodness;
goodness_constrained;
Final_value=[impedance, goodnessn];
 save('impedance.mat','Final_value');
Final_value_constrained=[impedance_constrained, goodnessc];
save('impedance_constrained.mat','Final_value_constrained');
if(plot_torque==1)
torque_comparison_plot;
end

if(plot_hist==1)
figure
hist(im,40);
saveas(gcf,'histogram.jpg');
end

if(plot_figs==1)
   fit_figures; 
end
if(plot_figs_constrained==1)
   fit_figures_constrained; 
end

plottingemgweight_old;
bootstrapping;
boot_goodness;
write_data_summary;
