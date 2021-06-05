close all
clear all



% Insert subject initial and name. Make sure it matches the format for naming
sub_initial='V';
sub_name='Vu';

DATA_FOLDER_REL_LOC = "./";

NUM_OF_BLOCKS = 6;
% insert lower limit of inertia of foot in the fit
% u_lim is the upper limit of the inertia and lim
% is the lower limit
lim= 0.007;
u_lim=0.02;

% change these flags to 1 for figures (normal fit and constrained fit)
plot_figs=1;
plot_figs_constrained=1;
% Change to 1 to get ank;le position histogram
plot_hist=1;
%  Change to get torque plot comparison figure
plot_torque=0;
shift=0;

%The PreTrial window starts 400 samples minus the peak start. Then the perturbation happens,
%900(?) samples after word -> 400+900 = 1300
%THIS WILL NEED TO CHANGE FOR EACH INDIVIDUAL SUBJECT
REGRESSION_WINDOW_MIN_INDEX = 1300;
REGRESSION_WINDOW_MAX_INDEX = 1500;



d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
%% Section to calculate goniometer gains
% t=gonio_values_func('D');
% DP_foot_gonio=t(1);
% DP_plat_gonio=t(2);
perturb_type = 'D';
[DP_foot_gonio, DP_plat_gonio] = get_gonio_sf(DATA_FOLDER_REL_LOC, perturb_type); %Finding gains of goniometer

mvc_evaluation;
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

for trials=1:NUM_OF_BLOCKS
    
    
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
    flag=Input1.data(:,17);
    foot_pos_data=filtfilt(d1,Input1.data(:,13));
    foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
    plat_pos_data=filtfilt(d1,Input1.data(:,14));
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
            p1=p1+1;
            
        end
        if test(i)==2
            ta_emg(p0,:)=ta(peaks(i)-400:peaks(i)+1800);
            sol_emg(p0,:)=sol(peaks(i)-400:peaks(i)+1800);
            pl_emg(p0,:)=pl(peaks(i)-400:peaks(i)+1800);
            gca_emg(p0,:)=gca(peaks(i)-400:peaks(i)+1800);
            weight4(p0,:)=w1(peaks(i)-400:peaks(i)+1800)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+1800)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+1800)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+1800)-w4(peaks(i)-20);
            p0_plat_torque(p0,:)=pert_torque(peaks(i)-400:peaks(i)+1800)-pert_torque(peaks(i)+50);
            p0_plat_pos(p0,:)=plat_pos_data(peaks(i)-179:peaks(i)+2021);
            p0_foot_pos(p0,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2021+shift);
            img0_pos(p0)=getmin(peaks(i),img_st,Img);
            p0=p0+1;
        end
        
        
        
    end
    
    fclose all;
    clear live_data Input1 peaks Img;
    
    
    
end
p0_er=0;
p1_er=0;

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

weight4m=trimmean(weight4,30);
p0_plat_torquem=trimmean(p0_plat_torque,30);
p0_plat_posm=trimmean(p0_plat_pos,30);
p0_foot_posm=nanmean(p0_foot_pos);
p1_plat_torquem=trimmean(p1_plat_torque,30);
p1_plat_posm=trimmean(p1_plat_pos,30);
p1_foot_posm=trimmean(p1_foot_pos,30);

for i=1:p0-1
    
    ta_emg(i,:)=filtfilt(d3,ta_emg(i,:));
    pl_emg(i,:)=filtfilt(d3,pl_emg(i,:));
    sol_emg(i,:)=filtfilt(d3,sol_emg(i,:));
    gca_emg(i,:)=filtfilt(d3,gca_emg(i,:));
    
end



ta_emgm=trimmean(ta_emg,30);
sol_emgm=trimmean(sol_emg,30);
pl_emgm=trimmean(pl_emg,30);
gca_emgm=trimmean(gca_emg,30);


%%
excluded=[p1];
analysis_value=min(excluded);
for i=1:analysis_value-1
    
    diff_p1_plat_pos(i,:)=p1_plat_pos(i,:);
   
    diff_p1_plat_torque(i,:)=p1_plat_torque(i,:)-p0_plat_torquem;
    
    diff_p1_foot_pos(i,:)=p1_foot_pos(i,:)-p0_foot_posm;
  
    
    diff_p1_foot_pos(i,:)=diff_p1_foot_pos(i,:)-diff_p1_foot_pos(i,REGRESSION_WINDOW_MIN_INDEX);
  
    p1_plat_vel(i,1)=0;
    for l=2:length(p1_plat_posm)
        p1_plat_vel(i,l)=(p1_plat_pos(i,l)-p1_plat_pos(i,l-1))/0.0005;
    end
   
    p1_foot_vel(i,1)=0;
    for l=2:length(p1_foot_posm)
        p1_foot_vel(i,l)=(p1_foot_pos(i,l)-p1_foot_pos(i,l-1))/0.0005;
    end
   
    
    diff_p1_plat_vel(i,1)=0;
    for l=2:length(p1_plat_posm)
        diff_p1_plat_vel(i,l)=(diff_p1_plat_pos(i,l)-diff_p1_plat_pos(i,l-1))/0.0005;
    end
   
    %       diff_p1_foot_pos(i,:)=filtfilt(d1,diff_p1_foot_pos(i,:));
    diff_p1_foot_vel(i,1)=0;
    for l=2:length(p1_foot_posm)
        diff_p1_foot_vel(i,l)=(diff_p1_foot_pos(i,l)-diff_p1_foot_pos(i,l-1))/0.0005;
    end
   
  
    
    p1_plat_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        p1_plat_acc(i,l)=(p1_plat_vel(i,l)-p1_plat_vel(i,l-1))/0.0005;
    end
    
    p1_foot_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        p1_foot_acc(i,l)=(p1_foot_vel(i,l)-p1_foot_vel(i,l-1))/0.0005;
    end
    
    
    diff_p1_plat_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        diff_p1_plat_acc(i,l)=(diff_p1_plat_vel(i,l)-diff_p1_plat_vel(i,l-1))/0.0005;
    end
    
    diff_p1_foot_acc(i,1)=0;
    for l=2:length(p1_plat_posm)
        diff_p1_foot_acc(i,l)=(diff_p1_foot_vel(i,l)-diff_p1_foot_vel(i,l-1))/0.0005;
    end
   
    
end

diff_p1_plat_torquem=trimmean(diff_p1_plat_torque,30);
diff_p1_plat_posm=trimmean(diff_p1_plat_pos,30);
diff_p1_foot_posm=trimmean(diff_p1_foot_pos,30);

diff_p1_plat_accm=trimmean(diff_p1_plat_acc,30);
diff_p1_foot_accm=trimmean(diff_p1_foot_acc,30);


%%
exc_rigid=[p1];
analysis_value=min(exc_rigid);
for i=1:analysis_value-1
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torque(i,:)-0.1945*diff_p1_plat_accm;%+0.02*diff_p1_foot_accm;%+1*diff_p1_foot_velm;
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torqueimp(i,:)-diff_p1_plat_torqueimp(i,REGRESSION_WINDOW_MIN_INDEX);
    
    diff_p1_foot_vel(i,:)=diff_p1_foot_vel(i,:);
    diff_p1_foot_acc(i,:)=diff_p1_foot_acc(i,:);
end




diff_p1_foot_accm=diff_p1_foot_accm;

diff_p1_plat_torqueimpm=trimmean(diff_p1_plat_torqueimp,30);
diff_p1_plat_torqueimpm=trimmean(diff_p1_plat_torqueimp,30);

diff_p1_foot_posm=trimmean(diff_p1_foot_pos,30);

diff_p1_plat_velm=trimmean(diff_p1_plat_vel,30);
diff_p1_foot_velm=trimmean(diff_p1_foot_vel,30);

for i=1:analysis_value-1
    
    p1imp(i,:)=regress(diff_p1_plat_torqueimp(i,REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)',[diff_p1_foot_pos(i,REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' diff_p1_foot_vel(i,REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' diff_p1_foot_acc(i,REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)']);
    p1impm=regress(diff_p1_plat_torqueimpm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)',[diff_p1_foot_posm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' diff_p1_foot_velm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' diff_p1_foot_accm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' ]);
    C=[diff_p1_foot_posm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' diff_p1_foot_velm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)' diff_p1_foot_accm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)'];
    d=diff_p1_plat_torqueimpm(REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
%     A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
%     B=[0 ;0 ;1000;1000;-1*lim;0.07];
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;u_lim];
    p1impm2=lsqlin(C,d,A,B);
    
    p1imp3(1)=trimmean(p1imp(:,1),30);
    p1imp3(2)=trimmean(p1imp(:,2),30);
    p1imp3(3)=trimmean(p1imp(:,3),30);
end




%%
impedance=[p1impm'];
impedance_constrained=[p1impm2'];
goodness;
goodness_constrained;
Final_value=[impedance, goodnessn];
Final_value_constrained=[impedance_constrained, goodnessc];


if(plot_figs==1)
   fit_figures; 
end
if(plot_figs_constrained==1)
   fit_figures_constrained; 
end

