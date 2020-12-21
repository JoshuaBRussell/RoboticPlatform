close all
clear all

%% ----Constants---- %%
%Insert subject initial and name.
%Make sure it matches the format for naming
SUB_INITIAL='M';
SUB_NAME='M';

NUM_PERT=40;                  %Number of perturbation you actually ran

INERTIAL_LOWER_LIM  = 0.01;   %Inertial Parametrs Lower/Upper Limits
INERTIAL_UPPER_LIMIT = 0.2;

% Flag == 1 -> Shows Figures
PLOT_FIGS_FLAG=1;
% plot_figs_constrained=0;
PLOT_HISTOGRAM_FLAG=1;        %Ankle position histogram
PLOT_TORQUE_FLAG=1;           %Torque plot comparison figure

SHIFT_VAL=0;

BOOTSTRAPPING_LOOPS=100;      %Number of bootstrapping loops

TRIAL_EXCLUDE_VEC=[];         %Trial number to be excluded from analysis

VOLTS_TO_NEWTONS_SCALAR = 53.4;

%DATA File Input Channels
TA_EMG_CHANN =  1;
SOL_EMG_CHANN = 2;
PL_EMG_CHANN =  3;
GCA_EMG_CHANN = 4;

FORCE_SENSOR_1_CHANN = 9;
FORCE_SENSOR_2_CHANN = 10;
FORCE_SENSOR_3_CHANN = 11;
FORCE_SENSOR_4_CHANN = 12;
FORCE_SENSOR_5_CHANN = 20;
FORCE_SENSOR_6_CHANN = 21;

ANN_RIGID_PHASE_CHANN = 15;
ANN_COMP_PHASE_CHANN = 16;

WEIGHT_CHANN = 18;
COP_CHANN =    19;

IMG_FLAG_CHANN = 6;
PERT_START_FLAG_CHANN = 22;
DP_TORQUE_CHANN = 7;
IE_TORQUE_CHANN = 8;

FOOT_POS_CHANN = 13;
PLAT_POS_CHANN = 14;

%Window/Timing Constants
PRE_PURT_TRIG_SAMPLES = 400;
POST_PURT_TRIG_SAMPLES = 2000;

FORCE_BASELINE_SHIFT = 360;
WEIGHT_BASELINE_SHIFT = 360;
PERT_TORQUE_SHIFT = 50;

GONIOMETER_SAMPLE_DELAY = 221; 

PRE_PERT_SAMPLES = 100;
POST_PERT_SAMPLES = 300;

REGRESSION_START_INDEX = 100;
REGRESSION_END_INDEX = 300;

PERT_POINT1_PERC = 0.18;
PERT_POINT2_PERC = 0.31;
PERT_POINT3_PERC = 0.44;
PERT_POINT4_PERC = 0.57;

FORCE_PLATE_DIST_X = 0.0;
FORCE_PLATE_DIST_Y = 0.0;
FORCE_PLATE_DIST_A = 0.0;
FORCE_PLATE_DIST_B = 0.0;

SAMPLE_TIME = 1/2000; 


PLAT_DP_INERTIA = 0.1945;
PLAT_IE_INERTIA = 0.0;


IMG_POS_CUTOFF_LIMIT = 500;



%% Setup, Filtering and Data Visualization settings

% Getting the MVC values
mvc_evaluation;
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',...
    5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',15,...
    'DesignMethod','butter','Samplerate',2000);
%% Section to calculate goniometer gains
t=gonio_values_func;
DP_foot_gonio=t(1);
DP_plat_gonio=t(2);

%close all    %uncomment to close all figures
%% Acquisition of trial data from different files
i=0;
p1=1;
p2=1;
p3=1;
p0=1;
p4=1;
p5=1;
p6=1;
p7=1;
% Make sure you only add the total number of
% trials you actually run
for trials=1:10
    
    if(ismember(trials,TRIAL_EXCLUDE_VEC)==0)
        if(trials<10)
            h = fopen(strcat(SUB_INITIAL,'W0',num2str(trials),'.dat'));
        else
            h = fopen(strcat(SUB_INITIAL,'W',num2str(trials),'.dat'));
        end
        
        live_data=fread(h);
        Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
        siz=size(Input1.data);
        Img_flag=Input1.data(:,IMG_FLAG_CHANN);
        [x,img_st]=findpeaks(diff(Img_flag));
        img_st=round(img_st(1)/20);
        if(trials<10)
            Img=csvread(strcat(SUB_NAME,'_00',num2str(trials),'.csv'));
        else
            Img=csvread(strcat(SUB_NAME,'_0',num2str(trials),'.csv'));
        end
        % obtaining data from channels
        pert_torque=filtfilt(d1,Input1.data(:,DP_TORQUE_CHANN));
        f1=Input1.data(:,FORCE_SENSOR_1_CHANN)*VOLTS_TO_NEWTONS_SCALAR;
        f1=filtfilt(d1,f1);
        f2=Input1.data(:,FORCE_SENSOR_2_CHANN)*VOLTS_TO_NEWTONS_SCALAR;
        f2=filtfilt(d1,f2);
        f3=Input1.data(:,FORCE_SENSOR_3_CHANN)*VOLTS_TO_NEWTONS_SCALAR;
        f3=filtfilt(d1,f3);
        f4=Input1.data(:,FORCE_SENSOR_4_CHANN)*VOLTS_TO_NEWTONS_SCALAR;
        f4=filtfilt(d1,f4);
        f5=Input1.data(:,FORCE_SENSOR_5_CHANN)*VOLTS_TO_NEWTONS_SCALAR/2;
        f5=filtfilt(d1,f5);
        f6=Input1.data(:,FORCE_SENSOR_6_CHANN)*VOLTS_TO_NEWTONS_SCALAR/2;
        f6=filtfilt(d1,f6);
        ta=Input1.data(:,TA_EMG_CHANN);
        ta=abs(ta-off_TA)*100/mvc_ta;
        sol=Input1.data(:,SOL_EMG_CHANN);
        sol=abs(sol-off_SOL)*100/mvc_sol;
        pl=Input1.data(:,PL_EMG_CHANN);
        pl=abs(pl-off_PL)*100/mvc_pl;
        gca=Input1.data(:,GCA_EMG_CHANN);
        gca=abs(gca-off_GCA)*100/mvc_gca;
        w1=filtfilt(d1,Input1.data(:,WEIGHT_CHANN));
        cop=filtfilt(d1,Input1.data(:,COP_CHANN));
        flag=Input1.data(:, 17);
        rigid_phase_tot=Input1.data(:,ANN_RIGID_PHASE_CHANN);

        perturb_start=Input1.data(:,PERT_START_FLAG_CHANN);
        foot_pos_data=filtfilt(d1,Input1.data(:,FOOT_POS_CHANN));
        foot_pos_data=((foot_pos_data-mean(foot_pos_data))...
            *DP_foot_gonio*pi/180);
        plat_pos_data=filtfilt(d1,Input1.data(:,PLAT_POS_CHANN));
        plat_pos_data=((plat_pos_data-mean(plat_pos_data))...
            *DP_plat_gonio*pi/180);
        % 17 records an impulse everytime a
        % perturbation occurs with diff amplitudes
        [test,peaks]=findpeaks(Input1.data(:,17));
        for i=1:length(peaks)
     % Each test value corresponds to a different pert       
            time=[-200:0.5:1000];
            if test(i)==1
                force1_1(p1,:)=f1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force1_2(p1,:)=f2(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force1_3(p1,:)=f3(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force1_4(p1,:)=f4(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force1_5(p1,:)=f5(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force1_6(p1,:)=f6(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                weight1(p1,:)=w1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-w1(peaks(i)-WEIGHT_BASELINE_SHIFT);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p1_plat_torque(p1,:)=pert_torque(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-pert_torque(peaks(i)+50);
                p1_plat_pos(p1,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p1_foot_pos(p1,:)=foot_pos_data(peaks(i)-179+SHIFT_VAL:peaks(i)+2221+SHIFT_VAL);
                p1_phase(p1,:)=rigid_phase_tot(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                p1_pert(p1,:)=perturb_start(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                img1_pos(p1)=getmin(peaks(i),img_st,Img);
                cop1(p1,:)=cop(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                [a,b]=findpeaks(diff(p1_pert(p1,:)));
                p1_peakst(p1)=b;
                p1=p1+1;
                
            end
            if test(i)==2
                force2_1(p2,:)=f1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force2_2(p2,:)=f2(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force2_3(p2,:)=f3(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force2_4(p2,:)=f4(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force2_5(p2,:)=f5(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force2_6(p2,:)=f6(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                
                weight2(p2,:)=w1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-w1(peaks(i)-WEIGHT_BASELINE_SHIFT);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p2_plat_torque(p2,:)=pert_torque(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-pert_torque(peaks(i)+50);
                p2_plat_pos(p2,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p2_foot_pos(p2,:)=foot_pos_data(peaks(i)-179+SHIFT_VAL:peaks(i)+2221+SHIFT_VAL);
                p2_phase(p2,:)=rigid_phase_tot(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                p2_pert(p2,:)=perturb_start(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                img2_pos(p2)=getmin(peaks(i),img_st,Img);
                cop2(p2,:)=cop(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                [a,b]=findpeaks(diff(p2_pert(p2,:)));

                p2_peakst(p2)=b;
                p2=p2+1;
            end
            if test(i)==3
                force3_1(p3,:)=f1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force3_2(p3,:)=f2(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force3_3(p3,:)=f3(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force3_4(p3,:)=f4(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force3_5(p3,:)=f5(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force3_6(p3,:)=f6(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                
                weight3(p3,:)=w1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-w1(peaks(i)-WEIGHT_BASELINE_SHIFT);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p3_plat_torque(p3,:)=1*(pert_torque(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-pert_torque(peaks(i)+50));
                p3_plat_pos(p3,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p3_foot_pos(p3,:)=foot_pos_data(peaks(i)-179+SHIFT_VAL:peaks(i)+2221+SHIFT_VAL);
                p3_phase(p3,:)=rigid_phase_tot(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                p3_pert(p3,:)=perturb_start(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                img3_pos(p3)=getmin(peaks(i),img_st,Img);
                cop3(p3,:)=cop(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                [a,b]=findpeaks(diff(p3_pert(p3,:)));
                p3_peakst(p3)=b;
                p3=p3+1;
            end
            if test(i)==4
                force0_1(p0,:)=f1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force0_2(p0,:)=f2(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force0_3(p0,:)=f3(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force0_4(p0,:)=f4(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force0_5(p0,:)=f5(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force0_6(p0,:)=f6(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                
                ta_emg(p0,:)=ta(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                sol_emg(p0,:)=sol(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                pl_emg(p0,:)=pl(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                gca_emg(p0,:)=gca(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                weight4(p0,:)=w1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-w1(peaks(i)-WEIGHT_BASELINE_SHIFT);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p0_plat_torque(p0,:)=pert_torque(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-pert_torque(peaks(i)+50);
                p0_plat_pos(p0,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p0_foot_pos(p0,:)=foot_pos_data(peaks(i)-179+SHIFT_VAL:peaks(i)+2221+SHIFT_VAL);
                p0_phase(p0,:)=rigid_phase_tot(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                p0_pert(p0,:)=perturb_start(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                img0_pos(p0)=getmin(peaks(i),img_st,Img);
                cop4(p0,:)=cop(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                [a,b]=findpeaks(diff(p0_pert(p0,:)));

                p0_peakst(p0)=b;

                [a,b]=min(diff(p0_phase(p0,:)));

                p0_peakend(p0)=b;
                p0=p0+1;
            end
            if test(i)==5
                force4_1(p4,:)=f1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force4_2(p4,:)=f2(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force4_3(p4,:)=f3(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force4_4(p4,:)=f4(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force4_5(p4,:)=f5(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                force4_6(p4,:)=f6(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-f1(peaks(i)-FORCE_BASELINE_SHIFT);
                weight3(p4,:)=w1(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-w1(peaks(i)-WEIGHT_BASELINE_SHIFT);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p4_plat_torque(p4,:)=pert_torque(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES)-pert_torque(peaks(i)+50);
                p4_plat_pos(p4,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p4_foot_pos(p4,:)=foot_pos_data(peaks(i)-179+SHIFT_VAL:peaks(i)+2221+SHIFT_VAL);
                %p4_phase(p4,:)=haptic_phase_tot(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                p4_pert(p4,:)=perturb_start(peaks(i)-PRE_PURT_TRIG_SAMPLES:peaks(i)+POST_PURT_TRIG_SAMPLES);
                img4_pos(p4)=getmin(peaks(i),img_st,Img);
                [a,b]=findpeaks(diff(p4_pert(p4,:)));

                p4_peakst(p4)=b;
        
                p4=p4+1;
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


p0_er=0;
p1_er=0;
p2_er=0;
p3_er=0;
p4_er=0;

%% removing outliers using motion capture data: Its redundant but I left it as it wont really change anything
for i=1:p1-1
    if(img1_pos(i)>IMG_POS_CUTOFF_LIMIT)
        img1_pos(i)=NaN;
    end
end
for i=1:p2-1
    if(img2_pos(i)>IMG_POS_CUTOFF_LIMIT)
        img2_pos(i)=NaN;
    end
end
for i=1:p3-1
    if(img3_pos(i)>IMG_POS_CUTOFF_LIMIT)
        img3_pos(i)=NaN;
    end
end
for i=1:p0-1
    if(img0_pos(i)>IMG_POS_CUTOFF_LIMIT)
        img0_pos(i)=NaN;
    end
end


if PLOT_HISTOGRAM_FLAG==1
    im=[img0_pos img1_pos img2_pos img3_pos img4_pos];
end


%% Obtaining perturbation torque and CoP data from different trials the act_torque is calculated by using the foot position obtained form motion capture
for i=1:p1-1
    if (isnan(p1_plat_torque(i,1))==1)
        p1_act_torque(i,1:2201)=NaN;
        ssss='Found one';
        disp(ssss);
    else
        p1_act_torque(i,:)=force1_1(i,:)*(0.315-(img1_pos(i)/100))+force1_4(i,:)*(0.315-(img1_pos(i)/100))-force1_2(i,:)*(0.105+(img1_pos(i)/100))-force1_3(i,:)*(0.105+(img1_pos(i)/100))-force1_5(i,:)*(0.095)-force1_6(i,:)*(0.095);
    end
end
for i=1:p2-1
    if (isnan(p2_plat_torque(i,1))==1)
        p2_act_torque(i,1:2201)=NaN;
    else
        p2_act_torque(i,:)=force2_1(i,:)*(0.315-(img2_pos(i)/100))+force2_4(i,:)*(0.315-(img2_pos(i)/100))-force2_2(i,:)*(0.105+(img2_pos(i)/100))-force2_3(i,:)*(0.105+(img2_pos(i)/100))-force2_5(i,:)*(0.095)-force2_6(i,:)*(0.095);
    end
end


for i=1:p3-1
    if (isnan(p3_plat_torque(i,1))==1)
        p3_act_torque(i,1:2201)=NaN;
    else
        p3_act_torque(i,:)=force3_1(i,:)*(0.315-(img3_pos(i)/100))+force3_4(i,:)*(0.315-(img3_pos(i)/100))-force3_2(i,:)*(0.105+(img3_pos(i)/100))-force3_3(i,:)*(0.105+(img3_pos(i)/100))-force3_5(i,:)*(0.095)-force3_6(i,:)*(0.095);
    end
end
for i=1:p0-1
    if (isnan(p0_plat_torque(i,1))==1)
        p0_act_torque(i,1:2201)=NaN;
        p0_cop_torque(i,1:2201)=NaN;
    else
        p0_act_torque(i,:)=force0_1(i,:)*(0.315-(img0_pos(i)/100))+force0_4(i,:)*(0.315-(img0_pos(i)/100))-force0_2(i,:)*(0.105+(img0_pos(i)/100))-force0_3(i,:)*(0.105+(img0_pos(i)/100))-force0_5(i,:)*(0.095)-force0_6(i,:)*(0.095);
        p0_cop_torque(i,:)=force0_1(i,:)*(0.315-(img0_pos(i)/100))+force0_4(i,:)*(0.315-(img0_pos(i)/100))-force0_2(i,:)*(0.105+(img0_pos(i)/100))-force0_3(i,:)*(0.105+(img0_pos(i)/100))-force0_5(i,:)*(0.025)-force0_6(i,:)*(0.025);
        
    end
end
for i=1:p4-1
    if (isnan(p4_plat_torque(i,1))==1)
        p4_act_torque(i,1:2201)=NaN;
    else
        p4_act_torque(i,:)=force4_1(i,:)*(0.315-(img4_pos(i)/100))+force4_4(i,:)*(0.315-(img4_pos(i)/100))-force4_2(i,:)*(0.105+(img4_pos(i)/100))-force4_3(i,:)*(0.105+(img4_pos(i)/100))-force4_5(i,:)*(0.095)-force4_6(i,:)*(0.095);
    end
end




%% Calculation of means for plotting...not really used in analysis
weight1m=nanmean(weight1);
weight2m=nanmean(weight2);
weight3m=nanmean(weight3);
weight4m=trimmean(weight4,30);

cop4m=trimmean(p0_cop_torque,30);

p0_plat_torquem=trimmean(p0_act_torque,30);
p0_plat_posm=trimmean(p0_plat_pos,30);
p0_foot_posm=nanmean(p0_foot_pos);
p1_plat_torquem=trimmean(p1_act_torque,30);
p1_plat_posm=trimmean(p1_plat_pos,30);
p1_foot_posm=trimmean(p1_foot_pos,30);
p2_plat_torquem=trimmean(p2_act_torque,30);
p2_plat_posm=trimmean(p2_plat_pos,30);
p2_foot_posm=trimmean(p2_foot_pos,30);
p3_plat_torquem=trimmean(p3_act_torque,30);
p3_plat_posm=trimmean(p3_plat_pos,30);
p3_foot_posm=trimmean(p3_foot_pos,30);
p4_plat_torquem=trimmean(p4_act_torque,30);
p4_plat_posm=trimmean(p4_plat_pos,30);
p4_foot_posm=trimmean(p4_foot_pos,30);


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

%% Rigid and Haptic no perturbation torque;
% the torque is collected from PRE_PERT_SAMPLES samples prior
% to POST_PERT_SAMPLES samples post perturbation (heel strike)
%It is then subtracted from observed torque to get
%differential torque.
%Excluded uses min number of perturbations so that
%analysis is unifrom
excluded=[p1,p2,p3,p4];
analysis_value=min(excluded);
for i=1:p0-1
    point_i=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT1_PERC)+p0_peakst(i)
    p0_plat_torque15(i,:)=p0_act_torque(i,point_i-PRE_PERT_SAMPLES:point_i+POST_PERT_SAMPLES);
    p0_cop_torque15(i,:)=p0_cop_torque(i,point_i-PRE_PERT_SAMPLES:point_i+POST_PERT_SAMPLES);
    p0_foot_pos15(i,:)=p0_foot_pos(i,point_i-PRE_PERT_SAMPLES:point_i+POST_PERT_SAMPLES);
    
    point_j=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT2_PERC)+p0_peakst(i)
    p0_plat_torque30(i,:)=p0_act_torque(i,point_j-PRE_PERT_SAMPLES:point_j+POST_PERT_SAMPLES);
    p0_cop_torque30(i,:)=p0_cop_torque(i,point_j-PRE_PERT_SAMPLES:point_j+POST_PERT_SAMPLES);
    p0_foot_pos30(i,:)=p0_foot_pos(i,point_j-PRE_PERT_SAMPLES:point_j+POST_PERT_SAMPLES);
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT3_PERC)+p0_peakst(i)
    p0_plat_torque45(i,:)=p0_act_torque(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    p0_cop_torque45(i,:)=p0_cop_torque(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    p0_foot_pos45(i,:)=p0_foot_pos(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT4_PERC)+p0_peakst(i)
    p0_plat_torque60(i,:)=p0_act_torque(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    p0_cop_torque60(i,:)=p0_cop_torque(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    p0_foot_pos60(i,:)=p0_foot_pos(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
end

% Weight during perturbation
for i=1:p0-1
    point_i=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT1_PERC)+p0_peakst(i)
    weightr15(i,:)=weight4(i,point_i-PRE_PERT_SAMPLES:point_i+POST_PERT_SAMPLES);
    point_j=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT2_PERC)+p0_peakst(i)
    weightr30(i,:)=weight4(i,point_j-PRE_PERT_SAMPLES:point_j+POST_PERT_SAMPLES);
    point_k=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT3_PERC)+p0_peakst(i)
    weightr45(i,:)=weight4(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    point_k=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT4_PERC)+p0_peakst(i)
    weightr60(i,:)=weight4(i,point_k-PRE_PERT_SAMPLES:point_k+POST_PERT_SAMPLES);
    
end


%CoP during perturbation
for i=1:p0-1
 for j=1:length(weightr15)
     copr15(i,j)=p0_cop_torque15(i,j)/weightr15(i,j);
 end
for j=1:length(weightr30)
     copr30(i,j)=p0_cop_torque30(i,j)/weightr30(i,j);
end
for j=1:length(weightr45)
     copr45(i,j)=p0_cop_torque45(i,j)/weightr45(i,j);
end
 for j=1:length(weightr60)
     copr60(i,j)=p0_cop_torque60(i,j)/weightr60(i,j);
 end
    
end


for i=1:p0-1
    point_i=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT1_PERC)+p0_peakst(i)
    ta15(i)=mean(ta_emg(i,point_i-2:point_i+2));
    pl15(i)=mean(pl_emg(i,point_i-25:point_i+25));
    sol15(i)=mean(sol_emg(i,point_i-25:point_i+25));
    gca15(i)=mean(gca_emg(i,point_i-25:point_i+25));
   
    point_j=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT2_PERC)+p0_peakst(i)
    ta30(i)=mean(ta_emg(i,point_j-25:point_j+25));
    pl30(i)=mean(pl_emg(i,point_j-25:point_j+25));
    sol30(i)=mean(sol_emg(i,point_j-25:point_j+25));
    gca30(i)=mean(gca_emg(i,point_j-25:point_j+25));
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT3_PERC)+p0_peakst(i)
    ta45(i)=mean(ta_emg(i,point_k-25:point_k+25));
    pl45(i)=mean(pl_emg(i,point_k-25:point_k+25));
    sol45(i)=mean(sol_emg(i,point_k-25:point_k+25));
    gca45(i)=mean(gca_emg(i,point_k-25:point_k+25));
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*PERT_POINT4_PERC)+p0_peakst(i)
    ta60(i)=mean(ta_emg(i,point_k-25:point_k+25));
    pl60(i)=mean(pl_emg(i,point_k-25:point_k+25));
    sol60(i)=mean(sol_emg(i,point_k-25:point_k+25));
    gca60(i)=mean(gca_emg(i,point_k-25:point_k+25));
end

p0_plat_torque15m=mean(p0_plat_torque15);
p0_plat_torque30m=mean(p0_plat_torque30);
p0_plat_torque45m=mean(p0_plat_torque45);
p0_plat_torque60m=mean(p0_plat_torque60);
p0_foot_pos15m=mean(p0_foot_pos15);
p0_foot_pos30m=mean(p0_foot_pos30);
p0_foot_pos45m=mean(p0_foot_pos45);
p0_foot_pos60m=mean(p0_foot_pos60);


ta15m=mean(ta15);
pl15m=mean(pl15);
sol15m=mean(sol15);
gca15m=mean(gca15);
ta30m=mean(ta30);
pl30m=mean(pl30);
sol30m=mean(sol30);
gca30m=mean(gca30);
ta45m=mean(ta45);
pl45m=mean(pl45);
sol45m=mean(sol45);
gca45m=mean(gca45);
ta60m=mean(ta60);
pl60m=mean(pl60);
sol60m=mean(sol60);
gca60m=mean(gca60);
emg_final=[ta15m,pl15m,sol15m,gca15m;ta30m,pl30m,sol30m,gca30m;ta45m,pl45m,sol45m,gca45m;ta60m,pl60m,sol60m,gca60m;];
%% Obtaining differential data for impedance analysis
% data is broken into chunks of 400 points:PRE_PERT_SAMPLES
% before pert and POST_PERT_SAMPLES after pert. 

for i=1:analysis_value-1
    p1_foot_pos(i,:)=p1_foot_pos(i,:)-p0_foot_posm;
    p2_foot_pos(i,:)=p2_foot_pos(i,:)-p0_foot_posm;
    p3_foot_pos(i,:)=p3_foot_pos(i,:)-p0_foot_posm;
    p4_foot_pos(i,:)=p4_foot_pos(i,:)-p0_foot_posm;
    
    diff_p1_plat_pos(i,:)=p1_plat_pos(i,p1_peakst(i)-PRE_PERT_SAMPLES:p1_peakst(i)+POST_PERT_SAMPLES);
    diff_p2_plat_pos(i,:)=p2_plat_pos(i,p2_peakst(i)-PRE_PERT_SAMPLES:p2_peakst(i)+POST_PERT_SAMPLES);
    diff_p3_plat_pos(i,:)=p3_plat_pos(i,p3_peakst(i)-PRE_PERT_SAMPLES:p3_peakst(i)+POST_PERT_SAMPLES);
    diff_p1_plat_torque(i,:)=p1_act_torque(i,p1_peakst(i)-PRE_PERT_SAMPLES:p1_peakst(i)+POST_PERT_SAMPLES);
    diff_p2_plat_torque(i,:)=p2_act_torque(i,p2_peakst(i)-PRE_PERT_SAMPLES:p2_peakst(i)+POST_PERT_SAMPLES);
    diff_p3_plat_torque(i,:)=p3_act_torque(i,p3_peakst(i)-PRE_PERT_SAMPLES:p3_peakst(i)+POST_PERT_SAMPLES);
    
    diff_p1_plat_torque(i,:)=diff_p1_plat_torque(i,:)-p0_plat_torque30m;
    diff_p2_plat_torque(i,:)=diff_p2_plat_torque(i,:)-p0_plat_torque45m;
    diff_p3_plat_torque(i,:)=diff_p3_plat_torque(i,:)-p0_plat_torque60m;
    
%     diff_p1_foot_pos(i,:)=p1_foot_pos(i,p1_peakst(i)-PRE_PERT_SAMPLES:p1_peakst(i)+300)-p0_foot_pos30m;
%     diff_p2_foot_pos(i,:)=p2_foot_pos(i,p2_peakst(i)-PRE_PERT_SAMPLES:p2_peakst(i)+300)-p0_foot_pos45m;
%     diff_p3_foot_pos(i,:)=p3_foot_pos(i,p3_peakst(i)-PRE_PERT_SAMPLES:p3_peakst(i)+300)-p0_foot_pos60m;

    diff_p1_foot_pos(i,:)=p1_foot_pos(i,p1_peakst(i)-PRE_PERT_SAMPLES:p1_peakst(i)+POST_PERT_SAMPLES);
    diff_p2_foot_pos(i,:)=p2_foot_pos(i,p2_peakst(i)-PRE_PERT_SAMPLES:p2_peakst(i)+POST_PERT_SAMPLES);
    diff_p3_foot_pos(i,:)=p3_foot_pos(i,p3_peakst(i)-PRE_PERT_SAMPLES:p3_peakst(i)+POST_PERT_SAMPLES);
    
%      diff_p1_foot_pos(i,:)=diff_p1_foot_pos(i,:)-p0_foot_pos30m;
%     diff_p2_foot_pos(i,:)=diff_p2_foot_pos(i,:)-p0_foot_pos45m;
%     diff_p3_foot_pos(i,:)=diff_p3_foot_pos(i,:)-p0_foot_pos60m;
    diff_p1_foot_pos(i,:)=diff_p1_foot_pos(i,:)-diff_p1_foot_pos(i,100);
    diff_p2_foot_pos(i,:)=diff_p2_foot_pos(i,:)-diff_p2_foot_pos(i,100);
    diff_p3_foot_pos(i,:)=diff_p3_foot_pos(i,:)-diff_p3_foot_pos(i,100);
   

    
    diff_p1_plat_pos(i,:)=diff_p1_plat_pos(i,:)-diff_p1_plat_pos(i,100);
    diff_p2_plat_pos(i,:)=diff_p2_plat_pos(i,:)-diff_p2_plat_pos(i,100);
    diff_p3_plat_pos(i,:)=diff_p3_plat_pos(i,:)-diff_p3_plat_pos(i,100);
    
    diff_p1_plat_vel(i,1)=0;
    for l=2:length(diff_p1_plat_pos(1,:))
        diff_p1_plat_vel(i,l)=(diff_p1_plat_pos(i,l)-diff_p1_plat_pos(i,l-1))/0.0005;
    end
    diff_p2_plat_vel(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p2_plat_vel(i,l)=(diff_p2_plat_pos(i,l)-diff_p2_plat_pos(i,l-1))/0.0005;
    end
    diff_p3_plat_vel(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p3_plat_vel(i,l)=(diff_p3_plat_pos(i,l)-diff_p3_plat_pos(i,l-1))/0.0005;
    end
    diff_p1_foot_vel(i,1)=0;
    for l=2:length(diff_p1_foot_pos(1,:))
        diff_p1_foot_vel(i,l)=(diff_p1_foot_pos(i,l)-diff_p1_foot_pos(i,l-1))/0.0005;
    end
    diff_p2_foot_vel(i,1)=0;
    for l=2:length(diff_p2_foot_pos(1,:))
        diff_p2_foot_vel(i,l)=(diff_p2_foot_pos(i,l)-diff_p2_foot_pos(i,l-1))/0.0005;
    end
    diff_p3_foot_vel(i,1)=0;
    for l=2:length(diff_p3_foot_pos(1,:))
        diff_p3_foot_vel(i,l)=(diff_p3_foot_pos(i,l)-diff_p3_foot_pos(i,l-1))/0.0005;
    end
    
      diff_p1_plat_acc(i,1)=0;
    for l=2:length(diff_p1_plat_pos(1,:))
        diff_p1_plat_acc(i,l)=(diff_p1_plat_vel(i,l)-diff_p1_plat_vel(i,l-1))/0.0005;
    end
    diff_p2_plat_acc(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p2_plat_acc(i,l)=(diff_p2_plat_vel(i,l)-diff_p2_plat_vel(i,l-1))/0.0005;
    end
    diff_p3_plat_acc(i,1)=0;
    for l=2:length(diff_p3_plat_pos(1,:))
        diff_p3_plat_acc(i,l)=(diff_p3_plat_vel(i,l)-diff_p3_plat_vel(i,l-1))/0.0005;
    end
    diff_p1_foot_acc(i,1)=0;
    for l=2:length(diff_p1_plat_pos(1,:))
        diff_p1_foot_acc(i,l)=(diff_p1_foot_vel(i,l)-diff_p1_foot_vel(i,l-1))/0.0005;
    end
    diff_p2_foot_acc(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p2_foot_acc(i,l)=(diff_p2_foot_vel(i,l)-diff_p2_foot_vel(i,l-1))/0.0005;
    end
    diff_p3_foot_acc(i,1)=0;
    for l=2:length(diff_p3_plat_pos(1,:))
        diff_p3_foot_acc(i,l)=(diff_p3_foot_vel(i,l)-diff_p3_foot_vel(i,l-1))/0.0005;
    end  
   
    
    diff_p4_plat_pos(i,:)=p4_plat_pos(i,p4_peakst(i)-PRE_PERT_SAMPLES:p4_peakst(i)+POST_PERT_SAMPLES);
    diff_p4_plat_torque(i,:)=p4_act_torque(i,p4_peakst(i)-PRE_PERT_SAMPLES:p4_peakst(i)+POST_PERT_SAMPLES);
    
    diff_p4_plat_torque(i,:)=diff_p4_plat_torque(i,:)-p0_plat_torque15m;
    
%     diff_p4_foot_pos(i,:)=p4_foot_pos(i,p4_peakst(i)-PRE_PERT_SAMPLES:p4_peakst(i)+POST_PERT_SAMPLES)-p0_foot_pos15m;
    diff_p4_foot_pos(i,:)=p4_foot_pos(i,p4_peakst(i)-PRE_PERT_SAMPLES:p4_peakst(i)+POST_PERT_SAMPLES);
    diff_p4_foot_pos(i,:)=diff_p4_foot_pos(i,:)-diff_p4_foot_pos(i,100);
    
    
    diff_p4_plat_pos(i,:)=diff_p4_plat_pos(i,:)-diff_p4_plat_pos(i,100);
   
    diff_p4_plat_vel(i,1)=0;
    for l=2:length(diff_p4_plat_pos(1,:))
        diff_p4_plat_vel(i,l)=(diff_p4_plat_pos(i,l)-diff_p4_plat_pos(i,l-1))/0.0005;
    end
    diff_p4_foot_vel(i,1)=0;
    for l=2:length(diff_p4_foot_pos(1,:))
        diff_p4_foot_vel(i,l)=(diff_p4_foot_pos(i,l)-diff_p4_foot_pos(i,l-1))/0.0005;
    end
    
    
      diff_p4_plat_acc(i,1)=0;
    for l=2:length(diff_p4_plat_pos(1,:))
        diff_p4_plat_acc(i,l)=(diff_p4_plat_vel(i,l)-diff_p4_plat_vel(i,l-1))/0.0005;
    end
    
 
    diff_p4_foot_acc(i,1)=0;
    for l=2:length(diff_p4_plat_pos(1,:))
        diff_p4_foot_acc(i,l)=(diff_p4_foot_vel(i,l)-diff_p4_foot_vel(i,l-1))/0.0005;
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

diff_p4_plat_accm=trimmean(diff_p4_plat_acc,30);
diff_p4_foot_accm=trimmean(diff_p4_foot_acc,30);


%% Removing platform dynamics and Impedance estimation

%Estimated platform dynamics are removed and
%individual impances are estimated.
exc_rigid=[p1,p2,p3,p4];
analysis_value=min(exc_rigid);
for i=1:analysis_value-1
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torque(i,:)-0.1945*diff_p1_plat_acc(i,:);%+0.02*diff_p1_foot_accm;%+1*diff_p1_foot_velm;
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torqueimp(i,:)-diff_p1_plat_torqueimp(i,100);
    
   
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torque(i,:)-0.1945*diff_p2_plat_acc(i,:);%+0.013*diff_p2_foot_accm;
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torqueimp(i,:)-diff_p2_plat_torqueimp(i,100);

    
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torque(i,:)-0.1945*diff_p3_plat_acc(i,:);%-0.09*diff_p3_foot_accm;
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torqueimp(i,:)-diff_p3_plat_torqueimp(i,100);
   
end





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
% constained linear regression
%Value is made NaN if foot position is not on
%platform
for i=1:analysis_value-1
    C1=[diff_p1_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p1_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p1_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d1=diff_p1_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B1=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    if(isnan(d1(1))==0)
        p1imp(i,:)=lsqlin(C1,d1,A1,B1);
        
    else
        p1imp(i,:)=[NaN NaN NaN];
    end
    vartor3=var(diff_p1_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX));
    varimp3=var(diff_p1_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p1_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p1imp(i,1)+diff_p1_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p1imp(i,2)+diff_p1_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p1imp(i,3)));
    p1impg(i)=100*(1-(varimp3/vartor3));
end

p1impm=regress(diff_p1_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)',[diff_p1_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p1_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p1_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' ]);
C=[diff_p1_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p1_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p1_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
d=diff_p1_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
p1impm2=lsqlin(C,d,A,B);


for i=1:analysis_value-1
    
    C1=[diff_p2_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p2_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p2_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d1=diff_p2_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B1=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    
    if(isnan(d1(1))==0)
        p2imp(i,:)=lsqlin(C1,d1,A1,B1);
        
    else
        p2imp(i,:)=[NaN NaN NaN];
    end
    
    vartor3=var(diff_p2_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX));
    varimp3=var(diff_p2_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p2_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p2imp(i,1)+diff_p2_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p2imp(i,2)+diff_p2_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p2imp(i,3)));
    p2impg(i)=100*(1-(varimp3/vartor3));
    
end
p2impm=regress(diff_p2_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)',[diff_p2_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p2_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p2_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)']);
C=[diff_p2_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p2_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p2_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
d=diff_p2_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
p2impm2=lsqlin(C,d,A,B);

for i=1:analysis_value-1
    
    C1=[diff_p3_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p3_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p3_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d1=diff_p3_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B1=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    
    if(isnan(d1(1))==0)
        p3imp(i,:)=lsqlin(C1,d1,A1,B1);
        
    else
        p3imp(i,:)=[NaN NaN NaN];
    end
    
    vartor3=var(diff_p3_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX));
    varimp3=var(diff_p3_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p3_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p3imp(i,1)+diff_p3_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p3imp(i,2)+diff_p3_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p3imp(i,3)));
    p3impg(i)=100*(1-(varimp3/vartor3));
end

p3impm=regress(diff_p3_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)',[diff_p3_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p3_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p3_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)']);

C=[diff_p3_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p3_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p3_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
d=diff_p3_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
p3impm2=lsqlin(C,d,A,B);

exc_haptic=[p4];
% analysis_value=min(exc_haptic);

for i=1:analysis_value-1
    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torque(i,:)-0.1945*diff_p4_plat_acc(i,:);%+0.02*diff_p4_foot_accm;%+1*diff_p4_foot_velm;
    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torqueimp(i,:)-diff_p4_plat_torqueimp(i,100);


end




diff_p4_foot_accm=diff_p4_foot_accm;
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);

diff_p4_foot_posm=trimmean(diff_p4_foot_pos,30);
diff_p4_foot_velm=trimmean(diff_p4_foot_vel,30);
diff_p4_plat_velm=trimmean(diff_p4_plat_vel,30);
diff_p4_foot_accm=trimmean(diff_p4_foot_acc,30);


for i=1:analysis_value-1
    
    %     p4imp(i,:)=regress(diff_p4_plat_torqueimp(i,100:REGRESSION_END_INEX)',[diff_p4_foot_pos(i,100:REGRESSION_END_INEX)' diff_p4_foot_vel(i,100:REGRESSION_END_INEX)' diff_p4_foot_acc(i,100:REGRESSION_END_INEX)']);
    C1=[diff_p4_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p4_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p4_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d1=diff_p4_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B1=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    
    if(isnan(d1(1))==0)
        p4imp(i,:)=lsqlin(C1,d1,A1,B1);
        
    else
        p4imp(i,:)=[NaN NaN NaN];
    end
    
    vartor3=var(diff_p4_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX));
    varimp4=var(diff_p4_plat_torqueimp(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p4_foot_pos(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p4imp(i,1)+diff_p4_foot_vel(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p4imp(i,2)+diff_p4_foot_acc(i,REGRESSION_START_INDEX:REGRESSION_END_INDEX)*p4imp(i,3)));
    p4impg(i)=100*(1-(varimp4/vartor3));
    
end
p4impm=regress(diff_p4_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)',[diff_p4_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p4_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p4_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)']);
C=[diff_p4_foot_posm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p4_foot_velm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' diff_p4_foot_accm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
d=diff_p4_plat_torqueimpm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
p4impm2=lsqlin(C,d,A,B);



%% Storing final data
% impedance=[p1impm';p2impm';p3impm';p4impm';p5impm';p6impm'];
% impedance_constrained=[p1impm2';p2impm2';p3impm2'...
%     ;p4impm2';p5impm2';p6impm2'];
% goodness;
% goodness_constrained;
% Final_value=[impedance, goodnessn];
% save('impedance.mat','Final_value');
% Final_value_constrained=[impedance_constrained, goodnessc];
% save('impedance_constrained.mat','Final_value_constrained');
% if(plot_torque==1)
%     torque_comparison_plot;
% end
% 
% if(plot_hist==1)
%     figure
%     hist(im,40);
%     saveas(gcf,'histogram.jpg');
% end
% 

% if(plot_figs_constrained==1)
%     fit_figures_constrained;
% end
% plottingemgweight;

individual_impedance
best_boot_no_avg
if(PLOT_FIGS_FLAG==1)
    %fit_figures_bootstrapping;
    plot_regression_results(diff_p1_foot_posv, diff_p1_foot_velv, diff_p1_foot_accv, diff_p1_plat_torqueimpv, ...
        boot_impv(1, 1:3), '31%');
    plot_regression_results(diff_p2_foot_posv, diff_p2_foot_velv, diff_p2_foot_accv, diff_p2_plat_torqueimpv, ...
        boot_impv(2, 1:3), '44%');
    plot_regression_results(diff_p3_foot_posv, diff_p3_foot_velv, diff_p3_foot_accv, diff_p3_plat_torqueimpv, ...
        boot_impv(3, 1:3), '57%');
    plot_regression_results(diff_p4_foot_posv, diff_p4_foot_velv, diff_p4_foot_accv, diff_p4_plat_torqueimpv, ...
        boot_impv(4, 1:3), '18%');
end