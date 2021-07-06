%----mastercode_onlyrigid_br.m----%
% Same as mastercode_onlyrigid.m except no bootstrapping. Just basic
% regression with outlier removal

%% Setup, Filtering and Data Visualization settings
close all
clearvars -except GROUP_DATA_FOLDER_REL_LOC multi_subj SUBJ_DATA_DIRS CURR_SUBJ_REL_LOC


%%Check to see if this is the main script or is being called by
%%mastercode_multiple_subjects.m
ise = evalin( 'base', 'exist(''multi_subj'',''var'') == 1' )
if ise == 1
   DATA_FOLDER_REL_LOC = CURR_SUBJ_REL_LOC; 

else
% Insert subject initial and name.
%Make sure it matches the format for naming
    DATA_FOLDER_REL_LOC = "./../../data/Vu_121420/" %Relative location for current code dir.
end

path_str = split(DATA_FOLDER_REL_LOC, ["/", "_"]);
sub_name =  path_str{5}; % 5 since the relative locations are considered parth of the path
sub_initial = sub_name(1);

RESULTS_DIR = strcat('./results/', sub_name, '/');
mkdir(RESULTS_DIR);
%Add number of perturbation you actually ran
num_pert=40;


BOOTSTRAP_COUNT = 100;

SAMPLE_RATE_HZ = 2000;
SAMPLE_PERIOD = 1/SAMPLE_RATE_HZ;


VOLTS_TO_NEWTONS_SCALER = 53.4; 
F1_SIG = 9;
F2_SIG = 10;
F3_SIG = 11;
F4_SIG = 12;
F5_SIG = 20;
F6_SIG = 21;

TA_EMG_SIG  = 1;
SOL_EMG_SIG = 2;
PL_EMG_SIG  = 3;
GCA_EMG_SIG = 4;

MOCAP_SAMPLE_INDEX_CONV_FACTOR = 20;
MOCAP_OUTLIER_LIMIT = 500;
IMG_ENABLED_SIG = 6;

PERT_TORQUE_SIG = 7;

FOOT_GON_POS_SIG = 13;
PLAT_GON_POS_SIG = 14;

RIGID_PHASE_PRED_SIG  = 15;
HAPTIC_PHASE_PRED_SIG = 16;

TRIAL_TYPE_FLAG = 17; %Also indicates start of trial

WEIGHT_SIG = 18;
COP_SIG = 19;

PERT_START_SIG = 22;


% insert lower limit of inertia of foot in the fit
% u_lim is the upper limit of the inertia and lim
% is the lower limit
lim= 0.007;
u_lim=0.02;
% change these flags to 1 for figures (normal fit and constrained fit)
plot_figs=1;
% plot_figs_constrained=0;
% Change to 1 to get ankle position histogram
plot_hist=1;
%  Change to get torque plot comparison figure
plot_torque=1;
shift=0;
%add the number of loops you want to run bootstrapping here
loops=100;

%Chunks intial data from the data files for each trial.
TRIAL_WINDOW_PRE_PERT  = -400;
TRIAL_WINDOW_POST_PERT = 2000;
TRIAL_BASELINE_INDEX   = -360; %Taken before the perturbation

%Data window surrounding the perturbation
PRE_PERT_WINDOW = -100;
POST_PERT_WINDOW = 300;

%Data window used for the regression
REGRESSION_WINDOW_MIN_INDEX = 100;
REGRESSION_WINDOW_MAX_INDEX = 300;

GON_DELAY_IN_SAMPLES = 221;

%Perturbation points are given in terms of the stance phase percentage.
PERT_POINT_1 = 0.18;
PERT_POINT_2 = 0.31;
PERT_POINT_3 = 0.44;
PERT_POINT_4 = 0.57;


outlier_criterion_std = 3.0; 

outlier_removal_type = "POS_TRQ";
%outlier_removal_type = "POS"; % POS | ACC_TOR
%outlier_removal_type = "ACC_TOR";
%outlier_removal_type = "POS_VEL_ACC_TRQ";
%outlier_removal_type = "NONE";
POS_REJECTION_LIMIT = 2.5;

VAF_REMOVAL_CRITERION = 80;  

% Add the trials you want to exclude in here
exclude=[];%exclude=[1,2];
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',...
    5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',15,...
    'DesignMethod','butter','Samplerate',2000);

% Getting the MVC values
mvc_evaluation;
%% Section to calculate goniometer gains
t=gonio_values_func(DATA_FOLDER_REL_LOC);
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
    
    if(ismember(trials,exclude)==0)
        if(trials<10)
            h = fopen(strcat(DATA_FOLDER_REL_LOC, sub_initial,'W0',num2str(trials),'.dat'));
        else
            h = fopen(strcat(DATA_FOLDER_REL_LOC, sub_initial,'W',num2str(trials),'.dat'));
        end
        
        live_data=fread(h);
        Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
        siz=size(Input1.data);
        Img_flag=Input1.data(:,IMG_ENABLED_SIG);
        [x,img_st]=findpeaks(diff(Img_flag));
        img_st=round(img_st(1)/MOCAP_SAMPLE_INDEX_CONV_FACTOR);
        if(trials<10)
            Img=csvread(strcat(DATA_FOLDER_REL_LOC, sub_name,'_00',num2str(trials),'.csv'));
        else
            Img=csvread(strcat(DATA_FOLDER_REL_LOC, sub_name,'_0',num2str(trials),'.csv'));
        end
        % obtaining data from channels
        pert_torque=filtfilt(d1,Input1.data(:,PERT_TORQUE_SIG));
        f1=Input1.data(:,F1_SIG)*VOLTS_TO_NEWTONS_SCALER;
        f1=filtfilt(d1,f1);
        f2=Input1.data(:,F2_SIG)*VOLTS_TO_NEWTONS_SCALER;
        f2=filtfilt(d1,f2);
        f3=Input1.data(:,F3_SIG)*VOLTS_TO_NEWTONS_SCALER;
        f3=filtfilt(d1,f3);
        f4=Input1.data(:,F4_SIG)*VOLTS_TO_NEWTONS_SCALER;
        f4=filtfilt(d1,f4);
        f5=Input1.data(:,F5_SIG)*VOLTS_TO_NEWTONS_SCALER/2;
        f5=filtfilt(d1,f5);
        f6=Input1.data(:,F6_SIG)*VOLTS_TO_NEWTONS_SCALER/2;
        f6=filtfilt(d1,f6);
        ta=Input1.data(:,TA_EMG_SIG);
        ta=abs(ta-off_TA)*100/mvc_ta;
        sol=Input1.data(:,SOL_EMG_SIG);
        sol=abs(sol-off_SOL)*100/mvc_sol;
        pl=Input1.data(:,PL_EMG_SIG);
        pl=abs(pl-off_PL)*100/mvc_pl;
        gca=Input1.data(:,GCA_EMG_SIG);
        gca=abs(gca-off_GCA)*100/mvc_gca;
        w1=filtfilt(d1,Input1.data(:,WEIGHT_SIG));
        cop=filtfilt(d1,Input1.data(:,COP_SIG));
        flag=Input1.data(:,TRIAL_TYPE_FLAG);
        rigid_phase_tot=Input1.data(:,RIGID_PHASE_PRED_SIG);
        haptic_phase_tot=Input1.data(:,HAPTIC_PHASE_PRED_SIG);
        perturb_start=Input1.data(:,PERT_START_SIG);
        foot_pos_data=filtfilt(d1,Input1.data(:,FOOT_GON_POS_SIG));
        foot_pos_data=((foot_pos_data-mean(foot_pos_data))...
            *DP_foot_gonio*pi/180);
        plat_pos_data=filtfilt(d1,Input1.data(:,PLAT_GON_POS_SIG));
        plat_pos_data=((plat_pos_data-mean(plat_pos_data))...
            *DP_plat_gonio*pi/180);
        % 17 records an impulse everytime a
        % perturbation occurs with diff amplitudes
        [test,peaks]=findpeaks(Input1.data(:,TRIAL_TYPE_FLAG));
        for i=1:length(peaks)
     % Each test value corresponds to a different pert       
            time=[-200:0.5:1000];
            if test(i)==1
                force1_1(p1,:)=f1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force1_2(p1,:)=f2(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force1_3(p1,:)=f3(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force1_4(p1,:)=f4(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force1_5(p1,:)=f5(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force1_6(p1,:)=f6(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                
                
                p1_EMG.TA(p1, :)  = ta(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p1_EMG.SOL(p1, :) = sol(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p1_EMG.PL(p1, :)  = pl(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p1_EMG.GCA(p1, :) = gca(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                weight1(p1,:)=w1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-w1(peaks(i)+TRIAL_BASELINE_INDEX);%+w2(peaks(i)-400:peaks(i)+TRIAL_WINDOW_POST_PERT)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p1_plat_torque(p1,:)=pert_torque(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-pert_torque(peaks(i)+50);
                p1_plat_pos(p1,:)=plat_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES);
                p1_foot_pos(p1,:)=foot_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES+shift:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES+shift);
                p1_phase(p1,:)=rigid_phase_tot(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p1_pert(p1,:)=perturb_start(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                img1_pos(p1)=getmin(peaks(i),img_st,Img);
                cop1(p1,:)=cop(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                [a,b]=findpeaks(abs(diff(p1_pert(p1,:))));
                p1_peakst(p1)=find(p1_phase(p1, :) > 31, 1);%min(b);
                p1=p1+1;
                
            end
            if test(i)==2
                force2_1(p2,:)=f1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force2_2(p2,:)=f2(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force2_3(p2,:)=f3(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force2_4(p2,:)=f4(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force2_5(p2,:)=f5(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force2_6(p2,:)=f6(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                
                
                p2_EMG.TA(p2, :)  = ta(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p2_EMG.SOL(p2, :) = sol(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p2_EMG.PL(p2, :)  = pl(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p2_EMG.GCA(p2, :) = gca(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                weight2(p2,:)=w1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-w1(peaks(i)+TRIAL_BASELINE_INDEX);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p2_plat_torque(p2,:)=pert_torque(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-pert_torque(peaks(i)+50);
                p2_plat_pos(p2,:)=plat_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES);
                p2_foot_pos(p2,:)=foot_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES+shift:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES+shift);
                p2_phase(p2,:)=rigid_phase_tot(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p2_pert(p2,:)=perturb_start(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                img2_pos(p2)=getmin(peaks(i),img_st,Img);
                cop2(p2,:)=cop(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                [a,b]=findpeaks(abs(diff(p2_pert(p2,:))));

                p2_peakst(p2)= find(p2_phase(p2, :) > 44, 1); %b
                p2=p2+1;
            end
            if test(i)==3
                force3_1(p3,:)=f1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force3_2(p3,:)=f2(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force3_3(p3,:)=f3(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force3_4(p3,:)=f4(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force3_5(p3,:)=f5(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force3_6(p3,:)=f6(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                
                p3_EMG.TA(p3,:)  = ta(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p3_EMG.SOL(p3,:) = sol(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p3_EMG.PL(p3,:)  = pl(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p3_EMG.GCA(p3,:) = gca(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                weight3(p3,:)=w1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-w1(peaks(i)+TRIAL_BASELINE_INDEX);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p3_plat_torque(p3,:)=1*(pert_torque(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-pert_torque(peaks(i)+50));
                p3_plat_pos(p3,:)=plat_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES);
                p3_foot_pos(p3,:)=foot_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES+shift:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES+shift);
                p3_phase(p3,:)=rigid_phase_tot(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p3_pert(p3,:)=perturb_start(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                img3_pos(p3)=getmin(peaks(i),img_st,Img);
                cop3(p3,:)=cop(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                [a,b]=findpeaks(abs(diff(p3_pert(p3,:))));
                p3_peakst(p3)=find(p3_phase(p3, :) > 57, 1);%b;
                p3=p3+1;
            end
            if test(i)==4
                force0_1(p0,:)=f1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force0_2(p0,:)=f2(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force0_3(p0,:)=f3(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force0_4(p0,:)=f4(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force0_5(p0,:)=f5(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force0_6(p0,:)=f6(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                
                ta_emg(p0,:)=ta(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                sol_emg(p0,:)=sol(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                pl_emg(p0,:)=pl(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                gca_emg(p0,:)=gca(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                weight0(p0,:)=w1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-w1(peaks(i)+TRIAL_BASELINE_INDEX);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p0_plat_torque(p0,:)=pert_torque(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-pert_torque(peaks(i)+50);
                p0_plat_pos(p0,:)=plat_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES);
                p0_foot_pos(p0,:)=foot_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES+shift:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES+shift);
                p0_phase(p0,:)=rigid_phase_tot(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p0_pert(p0,:)=perturb_start(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                img0_pos(p0)=getmin(peaks(i),img_st,Img);
                cop4(p0,:)=cop(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                [a,b]=findpeaks(abs(diff(p0_pert(p0,:))));

                p0_peakst(p0)=400;%b;

                [a,b]=min(diff(p0_phase(p0,:)));

                p0_peakend(p0)=b;
                p0=p0+1;
            end
            if test(i)==5
                force4_1(p4,:)=f1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force4_2(p4,:)=f2(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force4_3(p4,:)=f3(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force4_4(p4,:)=f4(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force4_5(p4,:)=f5(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                force4_6(p4,:)=f6(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-f1(peaks(i)+TRIAL_BASELINE_INDEX);
                
                p4_EMG.TA(p4,:)  = ta(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p4_EMG.SOL(p4,:) = sol(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p4_EMG.PL(p4,:)  = pl(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p4_EMG.GCA(p4,:) = gca(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                weight4(p4,:)=w1(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-w1(peaks(i)+TRIAL_BASELINE_INDEX);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p4_plat_torque(p4,:)=pert_torque(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT)-pert_torque(peaks(i)+50);
                p4_plat_pos(p4,:)=plat_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES);
                p4_foot_pos(p4,:)=foot_pos_data(peaks(i)+TRIAL_WINDOW_PRE_PERT+GON_DELAY_IN_SAMPLES+shift:peaks(i)+TRIAL_WINDOW_POST_PERT+GON_DELAY_IN_SAMPLES+shift);
                p4_phase(p4,:)=haptic_phase_tot(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                p4_pert(p4,:)=perturb_start(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
                img4_pos(p4)=getmin(peaks(i),img_st,Img);
                [a,b]=findpeaks(abs(diff(p4_pert(p4,:))));

                p4_peakst(p4)=find(p4_phase(p4, :) > 18, 1);%b;
        
                p4=p4+1;
            end
            
        end
        
%         fclose all;
%         clear live_data Input1 peaks Img;
        
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
img0_pos(img0_pos > MOCAP_OUTLIER_LIMIT) = NaN;
img1_pos(img1_pos > MOCAP_OUTLIER_LIMIT) = NaN;
img2_pos(img2_pos > MOCAP_OUTLIER_LIMIT) = NaN;
img3_pos(img3_pos > MOCAP_OUTLIER_LIMIT) = NaN;
img4_pos(img4_pos > MOCAP_OUTLIER_LIMIT) = NaN;

if plot_hist==1
    im=[img0_pos img1_pos img2_pos img3_pos img4_pos];
end

p0_removed_ind = (img0_pos < -0.5)' | (img0_pos > POS_REJECTION_LIMIT)';
p0_raw_data_ind = find(~p0_removed_ind);
p0_total_accepted = sum(~p0_removed_ind)



%% Obtaining perturbation torque and CoP data from different trials the act_torque is calculated by using the foot position obtained form motion capture

%The way this is implemented is only temporary, though the idea will remain
%the same: Find out what data needs to be removed, then create all
%proceeding variables from that selected data only. This is only being done
%currently for the No-Perturbation case, as it affects everything along the
%way. 

force_plate_data_p0.F1 = force0_1(p0_raw_data_ind, :);
force_plate_data_p0.F2 = force0_2(p0_raw_data_ind, :);
force_plate_data_p0.F3 = force0_3(p0_raw_data_ind, :);
force_plate_data_p0.F4 = force0_4(p0_raw_data_ind, :);
force_plate_data_p0.F5 = force0_5(p0_raw_data_ind, :);
force_plate_data_p0.F6 = force0_6(p0_raw_data_ind, :);

force_plate_data_p1.F1 = force1_1;
force_plate_data_p1.F2 = force1_2;
force_plate_data_p1.F3 = force1_3;
force_plate_data_p1.F4 = force1_4;
force_plate_data_p1.F5 = force1_5;
force_plate_data_p1.F6 = force1_6;

force_plate_data_p2.F1 = force2_1;
force_plate_data_p2.F2 = force2_2;
force_plate_data_p2.F3 = force2_3;
force_plate_data_p2.F4 = force2_4;
force_plate_data_p2.F5 = force2_5;
force_plate_data_p2.F6 = force2_6;

force_plate_data_p3.F1 = force3_1;
force_plate_data_p3.F2 = force3_2;
force_plate_data_p3.F3 = force3_3;
force_plate_data_p3.F4 = force3_4;
force_plate_data_p3.F5 = force3_5;
force_plate_data_p3.F6 = force3_6;

force_plate_data_p4.F1 = force4_1;
force_plate_data_p4.F2 = force4_2;
force_plate_data_p4.F3 = force4_3;
force_plate_data_p4.F4 = force4_4;
force_plate_data_p4.F5 = force4_5;
force_plate_data_p4.F6 = force4_6;

[p0_act_torque, p0_cop_torque] = get_act_and_cop_torque(img0_pos(p0_raw_data_ind),force_plate_data_p0);
[p1_act_torque, p1_cop_torque] = get_act_and_cop_torque(img1_pos,force_plate_data_p1);
[p2_act_torque, p2_cop_torque] = get_act_and_cop_torque(img2_pos,force_plate_data_p2);
[p3_act_torque, p3_cop_torque] = get_act_and_cop_torque(img3_pos,force_plate_data_p3);
[p4_act_torque, p4_cop_torque] = get_act_and_cop_torque(img4_pos,force_plate_data_p4);

TEMP_ANKLE_Y = 7.2;

max_weight = max(weight0, [], 'all');
for i = 1:size(p0_raw_data_ind, 1)
    start_index = min(find(weight0(p0_raw_data_ind(i),:) > 0.015*max_weight))
    end_index = max(find(weight0(p0_raw_data_ind(i),:) > 0.015*max_weight))
    start_index_vec(i) = start_index;
    end_index_vec(i) = end_index;
end

figure();
for i = 1:size(p0_raw_data_ind, 1)
    plot(weight0(p0_raw_data_ind(i), start_index_vec(i):end_index_vec(i))'); hold on;
end
hold off;

p0_sample_length = end_index_vec-start_index_vec;
stance_phase_duration = p0_sample_length * (1/2000);

cop_total = zeros(size(p0_raw_data_ind, 1), max(p0_sample_length)+1);
ankle_angle_total = zeros(size(p0_raw_data_ind, 1), max(p0_sample_length)+1);
weight_total = zeros(size(p0_raw_data_ind, 1), max(p0_sample_length)+1);
for i = 1:size(p0_raw_data_ind, 1)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   cop_torque_i = p0_cop_torque(i, start_index_vec(i):end_index_vec(i));
   weight_i = weight0(p0_raw_data_ind(i), start_index_vec(i):end_index_vec(i));
   cop_i = interp1(normalized_time_i', cop_torque_i./weight_i, 0:1/max(p0_sample_length):1);
   
   ankle_angle_i = p0_foot_pos(p0_raw_data_ind(i), start_index_vec(i):end_index_vec(i));
   ankle_angle_i = interp1(normalized_time_i', ankle_angle_i, 0:1/max(p0_sample_length):1);
   
   weight_i = interp1(normalized_time_i', weight_i, 0:1/max(p0_sample_length):1);
   
   cop_total(i, :) = cop_i;
   ankle_angle_total(i, :) = ankle_angle_i;
   weight_total(i, :) = weight_i;
end


plot(cop_total, -TEMP_ANKLE_Y*ones(size(cop_total)))
p_AF_x = zeros(size(cop_total));
p_AF_y = zeros(size(cop_total));
for trial_index = 1:size(p0_raw_data_ind, 1)
   for time_index = 1:size(ankle_angle_total, 2)
       theta = ankle_angle_total(trial_index, time_index);
       R = [cos(theta), -sin(theta);
            sin(theta),  cos(theta)];
   
       p_F = [100*cop_total(trial_index, time_index); - TEMP_ANKLE_Y];
       p_AF = R*p_F;
       
       p_AF_x(trial_index, time_index) = p_AF(1);
       p_AF_y(trial_index, time_index) = p_AF(2);
   
   end
end

figure();
plot(p_AF_x(:,                1:round(0.18*1554))', p_AF_y(:, 1:round(0.18*1554))', 'k'); hold on;
plot(p_AF_x(:, round(0.18*1554):round(0.33*1554))', p_AF_y(:, round(0.18*1554):round(0.33*1554))', 'r'); hold on;
plot(p_AF_x(:, round(0.33*1554):round(0.44*1554))', p_AF_y(:, round(0.33*1554):round(0.44*1554))', 'g'); hold on;
plot(p_AF_x(:, round(0.44*1554):round(0.57*1554))', p_AF_y(:, round(0.44*1554):round(0.57*1554))', 'b'); hold on;
plot(p_AF_x(:, round(0.57*1554):end)',              p_AF_y(:, round(0.57*1554):end)', 'k');
hold off;

figure();
plot(1:round(0.18*1554), ankle_angle_total(:, 1:round(0.18*1554))', 'k'); hold on;
plot(round(0.18*1554):round(0.33*1554) ,ankle_angle_total(:, round(0.18*1554):round(0.33*1554))', 'r'); hold on;
plot(round(0.33*1554):round(0.44*1554) ,ankle_angle_total(:, round(0.33*1554):round(0.44*1554))', 'g'); hold on;
plot(round(0.44*1554):round(0.57*1554) ,ankle_angle_total(:, round(0.44*1554):round(0.57*1554))', 'b'); hold on;
plot(round(0.57*1554):1554 ,ankle_angle_total(:, round(0.57*1554):end)', 'k');
hold off;

R_vec = [];
x_vec = [];
y_vec = [];

for trial_index = 1:size(p0_raw_data_ind, 1)
    [xc,yc,R,a] = circfit(p_AF_x(trial_index, round(0.018*1554):round(0.57*1554)), p_AF_y(trial_index, round(0.018*1554):round(0.57*1554)));
    R_vec(trial_index) = R;
    
    x_vec(trial_index) = xc;
    y_vec(trial_index) = yc;
end

th = 0:pi/50:2*pi;
xunit = median(R_vec) * cos(th) + median(x_vec);
yunit = median(R_vec) * sin(th) + median(y_vec);

figure();
plot(xunit, yunit); hold on;

plot(p_AF_x(:, round(0.18*1554):round(0.33*1554))', p_AF_y(:, round(0.18*1554):round(0.33*1554))', 'r'); hold on;
plot(p_AF_x(:, round(0.33*1554):round(0.44*1554))', p_AF_y(:, round(0.33*1554):round(0.44*1554))', 'g'); hold on;
plot(p_AF_x(:, round(0.44*1554):round(0.57*1554))', p_AF_y(:, round(0.44*1554):round(0.57*1554))', 'b'); hold on;
hold off;

R_matrix = [];
for trial_index = 1:size(p0_raw_data_ind, 1)
    [L, R, K] = curvature([p_AF_x(trial_index, round(0.18*1554):round(0.57*1554))', p_AF_y(trial_index, round(0.18*1554):round(0.57*1554))']);
    R_matrix(trial_index, :) = R;
end

figure();
th = 0:pi/50:2*pi;
for trial_index = 1:size(p0_raw_data_ind, 1)
    x_i = R_vec(trial_index) * cos(th) + mean(x_vec(trial_index));
    y_i = R_vec(trial_index) * sin(th) + mean(y_vec(trial_index));

    plot(x_i, y_i); hold on;
end