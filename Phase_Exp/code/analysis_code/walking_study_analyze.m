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
u_lim=0.25;
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




%% Calculation of means for plotting...not really used in analysis
% weight1m=nanmean(weight1);
% weight2m=nanmean(weight2);
% weight3m=nanmean(weight3);
% weight4m=trimmean(weight4,30);
% 
% cop4m=trimmean(p0_cop_torque,30);
% 
p0_plat_torquem=trimmean(p0_act_torque,30); %NOT created from raw data. No need to remove outlier indices
p0_plat_posm=trimmean(p0_plat_pos(~p0_removed_ind, :),30); %Created from raw data. Need to remove outlier indices
p0_foot_posm=nanmean(p0_foot_pos(~p0_removed_ind, :)); %Created from raw data. Need to remove outlier indices
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
% the torque is collected from 100 samples prior
% to 300 samples post perturbation (heel strike)
%It is then subtracted from observed torque to get
%differential torque.


%p0_peakend/st come from raw data, hence the weird indexing.
%p0_act_torque/p0_cop_torque come from outlier removed data. Hence no weird
%indexing
for i=1:p0_total_accepted
    point_i=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_1)+p0_peakst(p0_raw_data_ind(i))
    p0_plat_torque15(i,:)=p0_act_torque(i,point_i+PRE_PERT_WINDOW:point_i+POST_PERT_WINDOW);
    p0_cop_torque15(i,:)=p0_cop_torque(i,point_i+PRE_PERT_WINDOW:point_i+POST_PERT_WINDOW);
    p0_foot_pos15(i,:)=p0_foot_pos(p0_raw_data_ind(i),point_i+PRE_PERT_WINDOW:point_i+POST_PERT_WINDOW);
    
    point_j=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_2)+p0_peakst(p0_raw_data_ind(i))
    p0_plat_torque30(i,:)=p0_act_torque(i,point_j+PRE_PERT_WINDOW:point_j+POST_PERT_WINDOW);
    p0_cop_torque30(i,:)=p0_cop_torque(i,point_j+PRE_PERT_WINDOW:point_j+POST_PERT_WINDOW);
    p0_foot_pos30(i,:)=p0_foot_pos(p0_raw_data_ind(i),point_j+PRE_PERT_WINDOW:point_j+POST_PERT_WINDOW);
    
    point_k=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_3)+p0_peakst(p0_raw_data_ind(i))
    p0_plat_torque45(i,:)=p0_act_torque(i,point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    p0_cop_torque45(i,:)=p0_cop_torque(i,point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    p0_foot_pos45(i,:)=p0_foot_pos(p0_raw_data_ind(i),point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    
    point_k=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_4)+p0_peakst(p0_raw_data_ind(i))
    p0_plat_torque60(i,:)=p0_act_torque(i,point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    p0_cop_torque60(i,:)=p0_cop_torque(i,point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    p0_foot_pos60(i,:)=p0_foot_pos(p0_raw_data_ind(i),point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
end

% Weight during perturbation
for i=1:p0_total_accepted
    point_i=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_1)+p0_peakst(p0_raw_data_ind(i))
    weightr15(i,:)=weight0(p0_raw_data_ind(i),point_i+PRE_PERT_WINDOW:point_i+POST_PERT_WINDOW);
    point_j=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_2)+p0_peakst(p0_raw_data_ind(i))
    weightr30(i,:)=weight0(p0_raw_data_ind(i),point_j+PRE_PERT_WINDOW:point_j+POST_PERT_WINDOW);
    point_k=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_3)+p0_peakst(p0_raw_data_ind(i))
    weightr45(i,:)=weight0(p0_raw_data_ind(i),point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    point_k=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_4)+p0_peakst(p0_raw_data_ind(i))
    weightr60(i,:)=weight0(p0_raw_data_ind(i),point_k+PRE_PERT_WINDOW:point_k+POST_PERT_WINDOW);
    
end


%CoP during perturbation
for i=1:p0_total_accepted
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

%CoP, BW, EMG Vals, etc. BioMech vals just before/at perturbation onset
for i = 1:size(p1_cop_torque, 1)
    p1_cop_val(i,:)=p1_cop_torque(i, p1_peakst(i))./weight1(i, p1_peakst(i));
    p2_cop_val(i,:)=p2_cop_torque(i, p2_peakst(i))./weight2(i, p2_peakst(i));
    
    p1_weight_val(i, :) = weight1(i, p1_peakst(i));
    p2_weight_val(i, :) = weight2(i, p2_peakst(i));
    
    p1_FZ(i, :) = force1_5(i, p1_peakst(i)) + force1_6(i, p1_peakst(i));
    p2_FZ(i, :) = force2_5(i, p2_peakst(i)) + force2_6(i, p2_peakst(i));
    
    p1_EMG_vals.TA(i, 1)  = p1_EMG.TA(i, p1_peakst(i));
    p1_EMG_vals.SOL(i, 1) = p1_EMG.SOL(i, p1_peakst(i));
    p1_EMG_vals.PL(i, 1)  = p1_EMG.PL(i, p1_peakst(i));
    p1_EMG_vals.GCA(i, 1) = p1_EMG.GCA(i, p1_peakst(i));
    
    p2_EMG_vals.TA(i, 1)  = p2_EMG.TA(i, p2_peakst(i));
    p2_EMG_vals.SOL(i, 1) = p2_EMG.SOL(i, p2_peakst(i));
    p2_EMG_vals.PL(i, 1)  = p2_EMG.PL(i, p2_peakst(i));
    p2_EMG_vals.GCA(i, 1) = p2_EMG.GCA(i, p2_peakst(i));
    
    p1_time_since_healstrike = (p1_peakst - TRIAL_WINDOW_PRE_PERT)*SAMPLE_PERIOD;
    p2_time_since_healstrike = (p2_peakst - TRIAL_WINDOW_PRE_PERT)*SAMPLE_PERIOD;
    
end
for i = 1:size(p3_cop_torque, 1)
    p3_cop_val(i,:)=p3_cop_torque(i, p3_peakst(i))./weight3(i, p3_peakst(i));
    p4_cop_val(i,:)=p4_cop_torque(i, p4_peakst(i))./weight4(i, p4_peakst(i));

    p3_weight_val(i, :) = weight3(i, p3_peakst(i));
    p4_weight_val(i, :) = weight4(i, p4_peakst(i));
    
    p3_FZ(i, :) = force3_5(i, p3_peakst(i)) + force3_6(i, p3_peakst(i));
    p4_FZ(i, :) = force4_5(i, p4_peakst(i)) + force4_6(i, p4_peakst(i));

    p3_EMG_vals.TA(i, 1)  = p3_EMG.TA(i, p3_peakst(i));
    p3_EMG_vals.SOL(i, 1) = p3_EMG.SOL(i, p3_peakst(i));
    p3_EMG_vals.PL(i, 1)  = p3_EMG.PL(i, p3_peakst(i));
    p3_EMG_vals.GCA(i, 1) = p3_EMG.GCA(i, p3_peakst(i));
    
    p4_EMG_vals.TA(i, 1)  = p4_EMG.TA(i, p4_peakst(i));
    p4_EMG_vals.SOL(i, 1) = p4_EMG.SOL(i, p4_peakst(i));
    p4_EMG_vals.PL(i, 1)  = p4_EMG.PL(i, p4_peakst(i));
    p4_EMG_vals.GCA(i, 1) = p4_EMG.GCA(i, p4_peakst(i));

    p3_time_since_healstrike = (p3_peakst - TRIAL_WINDOW_PRE_PERT)*SAMPLE_PERIOD;
    p4_time_since_healstrike = (p4_peakst - TRIAL_WINDOW_PRE_PERT)*SAMPLE_PERIOD;
end


for i=1:p0_total_accepted
    point_i=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_1)+p0_peakst(p0_raw_data_ind(i))
    ta15(i)=mean(ta_emg(p0_raw_data_ind(i),point_i-25:point_i+25));
    pl15(i)=mean(pl_emg(p0_raw_data_ind(i),point_i-25:point_i+25));
    sol15(i)=mean(sol_emg(p0_raw_data_ind(i),point_i-25:point_i+25));
    gca15(i)=mean(gca_emg(p0_raw_data_ind(i),point_i-25:point_i+25));
   
    point_j=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_2)+p0_peakst(p0_raw_data_ind(i))
    ta30(i)=mean(ta_emg(p0_raw_data_ind(i),point_j-25:point_j+25));
    pl30(i)=mean(pl_emg(p0_raw_data_ind(i),point_j-25:point_j+25));
    sol30(i)=mean(sol_emg(p0_raw_data_ind(i),point_j-25:point_j+25));
    gca30(i)=mean(gca_emg(p0_raw_data_ind(i),point_j-25:point_j+25));
    
    point_k=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(i))*PERT_POINT_3)+p0_peakst(i)
    ta45(i)=mean(ta_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    pl45(i)=mean(pl_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    sol45(i)=mean(sol_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    gca45(i)=mean(gca_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    
    point_k=floor((p0_peakend(p0_raw_data_ind(i))-p0_peakst(p0_raw_data_ind(i)))*PERT_POINT_4)+p0_peakst(p0_raw_data_ind(i))
    ta60(i)=mean(ta_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    pl60(i)=mean(pl_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    sol60(i)=mean(sol_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
    gca60(i)=mean(gca_emg(p0_raw_data_ind(i),point_k-25:point_k+25));
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


%% Obtaining differential data for impedance analysis
% data is broken into chunks of 400 points: PRE_PERT_WINDOW
% before pert and POST_PERT_WINDOW after pert to create "Perturbation Window". 


for i = 1:size(p1_peakst , 2)
    p1_foot_pos_segment(i,:)=p1_foot_pos(i,p1_peakst(i)+PRE_PERT_WINDOW:p1_peakst(i)+POST_PERT_WINDOW);
    p2_foot_pos_segment(i,:)=p2_foot_pos(i,p2_peakst(i)+PRE_PERT_WINDOW:p2_peakst(i)+POST_PERT_WINDOW);
     

    p1_plat_torque_segment(i,:)=p1_act_torque(i,p1_peakst(i)+PRE_PERT_WINDOW:p1_peakst(i)+POST_PERT_WINDOW);
    p2_plat_torque_segment(i,:)=p2_act_torque(i,p2_peakst(i)+PRE_PERT_WINDOW:p2_peakst(i)+POST_PERT_WINDOW);
    
end
for i = 1:size(p3_peakst, 2)
    p3_foot_pos_segment(i,:)=p3_foot_pos(i,p3_peakst(i)+PRE_PERT_WINDOW:p3_peakst(i)+POST_PERT_WINDOW);
    p4_foot_pos_segment(i,:)=p4_foot_pos(i,p4_peakst(i)+PRE_PERT_WINDOW:p4_peakst(i)+POST_PERT_WINDOW);
    
    p3_plat_torque_segment(i,:)=p3_act_torque(i,p3_peakst(i)+PRE_PERT_WINDOW:p3_peakst(i)+POST_PERT_WINDOW);
    p4_plat_torque_segment(i,:)=p4_act_torque(i,p4_peakst(i)+PRE_PERT_WINDOW:p4_peakst(i)+POST_PERT_WINDOW);
end
 
%% ---- Outlier Rejection ---- %%
%Position Oulier Rejection
[p1_pos_seg_outliers_rm, p1_pos_seg_removed_ind] = rmoutliers(p1_foot_pos_segment, 'ThresholdFactor', outlier_criterion_std);
[p2_pos_seg_outliers_rm, p2_pos_seg_removed_ind] = rmoutliers(p2_foot_pos_segment, 'ThresholdFactor', outlier_criterion_std);
[p3_pos_seg_outliers_rm, p3_pos_seg_removed_ind] = rmoutliers(p3_foot_pos_segment, 'ThresholdFactor', outlier_criterion_std);
[p4_pos_seg_outliers_rm, p4_pos_seg_removed_ind] = rmoutliers(p4_foot_pos_segment, 'ThresholdFactor', outlier_criterion_std);

%Torque Outlier Rejection
[p1_trq_seg_outliers_rm, p1_trq_seg_removed_ind] = rmoutliers(p1_plat_torque_segment, 'ThresholdFactor', outlier_criterion_std);
[p2_trq_seg_outliers_rm, p2_trq_seg_removed_ind] = rmoutliers(p2_plat_torque_segment, 'ThresholdFactor', outlier_criterion_std);
[p3_trq_seg_outliers_rm, p3_trq_seg_removed_ind] = rmoutliers(p3_plat_torque_segment, 'ThresholdFactor', outlier_criterion_std);
[p4_trq_seg_outliers_rm, p4_trq_seg_removed_ind] = rmoutliers(p4_plat_torque_segment, 'ThresholdFactor', outlier_criterion_std);

%Foot Placement Outlier Rejection
p1_foot_placement_rej_ind = (img1_pos < -0.5)' | (img1_pos > POS_REJECTION_LIMIT)';
p2_foot_placement_rej_ind = (img2_pos < -0.5)' | (img2_pos > POS_REJECTION_LIMIT)';
p3_foot_placement_rej_ind = (img3_pos < -0.5)' | (img3_pos > POS_REJECTION_LIMIT)';
p4_foot_placement_rej_ind = (img4_pos < -0.5)' | (img4_pos > POS_REJECTION_LIMIT)';


%Combine effects of outier rejection
p1_seg_removed_ind = p1_pos_seg_removed_ind | p1_trq_seg_removed_ind | p1_foot_placement_rej_ind;
p2_seg_removed_ind = p2_pos_seg_removed_ind | p2_trq_seg_removed_ind | p2_foot_placement_rej_ind;
p3_seg_removed_ind = p3_pos_seg_removed_ind | p3_trq_seg_removed_ind | p3_foot_placement_rej_ind;
p4_seg_removed_ind = p4_pos_seg_removed_ind | p4_trq_seg_removed_ind | p4_foot_placement_rej_ind;


%% Find Platform Position Differentials
for i=1:size(p1_plat_pos, 1)

    %Differential Platform Position
    diff_p1_plat_pos(i,:)=p1_plat_pos(i,p1_peakst(i)+PRE_PERT_WINDOW:p1_peakst(i)+POST_PERT_WINDOW);
    diff_p2_plat_pos(i,:)=p2_plat_pos(i,p2_peakst(i)+PRE_PERT_WINDOW:p2_peakst(i)+POST_PERT_WINDOW);
    
    diff_p1_plat_pos(i,:)=diff_p1_plat_pos(i,:)-diff_p1_plat_pos(i,100);
    diff_p2_plat_pos(i,:)=diff_p2_plat_pos(i,:)-diff_p2_plat_pos(i,100);
        
end

for i=1:size(p3_plat_pos, 1)
    diff_p3_plat_pos(i,:)=p3_plat_pos(i,p3_peakst(i)+PRE_PERT_WINDOW:p3_peakst(i)+POST_PERT_WINDOW);
    diff_p3_plat_pos(i,:)=diff_p3_plat_pos(i,:)-diff_p3_plat_pos(i,100);
    
    diff_p4_plat_pos(i,:)=p4_plat_pos(i,p4_peakst(i)+PRE_PERT_WINDOW:p4_peakst(i)+POST_PERT_WINDOW);
    diff_p4_plat_pos(i,:)=diff_p4_plat_pos(i,:)-diff_p4_plat_pos(i,100);

end

%---- Put EMG Curves for Each Stance Percentage into a Struct ---- %% 
p0_EMG_struct1.TA = ta30;
p0_EMG_struct1.PL = pl30;
p0_EMG_struct1.SOL = sol30;
p0_EMG_struct1.GCA = gca30;

p0_EMG_struct2.TA = ta45;
p0_EMG_struct2.PL = pl45;
p0_EMG_struct2.SOL = sol45;
p0_EMG_struct2.GCA = gca45;

p0_EMG_struct3.TA = ta60;
p0_EMG_struct3.PL = pl60;
p0_EMG_struct3.SOL = sol60;
p0_EMG_struct3.GCA = gca60;

p0_EMG_struct4.TA = ta15;
p0_EMG_struct4.PL = pl15;
p0_EMG_struct4.SOL = sol15;
p0_EMG_struct4.GCA = gca15;

p1_EMG_vals.TA  = p1_EMG_vals.TA(~p1_seg_removed_ind, :);
p1_EMG_vals.PL  = p1_EMG_vals.PL(~p1_seg_removed_ind, :);
p1_EMG_vals.SOL = p1_EMG_vals.SOL(~p1_seg_removed_ind, :);
p1_EMG_vals.GCA = p1_EMG_vals.GCA(~p1_seg_removed_ind, :);

p2_EMG_vals.TA  = p2_EMG_vals.TA(~p2_seg_removed_ind, :);
p2_EMG_vals.PL  = p2_EMG_vals.PL(~p2_seg_removed_ind, :);
p2_EMG_vals.SOL = p2_EMG_vals.SOL(~p2_seg_removed_ind, :);
p2_EMG_vals.GCA = p2_EMG_vals.GCA(~p2_seg_removed_ind, :);

p3_EMG_vals.TA  = p3_EMG_vals.TA(~p3_seg_removed_ind, :);
p3_EMG_vals.PL  = p3_EMG_vals.PL(~p3_seg_removed_ind, :);
p3_EMG_vals.SOL = p3_EMG_vals.SOL(~p3_seg_removed_ind, :);
p3_EMG_vals.GCA = p3_EMG_vals.GCA(~p3_seg_removed_ind, :);

p4_EMG_vals.TA  = p4_EMG_vals.TA(~p4_seg_removed_ind, :);
p4_EMG_vals.PL  = p4_EMG_vals.PL(~p4_seg_removed_ind, :);
p4_EMG_vals.SOL = p4_EMG_vals.SOL(~p4_seg_removed_ind, :);
p4_EMG_vals.GCA = p4_EMG_vals.GCA(~p4_seg_removed_ind, :);

[diff_p1_foot_pos, diff_p1_plat_torque, diff_p1_plat_pos, bio_factors_p1] = get_pos_torque_diff(p0_foot_pos30, p1_foot_pos_segment(~p1_seg_removed_ind, :), p0_plat_torque30, p1_plat_torque_segment(~p1_seg_removed_ind, :), diff_p1_plat_pos(~p1_seg_removed_ind, :), p1_cop_val(~p1_seg_removed_ind), p1_weight_val(~p1_seg_removed_ind), p1_EMG_vals, p1_FZ(~p1_seg_removed_ind), p1_time_since_healstrike(~p1_seg_removed_ind));
[diff_p2_foot_pos, diff_p2_plat_torque, diff_p2_plat_pos, bio_factors_p2] = get_pos_torque_diff(p0_foot_pos45, p2_foot_pos_segment(~p2_seg_removed_ind, :), p0_plat_torque45, p2_plat_torque_segment(~p2_seg_removed_ind, :), diff_p2_plat_pos(~p2_seg_removed_ind, :), p2_cop_val(~p2_seg_removed_ind), p2_weight_val(~p2_seg_removed_ind), p2_EMG_vals, p2_FZ(~p2_seg_removed_ind), p2_time_since_healstrike(~p2_seg_removed_ind));
[diff_p3_foot_pos, diff_p3_plat_torque, diff_p3_plat_pos, bio_factors_p3] = get_pos_torque_diff(p0_foot_pos60, p3_foot_pos_segment(~p3_seg_removed_ind, :), p0_plat_torque60, p3_plat_torque_segment(~p3_seg_removed_ind, :), diff_p3_plat_pos(~p3_seg_removed_ind, :), p3_cop_val(~p3_seg_removed_ind), p3_weight_val(~p3_seg_removed_ind), p3_EMG_vals, p3_FZ(~p3_seg_removed_ind), p3_time_since_healstrike(~p3_seg_removed_ind));
[diff_p4_foot_pos, diff_p4_plat_torque, diff_p4_plat_pos, bio_factors_p4] = get_pos_torque_diff(p0_foot_pos15, p4_foot_pos_segment(~p4_seg_removed_ind, :), p0_plat_torque15, p4_plat_torque_segment(~p4_seg_removed_ind, :), diff_p4_plat_pos(~p4_seg_removed_ind, :), p4_cop_val(~p4_seg_removed_ind), p4_weight_val(~p4_seg_removed_ind), p4_EMG_vals, p4_FZ(~p4_seg_removed_ind), p4_time_since_healstrike(~p4_seg_removed_ind));




%% Get's Derivatives for ankle and platform
%Uses MATLAB's gradient function to find the derivatives, which uses the
%central difference method -> Should have better properties (Though, considering
%the sampling time vs. the foot dynamics, the difference should be very
%small). 
[diff_p1_foot_vel, diff_p1_foot_acc] = get_derivatives(diff_p1_foot_pos, SAMPLE_PERIOD);
[diff_p2_foot_vel, diff_p2_foot_acc] = get_derivatives(diff_p2_foot_pos, SAMPLE_PERIOD);
[diff_p3_foot_vel, diff_p3_foot_acc] = get_derivatives(diff_p3_foot_pos, SAMPLE_PERIOD);
[diff_p4_foot_vel, diff_p4_foot_acc] = get_derivatives(diff_p4_foot_pos, SAMPLE_PERIOD);

[diff_p1_plat_vel, diff_p1_plat_acc] = get_derivatives(diff_p1_plat_pos, SAMPLE_PERIOD);
[diff_p2_plat_vel, diff_p2_plat_acc] = get_derivatives(diff_p2_plat_pos, SAMPLE_PERIOD);
[diff_p3_plat_vel, diff_p3_plat_acc] = get_derivatives(diff_p3_plat_pos, SAMPLE_PERIOD);
[diff_p4_plat_vel, diff_p4_plat_acc] = get_derivatives(diff_p4_plat_pos, SAMPLE_PERIOD);



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
for i=1:BOOTSTRAP_COUNT
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torque(i,:)-0.1945*diff_p1_plat_acc(i,:);%+0.02*diff_p1_foot_accm;%+1*diff_p1_foot_velm;
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torqueimp(i,:)-diff_p1_plat_torqueimp(i,100);
    
   
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torque(i,:)-0.1945*diff_p2_plat_acc(i,:);%+0.013*diff_p2_foot_accm;
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torqueimp(i,:)-diff_p2_plat_torqueimp(i,100);

    
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torque(i,:)-0.1945*diff_p3_plat_acc(i,:);%-0.09*diff_p3_foot_accm;
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torqueimp(i,:)-diff_p3_plat_torqueimp(i,100);

    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torque(i,:)-0.1945*diff_p4_plat_acc(i,:);%+0.02*diff_p4_foot_accm;%+1*diff_p4_foot_velm;
    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torqueimp(i,:)-diff_p4_plat_torqueimp(i,100);
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

diff_p4_foot_accm=diff_p4_foot_accm;
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);

diff_p4_foot_posm=trimmean(diff_p4_foot_pos,30);
diff_p4_foot_velm=trimmean(diff_p4_foot_vel,30);
diff_p4_plat_velm=trimmean(diff_p4_plat_vel,30);
diff_p4_foot_accm=trimmean(diff_p4_foot_acc,30);
% constained linear regression
%Value is made NaN if foot position is not on
%platform

%Just temporary
if strcmp(outlier_removal_type, "NONE")
   p4_removed_ind =  logical(zeros(size(diff_p4_foot_pos, 1), 1));
   p1_removed_ind =  logical(zeros(size(diff_p1_foot_pos, 1), 1));
   p2_removed_ind =  logical(zeros(size(diff_p2_foot_pos, 1), 1));
   p3_removed_ind =  logical(zeros(size(diff_p3_foot_pos, 1), 1));
end


%%---P1---%%
% 
% %Outlier Rejection
% if strcmp(outlier_removal_type, "POS")
%     [p1_outliers_ind, p1_removed_ind] = rmoutliers(diff_p1_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "ACC_TOR")
%     [p1_outliers_ind, p1_removed_ind] = rmoutliers(diff_p1_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "POS_TRQ")
%       [p1_pos_outliers_ind, p1_pos_removed_ind] = rmoutliers(diff_p1_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       [p1_trq_outliers_ind, p1_trq_removed_ind] = rmoutliers(diff_p1_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       p1_removed_ind = p1_pos_removed_ind | p1_trq_removed_ind;
% elseif  strcmp(outlier_removal_type, "POS_VEL_ACC_TRQ")
%     [p1_pos_outliers_ind, p1_pos_removed_ind] = rmoutliers(diff_p1_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p1_vel_outliers_ind, p1_vel_removed_ind] = rmoutliers(diff_p1_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p1_acc_outliers_ind, p1_acc_removed_ind] = rmoutliers(diff_p1_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p1_trq_outliers_ind, p1_trq_removed_ind] = rmoutliers(diff_p1_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     p1_removed_ind = p1_pos_removed_ind | p1_vel_removed_ind | p1_acc_removed_ind | p1_trq_removed_ind;
% end

% p1_removed_ind = p1_removed_ind | (abs(img1_pos') > POS_REJECTION_LIMIT);


%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
% diff_p1_pos_data_vec = diff_p1_foot_pos(~p1_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p1_vel_data_vec = diff_p1_foot_vel(~p1_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p1_acc_data_vec = diff_p1_foot_acc(~p1_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p1_data_matrix = [diff_p1_pos_data_vec(:), diff_p1_vel_data_vec(:), diff_p1_acc_data_vec(:)];
% 
% p1_torque_data_matrix = diff_p1_plat_torqueimp(~p1_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p1_output_data_vec = p1_torque_data_matrix(:);
% 
 A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
 B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

[regress_coeffs_p1, ci_p1, GoF_p1, vaf_info_p1] = regress_after_bootstrap(diff_p1_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p1_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p1_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p1_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   A1, B1);

%%---P2---%%

% %Outlier Rejection
% if strcmp(outlier_removal_type, "POS")
%     [p2_outliers_ind, p2_removed_ind] = rmoutliers(diff_p2_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "ACC_TOR")
%     [p2_outliers_ind, p2_removed_ind] = rmoutliers(diff_p2_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "POS_TRQ")
%       [p2_pos_outliers_ind, p2_pos_removed_ind] = rmoutliers(diff_p2_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       [p2_trq_outliers_ind, p2_trq_removed_ind] = rmoutliers(diff_p2_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       p2_removed_ind = p2_pos_removed_ind | p2_trq_removed_ind;
% elseif  strcmp(outlier_removal_type, "POS_VEL_ACC_TRQ")
%     [p2_pos_outliers_ind, p2_pos_removed_ind] = rmoutliers(diff_p2_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p2_vel_outliers_ind, p2_vel_removed_ind] = rmoutliers(diff_p2_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p2_acc_outliers_ind, p2_acc_removed_ind] = rmoutliers(diff_p2_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p2_trq_outliers_ind, p2_trq_removed_ind] = rmoutliers(diff_p2_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     p2_removed_ind = p2_pos_removed_ind | p2_vel_removed_ind | p2_acc_removed_ind | p2_trq_removed_ind;
% end

% p2_removed_ind = p2_removed_ind | (abs(img2_pos') > POS_REJECTION_LIMIT);

%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
% diff_p2_pos_data_vec = diff_p2_foot_pos(~p2_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p2_vel_data_vec = diff_p2_foot_vel(~p2_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p2_acc_data_vec = diff_p2_foot_acc(~p2_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p2_data_matrix = [diff_p2_pos_data_vec(:), diff_p2_vel_data_vec(:), diff_p2_acc_data_vec(:)];
% 
% p2_torque_data_matrix = diff_p2_plat_torqueimp(~p2_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p2_output_data_vec = p2_torque_data_matrix(:);
% 
% A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
% B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

[regress_coeffs_p2, ci_p2, GoF_p2, vaf_info_p2] = regress_after_bootstrap(diff_p2_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p2_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p2_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p2_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   A1, B1);
%%---P3---%%
% 
% %Outlier Rejection
% if strcmp(outlier_removal_type, "POS")
%     [p3_outliers_ind, p3_removed_ind] = rmoutliers(diff_p3_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "ACC_TOR")
%     [p3_outliers_ind, p3_removed_ind] = rmoutliers(diff_p3_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "POS_TRQ")
%       [p3_pos_outliers_ind, p3_pos_removed_ind] = rmoutliers(diff_p3_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       [p3_trq_outliers_ind, p3_trq_removed_ind] = rmoutliers(diff_p3_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       p3_removed_ind = p3_pos_removed_ind | p3_trq_removed_ind;
% elseif  strcmp(outlier_removal_type, "POS_VEL_ACC_TRQ")
%     [p3_pos_outliers_ind, p3_pos_removed_ind] = rmoutliers(diff_p3_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p3_vel_outliers_ind, p3_vel_removed_ind] = rmoutliers(diff_p3_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p3_acc_outliers_ind, p3_acc_removed_ind] = rmoutliers(diff_p3_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p3_trq_outliers_ind, p3_trq_removed_ind] = rmoutliers(diff_p3_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     p3_removed_ind = p3_pos_removed_ind | p3_vel_removed_ind | p3_acc_removed_ind | p3_trq_removed_ind;
% end

% p3_removed_ind = p3_removed_ind | (abs(img3_pos') > POS_REJECTION_LIMIT);

%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
% diff_p3_pos_data_vec = diff_p3_foot_pos(~p3_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p3_vel_data_vec = diff_p3_foot_vel(~p3_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p3_acc_data_vec = diff_p3_foot_acc(~p3_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p3_data_matrix = [diff_p3_pos_data_vec(:), diff_p3_vel_data_vec(:), diff_p3_acc_data_vec(:)];
% 
% p3_torque_data_matrix = diff_p3_plat_torqueimp(~p3_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p3_output_data_vec = p3_torque_data_matrix(:);
% 
% A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
% B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

[regress_coeffs_p3, ci_p3, GoF_p3, vaf_info_p3] = regress_after_bootstrap(diff_p3_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p3_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p3_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p3_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   A1, B1);
%%---P4---%%

%Outlier Rejection
% if strcmp(outlier_removal_type, "POS")
%     [p4_outliers_ind, p4_removed_ind] = rmoutliers(diff_p4_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "ACC_TOR")
%     [p4_outliers_ind, p4_removed_ind] = rmoutliers(diff_p4_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
% elseif strcmp(outlier_removal_type, "POS_TRQ")
%       [p4_pos_outliers_ind, p4_pos_removed_ind] = rmoutliers(diff_p4_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       [p4_trq_outliers_ind, p4_trq_removed_ind] = rmoutliers(diff_p4_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%       p4_removed_ind = p4_pos_removed_ind | p4_trq_removed_ind;
% elseif  strcmp(outlier_removal_type, "POS_VEL_ACC_TRQ")
%     [p4_pos_outliers_ind, p4_pos_removed_ind] = rmoutliers(diff_p4_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p4_vel_outliers_ind, p4_vel_removed_ind] = rmoutliers(diff_p4_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p4_acc_outliers_ind, p4_acc_removed_ind] = rmoutliers(diff_p4_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     [p4_trq_outliers_ind, p4_trq_removed_ind] = rmoutliers(diff_p4_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), 'ThresholdFactor', outlier_criterion_std);
%     p4_removed_ind = p4_pos_removed_ind | p4_vel_removed_ind | p4_acc_removed_ind | p4_trq_removed_ind;
% end

% p4_removed_ind = p4_removed_ind | (abs(img4_pos') > POS_REJECTION_LIMIT);

%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
% diff_p4_pos_data_vec = diff_p4_foot_pos(~p4_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p4_vel_data_vec = diff_p4_foot_vel(~p4_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% diff_p4_acc_data_vec = diff_p4_foot_acc(~p4_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p4_data_matrix = [diff_p4_pos_data_vec(:), diff_p4_vel_data_vec(:), diff_p4_acc_data_vec(:)];
% 
% p4_torque_data_matrix = diff_p4_plat_torqueimp(~p4_removed_ind, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX)';
% p4_output_data_vec = p4_torque_data_matrix(:);
% 
% A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
% B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

[regress_coeffs_p4, ci_p4, GoF_p4, vaf_info_p4] = regress_after_bootstrap(diff_p4_foot_pos(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p4_foot_vel(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p4_foot_acc(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   diff_p4_plat_torqueimp(:, REGRESSION_WINDOW_MIN_INDEX:REGRESSION_WINDOW_MAX_INDEX), ...
                                                   A1, B1);


remaining_after_VAF1 = vaf_info_p1.remaining_after_VAF;
remaining_after_VAF2 = vaf_info_p2.remaining_after_VAF;                                               
remaining_after_VAF3 = vaf_info_p3.remaining_after_VAF;
remaining_after_VAF4 = vaf_info_p4.remaining_after_VAF;


%Remove the low %VAF bootstrap sample regression results from the corresponding bio_factor data
bio_factors_p1.CoP     = bio_factors_p1.CoP(vaf_info_p1.removal_indices);
bio_factors_p1.ankle_ang = bio_factors_p1.ankle_ang(vaf_info_p1.removal_indices);
bio_factors_p1.Weight  = bio_factors_p1.Weight(vaf_info_p1.removal_indices);
bio_factors_p1.EMG.TA  = bio_factors_p1.EMG.TA(vaf_info_p1.removal_indices);
bio_factors_p1.EMG.PL  = bio_factors_p1.EMG.PL(vaf_info_p1.removal_indices);
bio_factors_p1.EMG.SOL = bio_factors_p1.EMG.SOL(vaf_info_p1.removal_indices);
bio_factors_p1.EMG.GCA = bio_factors_p1.EMG.GCA(vaf_info_p1.removal_indices); 

bio_factors_p2.CoP     = bio_factors_p2.CoP(vaf_info_p2.removal_indices);
bio_factors_p2.ankle_ang = bio_factors_p2.ankle_ang(vaf_info_p2.removal_indices);
bio_factors_p2.Weight  = bio_factors_p2.Weight(vaf_info_p2.removal_indices);
bio_factors_p2.EMG.TA  = bio_factors_p2.EMG.TA(vaf_info_p2.removal_indices);
bio_factors_p2.EMG.PL  = bio_factors_p2.EMG.PL(vaf_info_p2.removal_indices);
bio_factors_p2.EMG.SOL = bio_factors_p2.EMG.SOL(vaf_info_p2.removal_indices);
bio_factors_p2.EMG.GCA = bio_factors_p2.EMG.GCA(vaf_info_p2.removal_indices); 

bio_factors_p3.CoP     = bio_factors_p3.CoP(vaf_info_p3.removal_indices);
bio_factors_p3.ankle_ang = bio_factors_p3.ankle_ang(vaf_info_p3.removal_indices);
bio_factors_p3.Weight  = bio_factors_p3.Weight(vaf_info_p3.removal_indices);
bio_factors_p3.EMG.TA  = bio_factors_p3.EMG.TA(vaf_info_p3.removal_indices);
bio_factors_p3.EMG.PL  = bio_factors_p3.EMG.PL(vaf_info_p3.removal_indices);
bio_factors_p3.EMG.SOL = bio_factors_p3.EMG.SOL(vaf_info_p3.removal_indices);
bio_factors_p3.EMG.GCA = bio_factors_p3.EMG.GCA(vaf_info_p3.removal_indices); 

bio_factors_p4.CoP     = bio_factors_p4.CoP(vaf_info_p4.removal_indices);
bio_factors_p4.ankle_ang = bio_factors_p4.ankle_ang(vaf_info_p4.removal_indices);
bio_factors_p4.Weight  = bio_factors_p4.Weight(vaf_info_p4.removal_indices);
bio_factors_p4.EMG.TA  = bio_factors_p4.EMG.TA(vaf_info_p4.removal_indices);
bio_factors_p4.EMG.PL  = bio_factors_p4.EMG.PL(vaf_info_p4.removal_indices);
bio_factors_p4.EMG.SOL = bio_factors_p4.EMG.SOL(vaf_info_p4.removal_indices);
bio_factors_p4.EMG.GCA = bio_factors_p4.EMG.GCA(vaf_info_p4.removal_indices); 



%% Storing final data


%Basic regression plots
% if(plot_figs==1)
%     plot_regression_results(diff_p4_foot_pos(~p4_removed_ind, :), diff_p4_foot_vel(~p4_removed_ind, :), ...
%                             diff_p4_foot_acc(~p4_removed_ind, :), diff_p4_plat_torqueimp(~p4_removed_ind, :), p4_regressors, '18')
%     plot_regression_results(diff_p1_foot_pos(~p1_removed_ind, :), diff_p1_foot_vel(~p1_removed_ind, :), ...
%                             diff_p1_foot_acc(~p1_removed_ind, :), diff_p1_plat_torqueimp(~p1_removed_ind, :), p1_regressors, '31')
%     plot_regression_results(diff_p2_foot_pos(~p2_removed_ind, :), diff_p2_foot_vel(~p2_removed_ind, :), ...
%                             diff_p2_foot_acc(~p2_removed_ind, :), diff_p2_plat_torqueimp(~p2_removed_ind, :), p2_regressors, '44')
%     plot_regression_results(diff_p3_foot_pos(~p3_removed_ind, :), diff_p3_foot_vel(~p3_removed_ind, :), ...
%                             diff_p3_foot_acc(~p3_removed_ind, :), diff_p3_plat_torqueimp(~p3_removed_ind, :), p3_regressors, '57')
% 
% end

% %% - Individual Trial Regression Low %VAF Removel - %
% [GoF_p4_vec, regression_coeffs_p4] = get_trial_VAF(diff_p4_pos_data_vec', diff_p4_vel_data_vec', diff_p4_acc_data_vec', ...
%                          p4_torque_data_matrix', 0.7, A1, B1);
% [GoF_p1_vec, regression_coeffs_p1] = get_trial_VAF(diff_p1_pos_data_vec', diff_p1_vel_data_vec', diff_p1_acc_data_vec', ...
%                          p1_torque_data_matrix', 0.7, A1, B1);
%           
% [GoF_p2_vec, regression_coeffs_p2] = get_trial_VAF(diff_p2_pos_data_vec', diff_p2_vel_data_vec', diff_p2_acc_data_vec', ...
%                          p2_torque_data_matrix', 0.7, A1, B1);
% [GoF_p3_vec, regression_coeffs_p3] = get_trial_VAF(diff_p3_pos_data_vec', diff_p3_vel_data_vec', diff_p3_acc_data_vec', ...
%                          p3_torque_data_matrix', 0.7, A1, B1);
% 
%                    
% p4_VAF_removed_indices =  (GoF_p4_vec < VAF_REMOVAL_CRITERION);                     
% diff_p4_pos_data_vec   =  diff_p4_pos_data_vec(:,  ~p4_VAF_removed_indices);
% diff_p4_vel_data_vec   =  diff_p4_vel_data_vec(:,  ~p4_VAF_removed_indices);
% diff_p4_acc_data_vec   =  diff_p4_acc_data_vec(:,  ~p4_VAF_removed_indices);
% p4_torque_data_matrix  =  p4_torque_data_matrix(:, ~p4_VAF_removed_indices);
% 
% p1_VAF_removed_indices =  (GoF_p1_vec < VAF_REMOVAL_CRITERION);                     
% diff_p1_pos_data_vec   =  diff_p1_pos_data_vec(:,  ~p1_VAF_removed_indices);
% diff_p1_vel_data_vec   =  diff_p1_vel_data_vec(:,  ~p1_VAF_removed_indices);
% diff_p1_acc_data_vec   =  diff_p1_acc_data_vec(:,  ~p1_VAF_removed_indices);
% p1_torque_data_matrix  =  p1_torque_data_matrix(:, ~p1_VAF_removed_indices);
% 
% p2_VAF_removed_indices =  (GoF_p2_vec < VAF_REMOVAL_CRITERION);                     
% diff_p2_pos_data_vec   =  diff_p2_pos_data_vec(:,  ~p2_VAF_removed_indices);
% diff_p2_vel_data_vec   =  diff_p2_vel_data_vec(:,  ~p2_VAF_removed_indices);
% diff_p2_acc_data_vec   =  diff_p2_acc_data_vec(:,  ~p2_VAF_removed_indices);
% p2_torque_data_matrix  =  p2_torque_data_matrix(:, ~p2_VAF_removed_indices);
% 
% p3_VAF_removed_indices =  (GoF_p3_vec < VAF_REMOVAL_CRITERION);                     
% diff_p3_pos_data_vec   =  diff_p3_pos_data_vec(:,  ~p3_VAF_removed_indices);
% diff_p3_vel_data_vec   =  diff_p3_vel_data_vec(:,  ~p3_VAF_removed_indices);
% diff_p3_acc_data_vec   =  diff_p3_acc_data_vec(:,  ~p3_VAF_removed_indices);
% p3_torque_data_matrix  =  p3_torque_data_matrix(:, ~p3_VAF_removed_indices);
% 

%% - Bootstrapping Regression - %
%The order is from lower to greater stance phase percentage.
%The odd numbering is due to how the connections are made in the Simulink
%file.

% [mean_p4, ci_p4, GoF_p4] = bootstrapping_regression(diff_p4_pos_data_vec', diff_p4_vel_data_vec', diff_p4_acc_data_vec', ...
%                          p4_torque_data_matrix', 0.7, A1, B1);
% [mean_p1, ci_p1, GoF_p1] = bootstrapping_regression(diff_p1_pos_data_vec', diff_p1_vel_data_vec', diff_p1_acc_data_vec', ...
%                          p1_torque_data_matrix', 0.7, A1, B1);
%           
% [mean_p2, ci_p2, GoF_p2] = bootstrapping_regression(diff_p2_pos_data_vec', diff_p2_vel_data_vec', diff_p2_acc_data_vec', ...
%                          p2_torque_data_matrix', 0.7, A1, B1);
% [mean_p3, ci_p3, GoF_p3] = bootstrapping_regression(diff_p3_pos_data_vec', diff_p3_vel_data_vec', diff_p3_acc_data_vec', ...
%                          p3_torque_data_matrix', 0.7, A1, B1);

% [mean_p4, ci_p4, GoF_p4] = bootstrapping_regression_testing(diff_p4_pos_data_vec', diff_p4_vel_data_vec', diff_p4_acc_data_vec', ...
%                          p4_torque_data_matrix', 0.7, A1, B1);
% [mean_p1, ci_p1, GoF_p1] = bootstrapping_regression_testing(diff_p1_pos_data_vec', diff_p1_vel_data_vec', diff_p1_acc_data_vec', ...
%                          p1_torque_data_matrix', 0.7, A1, B1);
%           
% [mean_p2, ci_p2, GoF_p2] = bootstrapping_regression_testing(diff_p2_pos_data_vec', diff_p2_vel_data_vec', diff_p2_acc_data_vec', ...
%                          p2_torque_data_matrix', 0.7, A1, B1);
% [mean_p3, ci_p3, GoF_p3] = bootstrapping_regression_testing(diff_p3_pos_data_vec', diff_p3_vel_data_vec', diff_p3_acc_data_vec', ...
%                          p3_torque_data_matrix', 0.7, A1, B1);
%                      
mean_p1 = mean(regress_coeffs_p1);
mean_p2 = mean(regress_coeffs_p2);
mean_p3 = mean(regress_coeffs_p3);
mean_p4 = mean(regress_coeffs_p4);

imp_vals_bs = [mean_p4; mean_p1; mean_p2; mean_p3];
imp_vals_bs_s = [ci_p4; ci_p1; ci_p2; ci_p3];
     
if(plot_figs==1)
    plot_regression_results(diff_p4_foot_pos(:, :), diff_p4_foot_vel(:, :), ...
                            diff_p4_foot_acc(:, :), diff_p4_plat_torqueimp(:, :), mean_p4, strcat(RESULTS_DIR, '18'))
    plot_regression_results(diff_p1_foot_pos(:, :), diff_p1_foot_vel(:, :), ...
                            diff_p1_foot_acc(:, :), diff_p1_plat_torqueimp(:, :), mean_p1, strcat(RESULTS_DIR, '31'))
    plot_regression_results(diff_p2_foot_pos(:, :), diff_p2_foot_vel(:, :), ...
                            diff_p2_foot_acc(:, :), diff_p2_plat_torqueimp(:, :), mean_p2, strcat(RESULTS_DIR, '44'))
    plot_regression_results(diff_p3_foot_pos(:, :), diff_p3_foot_vel(:, :), ...
                            diff_p3_foot_acc(:, :), diff_p3_plat_torqueimp(:, :), mean_p3, strcat(RESULTS_DIR, '57'))

end                     
remaining_after_VAF = [remaining_after_VAF4; remaining_after_VAF1; remaining_after_VAF2; remaining_after_VAF3]
GoF_matrix = [GoF_p4; GoF_p1; GoF_p2; GoF_p3]


%Finds Mean/CI from the bootstrapped data using percentiles.
[weight_m, weight_s] = process_weight(bio_factors_p4.Weight', bio_factors_p1.Weight', bio_factors_p2.Weight', bio_factors_p3.Weight');
[cop_m, cop_s] = process_cop(bio_factors_p4.CoP', bio_factors_p1.CoP', bio_factors_p2.CoP', bio_factors_p3.CoP');
[emg_final] = process_EMG(bio_factors_p4.EMG, bio_factors_p1.EMG, bio_factors_p2.EMG, bio_factors_p3.EMG); 


pre_diff_selection_count = [sum(~p4_seg_removed_ind); sum(~p1_seg_removed_ind); sum(~p2_seg_removed_ind); sum(~p3_seg_removed_ind);];
% step_placement_number_remaining = [sum(~p4_removed_ind); sum(~p1_removed_ind); sum(~p2_removed_ind); sum(~p3_removed_ind)];
% VAF_number_remaining = [sum(~p4_VAF_removed_indices); sum(~p1_VAF_removed_indices); sum(~p2_VAF_removed_indices); sum(~p3_VAF_removed_indices)];


write_data_summary;

figure();
plot((diff_p4_plat_torqueimp(:, :))', 'k'); hold on;
plot((diff_p1_plat_torqueimp(:, :))', 'r'); hold on;
plot((diff_p2_plat_torqueimp(:, :))', 'g'); hold on;
plot((diff_p3_plat_torqueimp(:, :))', 'b'); hold off;
figure();
plot(mean(diff_p4_plat_torqueimp(:, :))', 'k'); hold on;
plot(mean(diff_p1_plat_torqueimp(:, :))', 'r'); hold on;
plot(mean(diff_p2_plat_torqueimp(:, :))', 'g'); hold on;
plot(mean(diff_p3_plat_torqueimp(:, :))', 'b'); hold off;

figure();
plot((diff_p4_foot_pos(:, :))', 'k'); hold on;
plot((diff_p1_foot_pos(:, :))', 'r'); hold on;
plot((diff_p2_foot_pos(:, :))', 'g'); hold on;
plot((diff_p3_foot_pos(:, :))', 'b'); hold off;
figure();
plot(mean(diff_p4_foot_pos(:, :))', 'k'); hold on;
plot(mean(diff_p1_foot_pos(:, :))', 'r'); hold on;
plot(mean(diff_p2_foot_pos(:, :))', 'g'); hold on;
plot(mean(diff_p3_foot_pos(:, :))', 'b'); hold off;

diff2 = mean(p2_plat_torque_segment(:, :)) - mean(p0_plat_torque45);
diff3 = mean(p3_plat_torque_segment(:, :)) - mean(p0_plat_torque60);
