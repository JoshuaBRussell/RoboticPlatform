function [TA_NORM, SOL_NORM, PL_NORM, GCA_NORM] = find_emg_normalization(DATA_DIR, DATA_FILE_NAME)

TA_EMG_SIG  = 1;
SOL_EMG_SIG = 2;
PL_EMG_SIG  = 3;
GCA_EMG_SIG = 4;

TRIAL_WINDOW_PRE_PERT = -400;
TRIAL_WINDOW_POST_PERT = 3500;


WEIGHT_PERCENTAGE_CUTOFF = 0.025;

NUM_TRAINING_TRIALS = 5;

[num,den] = butter(2,10/1000);
[filt_num, filt_den] = butter(2,10/1000);
% Add the trials you want to exclude in here
exclude=[0];%exclude=[1,2];
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);

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

for trials=1:1
    
  if(ismember(trials,exclude)==0)  
    if(trials<10)
        h = fopen(strcat(DATA_DIR, DATA_FILE_NAME));
    end
    
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    siz=size(Input1.data);
    
    %             %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','CWamplerate',2000);
    %d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 20, 'CWampleRate', 2000, 'DesignMethod', 'window');
    pert_torque=filtfilt(d1,Input1.data(:,7));
    f1=Input1.data(:,9)*53.4;
    
    f2=Input1.data(:,10)*53.4;

    f3=Input1.data(:,11)*53.4;

    f4=Input1.data(:,12)*53.4;
 
    f5=Input1.data(:,20)*53.4/2;
  
    f6=Input1.data(:,21)*53.4/2;
    
    ta=Input1.data(:,TA_EMG_SIG);
    sol=Input1.data(:,SOL_EMG_SIG);
    pl=Input1.data(:,PL_EMG_SIG);
    gca=Input1.data(:,GCA_EMG_SIG); 
    
    off_TA = mean(ta);
    off_SOL = mean(sol);
    off_PL = mean(pl);
    off_GCA = mean(gca);
  
    w1=filtfilt(d1,Input1.data(:,18));
  
    flag=Input1.data(:,17);
%     foot_pos_data=Input1.data(:,13);
%     foot_const=mean(foot_pos_data);
%     foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
%     plat_pos_data=Input1.data(:,14);
%     plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
    ramp=Input1.data(:,22);
    
    [test,peaks]=findpeaks(Input1.data(:,17));
    
    for i=1:length(peaks)
        
        time=[-200:0.5:2500];
       
        if test(i)==2
            force1r_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f2(peaks(i)-360);
            force1r_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f3(peaks(i)-360);
            force1r_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f4(peaks(i)-360);
            force1r_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f5(peaks(i)-360);
            force1r_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f6(peaks(i)-360);
            
            weight1r(p0,:)=w1(peaks(i)-400:peaks(i)+TRIAL_WINDOW_POST_PERT)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
            
            ta_emg(p0, :) = ta(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
            ta_emg(p0, :) = abs(ta_emg(p0, :)-off_TA);
            ta_emg(p0, :) = filtfilt(d3, ta_emg(p0, :));
            
            sol_emg(p0, :) = sol(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
            sol_emg(p0, :) = abs(sol_emg(p0, :)-off_SOL);
            sol_emg(p0, :) = filtfilt(d3, sol_emg(p0, :));
            
            pl_emg(p0, :)  = pl(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
            pl_emg(p0, :) = abs(pl_emg(p0, :)-off_TA);
            pl_emg(p0, :) = filtfilt(d3, pl_emg(p0, :));
            
            gca_emg(p0, :) = gca(peaks(i)+TRIAL_WINDOW_PRE_PERT:peaks(i)+TRIAL_WINDOW_POST_PERT);
            gca_emg(p0, :) = abs(gca_emg(p0, :)-off_TA);
            gca_emg(p0, :) = filtfilt(d3, gca_emg(p0, :));
            
%             p1r_plat_pos(p0,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
%             p1r_foot_pos(p0,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
            phase1r(p0,:)=ramp(peaks(i)-400:peaks(i)+2000);
            p0=p0+1;
        end        
    end
    
  end
    %True Gait Phase
    max_weight = mean(max(weight1r'));
    start_index_vec_r = zeros(1, NUM_TRAINING_TRIALS); 
    end_index_vec_r = zeros(1, NUM_TRAINING_TRIALS);

    for i = 1:NUM_TRAINING_TRIALS
    start_index = min(find(weight1r(i,:) > WEIGHT_PERCENTAGE_CUTOFF*max_weight)); %This is
    end_index = max(find(weight1r(i,:) > WEIGHT_PERCENTAGE_CUTOFF*max_weight));
    start_index_vec_r(i) = start_index;
    end_index_vec_r(i) = end_index;
    end

    stance_phase_duration_vec = end_index_vec_r - start_index_vec_r;

    mean_start_index = mean(start_index_vec_r);
    mean_stance_phase_duration = mean(stance_phase_duration_vec);
    median_stance_phase_duration = median(stance_phase_duration_vec);
end

TA_NORM  = max(max(ta_emg(1:NUM_TRAINING_TRIALS, -TRIAL_WINDOW_PRE_PERT:-TRIAL_WINDOW_PRE_PERT+median_stance_phase_duration)'));
SOL_NORM = max(max(sol_emg(1:NUM_TRAINING_TRIALS, -TRIAL_WINDOW_PRE_PERT:-TRIAL_WINDOW_PRE_PERT+median_stance_phase_duration)'));
PL_NORM  = max(max(pl_emg(1:NUM_TRAINING_TRIALS, -TRIAL_WINDOW_PRE_PERT:-TRIAL_WINDOW_PRE_PERT+median_stance_phase_duration)'));
GCA_NORM = max(max(gca_emg(1:NUM_TRAINING_TRIALS, -TRIAL_WINDOW_PRE_PERT:-TRIAL_WINDOW_PRE_PERT+median_stance_phase_duration)'));

end


