close all
clear all
mvc_evaluation;
% Insert subject initial and name. Make sure it matches the format for naming

WEIGHT_PERCENTAGE_CUTOFF = 0.025;

NUM_TRAINING_TRIALS = 10;
[num,den] = butter(2,10/1000);
[filt_num, filt_den] = butter(2,10/1000);
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

for trials=1:1
    
  if(ismember(trials,exclude)==0)  
    if(trials<10)
        h = fopen(strcat('GRF','.dat'));
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
 
  
    w1=filtfilt(d1,Input1.data(:,18));
  
    flag=Input1.data(:,17);
    foot_pos_data=Input1.data(:,13);
    foot_const=mean(foot_pos_data);
    foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
    plat_pos_data=Input1.data(:,14);
    plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
    ramp=Input1.data(:,22);
    
    [test,peaks]=findpeaks(Input1.data(:,17));
    
    for i=1:length(peaks)
        
        time=[-200:0.5:1000];
       
        if test(i)==2
            force1r_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f2(peaks(i)-360);
            force1r_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f3(peaks(i)-360);
            force1r_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f4(peaks(i)-360);
            force1r_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f5(peaks(i)-360);
            force1r_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f6(peaks(i)-360);
            
            weight1r(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
            
            p1r_plat_pos(p0,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
            p1r_foot_pos(p0,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
            phase1r(p0,:)=ramp(peaks(i)-400:peaks(i)+2000);
            p0=p0+1;
        end        
    end
    
    fclose all;
    clear live_data Input1 peaks Img;
    
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



%% Plot for Verification %%
figure();

plot(weight1r'); hold on;
xline(mean_start_index);
xline(mean_start_index + mean_stance_phase_duration);


disp("Mean Stance Phase Duration:");
disp(["Samples: ", num2str(mean_stance_phase_duration)]);
disp(["Milliseconds: ", num2str(mean_stance_phase_duration/2)]); %1/2 ms per sample

disp(["50% Point (Samples): ", num2str(mean_stance_phase_duration/2)]);
disp(["50% Point (ms): ", num2str(mean_stance_phase_duration/4)]);






