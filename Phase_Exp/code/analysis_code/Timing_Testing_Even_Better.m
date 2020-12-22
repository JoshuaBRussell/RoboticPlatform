
%% Setup, Filtering and Data Visualization settings
close all
clear all
% Getting the MVC values
mvc_evaluation;
% Insert subject initial and name.
%Make sure it matches the format for naming
sub_initial='M';
sub_name='Morgan';
SUB_NAME = sub_name;
%Add number of perturbation you actually ran
num_pert=40;
% insert lower limit of inertia of foot in the fit
% u_lim is the upper limit of the inertia and lim
% is the lower limit
lim=0.01;
u_lim=0.2;
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
% Add the trials you want to exclude in here
exclude=[];%exclude=[1,2];
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
        % obtaining data from channels
        pert_torque=filtfilt(d1,Input1.data(:,7));
        f1=Input1.data(:,9)*53.4;
        f1=filtfilt(d1,f1);
        f2=Input1.data(:,10)*53.4;
        f2=filtfilt(d1,f2);
        f3=Input1.data(:,11)*53.4;
        f3=filtfilt(d1,f3);
        f4=Input1.data(:,12)*53.4;
        f4=filtfilt(d1,f4);
        f5=Input1.data(:,20)*53.4/2;
        f5=filtfilt(d1,f5);
        f6=Input1.data(:,21)*53.4/2;
        f6=filtfilt(d1,f6);
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
        rigid_phase_tot=Input1.data(:,15);
        haptic_phase_tot=Input1.data(:,16);
        perturb_start=Input1.data(:,22);
        foot_pos_data=filtfilt(d1,Input1.data(:,13));
        foot_pos_data=((foot_pos_data-mean(foot_pos_data))...
            *DP_foot_gonio*pi/180);
        plat_pos_data=filtfilt(d1,Input1.data(:,14));
        plat_pos_data=((plat_pos_data-mean(plat_pos_data))...
            *DP_plat_gonio*pi/180);
        % 17 records an impulse everytime a
        % perturbation occurs with diff amplitudes
        [test,peaks]=findpeaks(Input1.data(:,17));
        for i=1:length(peaks)
     % Each test value corresponds to a different pert       
            time=[-200:0.5:1000];
            if test(i)==1
                force1_1(p1,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_2(p1,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_3(p1,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_4(p1,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_5(p1,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_6(p1,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                weight1(p1,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p1_plat_torque(p1,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p1_plat_pos(p1,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p1_foot_pos(p1,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p1_phase(p1,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p1_pert(p1,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img1_pos(p1)=getmin(peaks(i),img_st,Img);
                cop1(p1,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p1_pert(p1,:)));
                p1_peakst(p1)=b;
                p1=p1+1;
                
            end
            if test(i)==2
                force2_1(p2,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_2(p2,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_3(p2,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_4(p2,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_5(p2,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_6(p2,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                
                weight2(p2,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p2_plat_torque(p2,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p2_plat_pos(p2,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p2_foot_pos(p2,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p2_phase(p2,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p2_pert(p2,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img2_pos(p2)=getmin(peaks(i),img_st,Img);
                cop2(p2,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p2_pert(p2,:)));

                p2_peakst(p2)=b;
                p2=p2+1;
            end
            if test(i)==3
                force3_1(p3,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_2(p3,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_3(p3,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_4(p3,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_5(p3,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_6(p3,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                
                weight3(p3,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p3_plat_torque(p3,:)=1*(pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50));
                p3_plat_pos(p3,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p3_foot_pos(p3,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p3_phase(p3,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p3_pert(p3,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img3_pos(p3)=getmin(peaks(i),img_st,Img);
                cop3(p3,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p3_pert(p3,:)));
                p3_peakst(p3)=b;
                p3=p3+1;
            end
            if test(i)==4
                force0_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                
                ta_emg(p0,:)=ta(peaks(i)-400:peaks(i)+2000);
                sol_emg(p0,:)=sol(peaks(i)-400:peaks(i)+2000);
                pl_emg(p0,:)=pl(peaks(i)-400:peaks(i)+2000);
                gca_emg(p0,:)=gca(peaks(i)-400:peaks(i)+2000);
                weight4(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p0_plat_torque(p0,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p0_plat_pos(p0,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p0_foot_pos(p0,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p0_phase(p0,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p0_pert(p0,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img0_pos(p0)=getmin(peaks(i),img_st,Img);
                cop4(p0,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p0_pert(p0,:)));

                p0_peakst(p0)=400;

                [a,b]=min(diff(p0_phase(p0,:)));

                p0_peakend(p0)=b;
                p0=p0+1;
            end
            if test(i)==5
                force4_1(p4,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_2(p4,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_3(p4,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_4(p4,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_5(p4,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_6(p4,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                weight3(p4,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p4_plat_torque(p4,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p4_plat_pos(p4,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p4_foot_pos(p4,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p4_phase(p4,:)=haptic_phase_tot(peaks(i)-400:peaks(i)+2000);
                p4_pert(p4,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img4_pos(p4)=getmin(peaks(i),img_st,Img);
                [a,b]=findpeaks(diff(p4_pert(p4,:)));

                p4_peakst(p4)=b;
        
                p4=p4+1;
            end
            
        end
        
%         fclose all;
%         clear live_data Input1 peaks Img;
        
    end
    
end

%% Perturbation Comparisons
SAMPLES_PER_MS = 2;

%Get Perturbation Timings
p1_pert_start_in_samples = p1_peakst - 400;
p2_pert_start_in_samples = p2_peakst - 400;
p3_pert_start_in_samples = p3_peakst - 400;
p4_pert_start_in_samples = p4_peakst - 400;

p1_pert_start_in_ms = (1/SAMPLES_PER_MS)*p1_peakst;
p2_pert_start_in_ms = (1/SAMPLES_PER_MS)*p2_peakst;
p3_pert_start_in_ms = (1/SAMPLES_PER_MS)*p3_peakst;
p4_pert_start_in_ms = (1/SAMPLES_PER_MS)*p4_peakst;

%Check for Rigid Percentage Timings


%Compare the Two


%% No-Perturbation Check

%True Gait Phase
nominal_weight = 789;
start_index_vec = zeros(1, 4); 
end_index_vec = zeros(1, 4);
for i = 1:40
    start_index = min(find(weight4(i,:) > 0.01*nominal_weight))
    end_index = max(find(weight4(i,:) > 0.01*nominal_weight))
    start_index_vec(i) = start_index;
    end_index_vec(i) = end_index;
    
    point_1(i) = floor(end_index-start_index)*0.18 ;
    point_2(i) = floor(end_index-start_index)*0.31 ;
    point_3(i) = floor(end_index-start_index)*0.44 ;
    point_4(i) = floor(end_index-start_index)*0.57 ;
end
%A.N.N Predicted Gait Phase
for i = 1:p0-1
    pred_point_1(i) = min(find(p0_phase(i, 400:end) > 18));
    pred_point_2(i) = min(find(p0_phase(i, 400:end) > 31));
    pred_point_3(i) = min(find(p0_phase(i, 400:end) > 44));
    pred_point_4(i) = min(find(p0_phase(i, 400:end) > 57));
end



figure();
%Since the trigger point isn't the true start of stance, we need to find
%the time difference (in samples) between the true start and the trigger
%point.
stance_start_to_trigger = 400 - start_index_vec;
sim_pert_1_point = stance_start_to_trigger + pred_point_1;
sim_pert_2_point = stance_start_to_trigger + pred_point_2;
sim_pert_3_point = stance_start_to_trigger + pred_point_3;
sim_pert_4_point = stance_start_to_trigger + pred_point_4;

%Assumptions:
%Stride cycle time to be 100bpm -> 1.6667Hz -> 0.6s
%Stance phase is 60% of the cycle time
%2000Hz sampling freq.
time_pert_1_point = stance_start_to_trigger + 260;
time_pert_2_point = stance_start_to_trigger + 446;
time_pert_3_point = stance_start_to_trigger + 634;
time_pert_4_point = stance_start_to_trigger + 821;

p0_sample_length = end_index_vec-start_index_vec;
stance_phase_duration = p0_sample_length * (1/2000);
%% Normalized (Time) Weight Plot
for i = 1:length(start_index_vec)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', weight4(i, start_index_vec(i):end_index_vec(i)), 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;

end
figure();
plot(0:1/max(p0_sample_length):1, mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

formatSpec = 'Weight Profile %s';
str = sprintf(formatSpec,SUB_NAME);
title(str)
xlabel('Stance Phase (%)');
ylabel('Weight');
hold on;
plot(sim_pert_1_point./p0_sample_length, 500*ones(1, length(sim_pert_1_point)), 'rx'); hold on;
plot(sim_pert_2_point./p0_sample_length, 500*ones(1, length(sim_pert_2_point)), 'gx'); hold on;
plot(sim_pert_3_point./p0_sample_length, 500*ones(1, length(sim_pert_3_point)), 'bx'); hold on;
plot(sim_pert_4_point./p0_sample_length, 500*ones(1, length(sim_pert_4_point)), 'rx'); hold on;
%legend(["Mean (Time Normalized) Weight", "Predicted Perturbation @ 18%", "Predicted Perturbation @ 31%", "Predicted Perturbation @ 44%", "Predicted Perturbation @ 57%"]);

plot(time_pert_1_point./p0_sample_length, 400*ones(1, length(sim_pert_1_point)), 'ro'); hold on;
plot(time_pert_2_point./p0_sample_length, 400*ones(1, length(sim_pert_2_point)), 'go'); hold on;
plot(time_pert_3_point./p0_sample_length, 400*ones(1, length(sim_pert_3_point)), 'bo'); hold on;
plot(time_pert_4_point./p0_sample_length, 400*ones(1, length(sim_pert_4_point)), 'ro'); hold on;

legend(["- GRF","x ANN","o Timing Based"]);

% figure()
% stem(pred_point_1, 500*ones(size(pred_point_1))); hold on;
% stem(pred_point_2, 500*ones(size(pred_point_2))); hold on;
% stem(pred_point_3, 500*ones(size(pred_point_3))); hold on;
% stem(pred_point_4, 500*ones(size(pred_point_4))); hold on;
% plot(weight4(:, 250:end)');
% xlabel("Time (0.5 ms)");


% pred_err_1 = point_1 - pred_point_1;
% pred_err_2 = point_2 - pred_point_2;
% pred_err_3 = point_3 - pred_point_3;
% pred_err_4 = point_4 - pred_point_4;
% 
% pred_err_1_in_ms = (1/SAMPLES_PER_MS)*pred_err_1;
% pred_err_2_in_ms = (1/SAMPLES_PER_MS)*pred_err_2;
% pred_err_3_in_ms = (1/SAMPLES_PER_MS)*pred_err_3;
% pred_err_4_in_ms = (1/SAMPLES_PER_MS)*pred_err_4;
% 
% pred_err_1_mean = mean(pred_err_1_in_ms);
% pred_err_2_mean = mean(pred_err_2_in_ms);
% pred_err_3_mean = mean(pred_err_3_in_ms);
% pred_err_4_mean = mean(pred_err_4_in_ms);
% 
% pred_err_1_std = std(pred_err_1_in_ms);
% pred_err_2_std = std(pred_err_2_in_ms);
% pred_err_3_std = std(pred_err_3_in_ms);
% pred_err_4_std = std(pred_err_4_in_ms);
% 
% pred_err_means = [pred_err_1_mean, pred_err_2_mean, pred_err_3_mean, pred_err_4_mean]
% pred_err_stds = [pred_err_1_std, pred_err_2_std, pred_err_3_std, pred_err_4_std]



% %% Stance Phase Duration Plots
% nominal_weight = 789;
% for i = 1:40
%     
%     t = [-400:1:2000];
%     t = t*(1/2000);
%     start_index = min(find(weight4(i,:) > 0.01*nominal_weight))
%     end_index = max(find(weight4(i,:) > 0.01*nominal_weight))
%  
%     %figure();
%     plot(t, weight4(i,:))
%     hold on;
%     stem(t(start_index), 500, 'r')
%     hold on;
%     stem(t(end_index), 500, 'r'); hold on;
%     
%     stem(t(p0_peakst(i)), 500, 'b'); hold on;
%     stem(t(p0_peakend(i)), 500, 'b');
%     
%     legend(["Blue - Current Method", "Red - 1% Nominal Weight"]);
%     ylabel("Weight (N)");
%     xlabel("time (s)");
%     
% end 