%% Run this to have the analysis code run for multiple subjects.
%% Analyzing all of the subjects after a code update becomes much 
%% easier by just calling this. 

clc;
clear;


GROUP_DATA_FOLDER_REL_LOC = "./../../data/";
RESULTS_DIR = './results/';
multi_subj = 1;

SUBJ_DATA_DIRS = {'Carl_121720', ...
                  'Vu_121420', ...
                  'Emily_122020', ...
                  'Ashley_011721', ...
                  'Carly_020121', ...
                  'Lily_020221', ...
                  'Matt_022821', ...
                  'James_030221', ...
                  'Kwanghee_030321', ...
                  'Anna_042321'};%, ...
                  %'Ian_041021'};

for subjects = 1:length(SUBJ_DATA_DIRS)
    
    CURR_SUBJ_REL_LOC = strcat(GROUP_DATA_FOLDER_REL_LOC, SUBJ_DATA_DIRS{subjects}, '/');
    walking_study_analyze;
end

clear;
clc;

%% ---- Process Saved Data From Each Subject ---- %%
RESULTS_DIR = './results/';

SUBJ_DATA_DIRS = {'Carl_121720', ...
                  'Vu_121420', ...
                  'Emily_122020', ...
                  'Ashley_011721', ...
                  'Carly_020121', ...
                  'Lily_020221', ...
                  'Matt_022821', ...
                  'James_030221', ...
                  'Kwanghee_030321', ...
                  'Anna_042321'};%, ...
                  %'Ian_041021'};

figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(sign(bio_factors_p1.CoP).*sqrt(abs(bio_factors_p1.CoP)), regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(sign(bio_factors_p2.CoP).*sqrt(abs(bio_factors_p2.CoP)), regress_coeffs_p2(:, 1), 'r');
    scatter(sign(bio_factors_p3.CoP).*sqrt(abs(bio_factors_p3.CoP)), regress_coeffs_p3(:, 1), 'g');
    scatter(sign(bio_factors_p4.CoP).*sqrt(abs(bio_factors_p4.CoP)), regress_coeffs_p4(:, 1), 'b'); %hold off;

end

legend(["31";"44"; "57"; "18"])
title("K vs Sqrt(CoP)");
xlabel("Sqrt(CoP)");
ylabel("Stiffness (Nm/rad)");
saveas(gcf,strcat('sqrt_copf_plot.jpg'));

%% K vs Ankle Angle
figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(bio_factors_p1.ankle_ang', regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(bio_factors_p2.ankle_ang', regress_coeffs_p2(:, 1), 'r');
    scatter(bio_factors_p3.ankle_ang', regress_coeffs_p3(:, 1), 'g');
    scatter(bio_factors_p4.ankle_ang', regress_coeffs_p4(:, 1), 'b'); 
end

legend(["31";"44"; "57"; "18"])
title("K vs Ankle Angle");
xlabel("Ankle Angle (rads)");
ylabel("Stiffness (Nm/rad)");
saveas(gcf,strcat('ankle_angle_plot.jpg'));


%% K vs EMG TA
figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(bio_factors_p1.EMG.TA, regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(bio_factors_p2.EMG.TA, regress_coeffs_p2(:, 1), 'r');
    scatter(bio_factors_p3.EMG.TA, regress_coeffs_p3(:, 1), 'g');
    scatter(bio_factors_p4.EMG.TA, regress_coeffs_p4(:, 1), 'b'); 
end

legend(["31";"44"; "57"; "18"])
title("K vs Tibialis Anterior");
xlabel("Tibialis Anterior(%MVC)");
ylabel("Stiffness (Nm/rad)");
saveas(gcf,strcat('EMG_TA_plot.jpg'));

%% K vs EMG Triceps Surae
figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    x1 = (bio_factors_p1.EMG.GCA + bio_factors_p1.EMG.SOL);
    x2 = (bio_factors_p2.EMG.GCA + bio_factors_p2.EMG.SOL);
    x3 = (bio_factors_p3.EMG.GCA + bio_factors_p3.EMG.SOL);
    x4 = (bio_factors_p4.EMG.GCA + bio_factors_p4.EMG.SOL);
    
    %figure();
    scatter(x1, regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(x2, regress_coeffs_p2(:, 1), 'r');
    scatter(x3, regress_coeffs_p3(:, 1), 'g');
    scatter(x4, regress_coeffs_p4(:, 1), 'b'); 
end

legend(["31";"44"; "57"; "18"])
title("K vs (Triceps Surae)");
xlabel("(Triceps Surae (SOL + GCA)) (%MVC)");
ylabel("Stiffness (Nm/rad)");
saveas(gcf,strcat('EMG_TS_plot.jpg'));

%% K vs BW
figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(bio_factors_p1.Weight, regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(bio_factors_p2.Weight, regress_coeffs_p2(:, 1), 'r');
    scatter(bio_factors_p3.Weight, regress_coeffs_p3(:, 1), 'g');
    scatter(bio_factors_p4.Weight, regress_coeffs_p4(:, 1), 'b'); 
end

legend(["31";"44"; "57"; "18"])
title("K vs BodyWeight");
xlabel("Body Weight (Newtons)");
ylabel("Stiffness (Nm/rad)");
saveas(gcf,strcat('BW_plot.jpg'));

%% K vs Time Since Healstrike
figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(bio_factors_p1.time_since_healstrike-0.2, regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(bio_factors_p2.time_since_healstrike-0.2, regress_coeffs_p2(:, 1), 'r');
    scatter(bio_factors_p3.time_since_healstrike-0.2, regress_coeffs_p3(:, 1), 'g');
    scatter(bio_factors_p4.time_since_healstrike-0.2, regress_coeffs_p4(:, 1), 'b'); %hold off;

end

legend(["31";"44"; "57"; "18"])
title("K vs Time");
xlabel("Time (s)");
ylabel("Stiffness (Nm/rad)");
saveas(gcf,strcat('K_vs_Time_plot.jpg'));

%% K vs Sex
 
SUBJ_SEX = {"M"; "M"; "F"; "F"; "F"; "F"; "M"; "M"; "M"; "F"};
figure();
for subj_index = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subj_index}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    plot_color = 'r';
    if SUBJ_SEX{subj_index} == "M"
        plot_color = 'b';
    end
    
    %figure();
    scatter(bio_factors_p1.time_since_healstrike-0.2, regress_coeffs_p1(:, 1), plot_color); hold on;
    scatter(bio_factors_p2.time_since_healstrike-0.2, regress_coeffs_p2(:, 1), plot_color);
    scatter(bio_factors_p3.time_since_healstrike-0.2, regress_coeffs_p3(:, 1), plot_color);
    scatter(bio_factors_p4.time_since_healstrike-0.2, regress_coeffs_p4(:, 1), plot_color); %hold off;

end

legend(["31";"44"; "57"; "18"])
title("K vs Time");
xlabel("Time (s)");
ylabel("Stiffness (Nm/rad)");
%saveas(gcf,strcat('K_vs_Time_plot.jpg'));

%% K (BodyWeight Normalized) vs Time Since Healstrike 
DATA_FOLDER_REL_LOC = "./../../data/";

subj_weight_vec = get_subj_weight(DATA_FOLDER_REL_LOC, "WEIGHT.DAT", SUBJ_DATA_DIRS);

figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(bio_factors_p1.time_since_healstrike-0.2, regress_coeffs_p1(:, 1)/subj_weight_vec(subjects), 'k'); hold on;
    scatter(bio_factors_p2.time_since_healstrike-0.2, regress_coeffs_p2(:, 1)/subj_weight_vec(subjects), 'r');
    scatter(bio_factors_p3.time_since_healstrike-0.2, regress_coeffs_p3(:, 1)/subj_weight_vec(subjects), 'g');
    scatter(bio_factors_p4.time_since_healstrike-0.2, regress_coeffs_p4(:, 1)/subj_weight_vec(subjects), 'b'); %hold off;

end

legend(["31";"44"; "57"; "18"])
title("K/BW vs Time");
xlabel("Time (s)");
ylabel("Stiffness/BW ((Nm/rad)/Newtons)");
saveas(gcf,strcat('K_BW_norm_vs_Time_plot.jpg'));

%% K (Subj. Range Normalized) vs Time Since Healstrike 
DATA_FOLDER_REL_LOC = "./../../data/";

subj_weight_vec = get_subj_weight(DATA_FOLDER_REL_LOC, "WEIGHT.DAT", SUBJ_DATA_DIRS);

figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    subj_K = [regress_coeffs_p1(:, 1); regress_coeffs_p2(:, 1); regress_coeffs_p3(:, 1); regress_coeffs_p4(:, 1)]; 
    norm_factor = (max(subj_K) - min(subj_K));
    
    %figure();
    scatter(bio_factors_p1.time_since_healstrike-0.2, regress_coeffs_p1(:, 1)/norm_factor, 'k'); hold on;
    scatter(bio_factors_p2.time_since_healstrike-0.2, regress_coeffs_p2(:, 1)/norm_factor, 'r');
    scatter(bio_factors_p3.time_since_healstrike-0.2, regress_coeffs_p3(:, 1)/norm_factor, 'g');
    scatter(bio_factors_p4.time_since_healstrike-0.2, regress_coeffs_p4(:, 1)/norm_factor, 'b'); %hold off;

end

legend(["31";"44"; "57"; "18"])
title("K subject range noramlized vs Time");
xlabel("Time (s)");
ylabel("Stiffness Range Normalized");
saveas(gcf,strcat('K_range_norm_vs_Time_plot.jpg'));

%% K (Subj. Range Normalized) vs Time Since Healstrike 
DATA_FOLDER_REL_LOC = "./../../data/";

subj_weight_vec = get_subj_weight(DATA_FOLDER_REL_LOC, "WEIGHT.DAT", SUBJ_DATA_DIRS);

figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    subj_K = [regress_coeffs_p1(:, 1); regress_coeffs_p2(:, 1); regress_coeffs_p3(:, 1); regress_coeffs_p4(:, 1)]; 
    norm_factor = (744); %Find manually
    
    %figure();
    scatter(bio_factors_p1.time_since_healstrike-0.2, regress_coeffs_p1(:, 1)/norm_factor, 'k'); hold on;
    scatter(bio_factors_p2.time_since_healstrike-0.2, regress_coeffs_p2(:, 1)/norm_factor, 'r');
    scatter(bio_factors_p3.time_since_healstrike-0.2, regress_coeffs_p3(:, 1)/norm_factor, 'g');
    scatter(bio_factors_p4.time_since_healstrike-0.2, regress_coeffs_p4(:, 1)/norm_factor, 'b'); %hold off;

end

legend(["31";"44"; "57"; "18"])
title("K Population range noramlized vs Time");
xlabel("Time (s)");
ylabel("Stiffness Pop. Range Normalized");
saveas(gcf,strcat('K_pop_range_norm_vs_Time_plot.jpg'));


%% ---- Regression Procedure ---- %%
RESULTS_DIR = './results/';
DATA_FOLDER_REL_LOC = "./../../data/";

SUBJ_DATA_DIRS = {'Carl_121720', ...
                  'Vu_121420', ...
                  'Emily_122020', ...
                  'Ashley_011721', ...
                  'Carly_020121', ...
                  'Lily_020221', ...
                  'Matt_022821', ...
                  'James_030221', ...
                  'Kwanghee_030321', ...
                  'Anna_042321'};%, ...
                  %'Ian_041021'};

                  
% Collect Data to Perform Regression with 
K = [];
BW_norm_K_data = [];
stiff_range_norm_K_data = [];

cop_data = [];
ang_data = [];
emg_data_TA = []; %Tibialis Anterior
emg_data_TS = []; %Triceps Surae
bw_data = [];

subj_weight_vec = get_subj_weight(DATA_FOLDER_REL_LOC, "WEIGHT.DAT", SUBJ_DATA_DIRS);

for subj_index = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subj_index}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    subj_K = [regress_coeffs_p1(:, 1); regress_coeffs_p2(:, 1); regress_coeffs_p3(:, 1); regress_coeffs_p4(:, 1)]; 
    K = [K; subj_K];
    
    subj_BW_norm_K = subj_K/subj_weight_vec(subj_index);
    BW_norm_K_data = [BW_norm_K_data; subj_BW_norm_K];
    
    subj_stiff_range_norm_K = subj_K/(max(subj_K) - min(subj_K));
    stiff_range_norm_K_data = [stiff_range_norm_K_data; subj_stiff_range_norm_K];
    
    subj_cop = [bio_factors_p1.CoP; bio_factors_p2.CoP; bio_factors_p3.CoP; bio_factors_p4.CoP];
    cop_data = [cop_data; subj_cop];
    
    subj_ang = [bio_factors_p1.ankle_ang, bio_factors_p2.ankle_ang, bio_factors_p3.ankle_ang, bio_factors_p4.ankle_ang]';
    ang_data = [ang_data; subj_ang];
    
    subj_emg_TA = [bio_factors_p1.EMG.TA; bio_factors_p2.EMG.TA; bio_factors_p3.EMG.TA; bio_factors_p4.EMG.TA];
    emg_data_TA = [emg_data_TA; subj_emg_TA];
    
    subj_emg_TS = [bio_factors_p1.EMG.SOL + bio_factors_p1.EMG.GCA; ...
                bio_factors_p2.EMG.SOL + bio_factors_p2.EMG.GCA; ...
                bio_factors_p3.EMG.SOL + bio_factors_p3.EMG.GCA; ...
                bio_factors_p4.EMG.SOL + bio_factors_p4.EMG.GCA];
    emg_data_TS = [emg_data_TS; subj_emg_TS];
    
    subj_bw = [bio_factors_p1.Weight; bio_factors_p2.Weight; bio_factors_p3.Weight; bio_factors_p4.Weight];
    bw_data = [bw_data; subj_bw];

end

pop_stiff_range_norm_K_data = K/(max(K) - min(K));

%Am going to take the sqrt of the CoP data and EMG Triceps Surae data
cop_sign_data = sign(cop_data);
cop_data = cop_sign_data.*(sqrt(abs(cop_data)));
emg_data_TS = log(emg_data_TS);

partialcorr([K, cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data]);
corr([K, cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data]);

% ---- Organize Data Bunches into Tables ---- %
all_data_K = table(cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data, K);
all_indep_wo_angle_K = table(cop_data, emg_data_TA, emg_data_TS, bw_data, K);
% ---- Create Regression Models ---- %
% Y = K
mdl_K_1 = fitlm(all_data_K); 
mdl_K_2 = fitlm(all_indep_wo_angle_K); %same as 1 above but missing ankle angle

%Y = BW_norm_K_data
mdl_BW_norm_K_1 = fitlm([cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data], BW_norm_K_data); 
mdl_BW_norm_K_2 = fitlm([cop_data, emg_data_TS, emg_data_TS, bw_data], BW_norm_K_data); %same as 1 above but missing ankle angle

%Y = stiff_range_norm_K_data
mdl_stiff_range_norm_K_1 = fitlm([cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data], stiff_range_norm_K_data); 
mdl_stiff_range_norm_K_2 = fitlm([cop_data, emg_data_TS, emg_data_TS, bw_data], stiff_range_norm_K_data); %same as 1 above but missing ankle angle

%Y = pop_stiff_range_norm_K_data
mdl_pop_stiff_range_norm_K_1 = fitlm([cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data], pop_stiff_range_norm_K_data); 
mdl_pop_stiff_range_norm_K_2 = fitlm([cop_data, emg_data_TS, emg_data_TS, bw_data], pop_stiff_range_norm_K_data); %same as 1 above but missing ankle angle
