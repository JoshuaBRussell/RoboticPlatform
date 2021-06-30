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
    scatter(sqrt(bio_factors_p1.CoP)', regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(sqrt(bio_factors_p2.CoP)', regress_coeffs_p2(:, 1), 'r');
    scatter(sqrt(bio_factors_p3.CoP)', regress_coeffs_p3(:, 1), 'g');
    scatter(sqrt(bio_factors_p4.CoP)', regress_coeffs_p4(:, 1), 'b'); %hold off;
    %title(sub_name;
    legend(["31";"44"; "57"; "18"])
    
    %saveas(gcf,strcat(RESULTS_DIR, sub_name,'copf_plot.jpg'));


end

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
    scatter(bio_factors_p4.ankle_ang', regress_coeffs_p4(:, 1), 'b'); %hold off;
    %title(sub_name;
    legend(["31";"44"; "57"; "18"])
    
    %saveas(gcf,strcat(RESULTS_DIR, sub_name,'copf_plot.jpg'));


end

%% K vs Time Since Healstrike
figure();
for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    %figure();
    scatter(bio_factors_p1.time_since_healstrike, regress_coeffs_p1(:, 1), 'k'); hold on;
    scatter(bio_factors_p2.time_since_healstrike, regress_coeffs_p2(:, 1), 'r');
    scatter(bio_factors_p3.time_since_healstrike, regress_coeffs_p3(:, 1), 'g');
    scatter(bio_factors_p4.time_since_healstrike, regress_coeffs_p4(:, 1), 'b'); %hold off;
    %title(sub_name;
    legend(["31";"44"; "57"; "18"])
    
    %saveas(gcf,strcat(RESULTS_DIR, sub_name,'copf_plot.jpg'));


end

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
    scatter(bio_factors_p1.time_since_healstrike, regress_coeffs_p1(:, 1)/subj_weight_vec(subjects), 'k'); hold on;
    scatter(bio_factors_p2.time_since_healstrike, regress_coeffs_p2(:, 1)/subj_weight_vec(subjects), 'r');
    scatter(bio_factors_p3.time_since_healstrike, regress_coeffs_p3(:, 1)/subj_weight_vec(subjects), 'g');
    scatter(bio_factors_p4.time_since_healstrike, regress_coeffs_p4(:, 1)/subj_weight_vec(subjects), 'b'); %hold off;
    %title(sub_name;
    legend(["31";"44"; "57"; "18"])
    
    %saveas(gcf,strcat(RESULTS_DIR, sub_name,'copf_plot.jpg'));


end


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
cop_data = [];
ang_data = [];
emg_data_TA = []; %Tibialis Anterior
emg_data_TS = []; %Triceps Surae
bw_data = [];

figure();

for subjects = 1:length(SUBJ_DATA_DIRS)
    temp_cell_array = split(SUBJ_DATA_DIRS{subjects}, '_');
    sub_name = temp_cell_array{1};
    curr_results_dir = strcat(RESULTS_DIR, sub_name, '/');

    load(strcat(curr_results_dir, sub_name, "_bootstrap_vars.mat"));
    
    subj_K = [regress_coeffs_p1(:, 1); regress_coeffs_p2(:, 1); regress_coeffs_p3(:, 1); regress_coeffs_p4(:, 1)]; 
    K = [K; subj_K];
    
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

partialcorr([K, cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data]);
corr([K, cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data]);

mdl_1 = fitlm([cop_data, emg_data_TA, emg_data_TS, bw_data, ang_data], K);
mdl_2 = fitlm([cop_data, emg_data_TS, emg_data_TS, bw_data], K);