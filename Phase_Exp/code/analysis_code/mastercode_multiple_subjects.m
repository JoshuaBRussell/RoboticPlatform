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
    scatter(sqrt(bio_factors_p1.CoP)', regress_coeffs_p1(:, 1)); hold on;
    scatter(sqrt(bio_factors_p2.CoP)', regress_coeffs_p2(:, 1));
    scatter(sqrt(bio_factors_p3.CoP)', regress_coeffs_p3(:, 1));
    scatter(sqrt(bio_factors_p4.CoP)', regress_coeffs_p4(:, 1)); %hold off;
    %title(sub_name);
    
    %saveas(gcf,strcat(RESULTS_DIR, sub_name,'copf_plot.jpg'));


end
