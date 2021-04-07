%% Run this to have the analysis code run for multiple subjects.
%% Analyzing all of the subjects after a code update becomes much 
%% easier by just calling this. 

clc;
clear;


GROUP_DATA_FOLDER_REL_LOC = "./../../data/";

multi_subj = 1;

SUBJ_DATA_DIRS = {'Carl_121720', ...
                  'Vu_121420', ...
                  'Emily_122020', ...
                  'Ashley_011721', ...
                  'Carly_020121', ...
                  'Lily_020221', ...
                  'Matt_022821', ...
                  'James_030221', ...
                  'Kwanghee_030321'};

for subjects = 1:length(SUBJ_DATA_DIRS)
    
    CURR_SUBJ_REL_LOC = strcat(GROUP_DATA_FOLDER_REL_LOC, SUBJ_DATA_DIRS{subjects}, '/');
    walking_study_analyze;
end
    