function [diff_pos,diff_torque, mean_plat_pos_profiles, bio_factors_struct] = get_pos_torque_diff(p0_pos_profiles, pert_pos_profiles, ...
                                                      p0_torque_profiles, pert_torque_profiles, ...
                                                      diff_plat_pos_profiles, ...
                                                      cop_vals, ...
                                                      weight_vals, ...
                                                      emg_vals) %This will be a struct since passing in each
                                                                    %individual set of EMG curves would create a 
                                                                    %long function
                                                                    
                                                                    
                                                                   
                                                                    
                                                                    
                                                                    



                                                  
%% New Section Here

NUM_OF_CLOSEST_MSE_PREPURT = 5;


%% ---- Nominal Profiles and BioMech Factors ---- %%
mean_p0_pos_profiles = [];
mean_p0_torque_profiles = [];

mean_pert_cop_vals = [];
mean_pert_weight_vals = [];
mean_pert_TA_profiles = [];
mean_pert_SOL_profiles = [];
mean_pert_PL_profiles = [];
mean_pert_GCA_profiles = [];
                                                  
RESAMPLE_COUNT = 100;
bs_selection_count = round(0.6*size(p0_pos_profiles, 1));
for i = 1:RESAMPLE_COUNT
   %Random selection of Foot Position Curves 
    [pos_data_sample, sample_ind] = datasample(p0_pos_profiles, bs_selection_count);
    
    pos_mean = mean(pos_data_sample, 1);
    torque_mean = mean(p0_torque_profiles(sample_ind, :), 1);
    
    mean_p0_pos_profiles(i, :) = pos_mean;
    mean_p0_torque_profiles(i, :) = torque_mean;
   
end


%% ---- Perturbation Profiles ---- %%
bs_selection_count = round(0.6*size(pert_pos_profiles, 1));
for i = 1:RESAMPLE_COUNT
   %Random selection of Foot Position Curves 
    [pos_data_sample, sample_ind] = datasample(pert_pos_profiles, bs_selection_count);
    
    pos_mean = mean(pos_data_sample, 1);
    torque_mean = mean(pert_torque_profiles(sample_ind, :), 1);
    diff_plat_pos_mean = mean(diff_plat_pos_profiles(sample_ind, :), 1);
    
    mean_pert_pos_profiles(i, :) = pos_mean;
    mean_pert_torque_profiles(i, :) = torque_mean;
    mean_plat_pos_profiles(i, :) = diff_plat_pos_mean;
    
    cop_mean = mean(cop_vals(sample_ind, :));
    weight_mean = mean(weight_vals(sample_ind));
    emg_TA_mean = mean(emg_vals.TA(sample_ind, :));
    emg_PL_mean = mean(emg_vals.PL(sample_ind, :));
    emg_SOL_mean = mean(emg_vals.SOL(sample_ind, :));
    emg_GCA_mean = mean(emg_vals.GCA(sample_ind, :));
    
    
    mean_pert_cop_profiles(i, :) = cop_mean;
    mean_pert_weight_vals(i, :) = weight_mean;
    mean_pert_TA_profiles(i, :) = emg_TA_mean;
    mean_pert_SOL_profiles(i, :) = emg_SOL_mean;
    mean_pert_PL_profiles(i, :) = emg_PL_mean;
    mean_pert_GCA_profiles(i, :) = emg_GCA_mean;
    mean_pert_ang(i) = pos_mean(100);
end


%% ---- Find NUM_OF_CLOSEST_MSE_PREPURT# of closest Curves ---- %%
p0_pre_perturbation_pos_segs = mean_p0_pos_profiles(:, 1:100) - mean_p0_pos_profiles(:, 1);
p0_pre_perturbation_trq_segs = mean_p0_torque_profiles(:, 1:100) - mean_p0_torque_profiles(:, 1);

pre_pert_pos_segments = mean_pert_pos_profiles(:, 1:100) - mean_pert_pos_profiles(:, 1);
pre_pert_trq_segments = mean_pert_torque_profiles(:, 1:100) - mean_pert_torque_profiles(:, 1);


MSE_NP_trials_vec = zeros(size(mean_pert_pos_profiles, 1), NUM_OF_CLOSEST_MSE_PREPURT);

for pert_trial = 1:size(mean_pert_pos_profiles, 1)
    pos_err_vecs = (p0_pre_perturbation_pos_segs - pre_pert_pos_segments(pert_trial, :));
    MSE_squared_pos = dot(pos_err_vecs, pos_err_vecs, 2);
    
    trq_err_vecs = (p0_pre_perturbation_trq_segs - pre_pert_trq_segments(pert_trial, :));
    MSE_squared_trq = dot(trq_err_vecs, trq_err_vecs, 2);
    
    MSE_squared_sum = MSE_squared_pos + MSE_squared_trq;
    
    %[min_err, ideal_non_pert_index] = min(MSE_squared);
    [sorted_values, sorted_indices] = sort(MSE_squared_sum);
    
    MSE_NP_trials_vec(pert_trial, :) = sorted_indices(1:NUM_OF_CLOSEST_MSE_PREPURT)';
end


%% ---- Find the Differences between the NonPerturbation and Perturbation Curves ---- %%
%Differential Position
for i = 1:size(mean_pert_pos_profiles, 1)
    pert_curve = mean_pert_pos_profiles(i, :) - mean_pert_pos_profiles(i, 100);
    non_pert_curve = mean(mean_p0_pos_profiles(MSE_NP_trials_vec(i, :), :), 1) - mean(mean_p0_pos_profiles(MSE_NP_trials_vec(i, :), 100), 1);
    diff_pos(i, :) = pert_curve - non_pert_curve;
end

%Differential Torque
for i = 1:size(mean_pert_pos_profiles, 1)
    pert_curve = mean_pert_torque_profiles(i, :) - mean_pert_torque_profiles(i, 100);
    non_pert_curve = mean(mean_p0_torque_profiles(MSE_NP_trials_vec(i, :), :), 1) - mean(mean_p0_torque_profiles(MSE_NP_trials_vec(i, :), 100), 1);
    diff_torque(i, :) = pert_curve - non_pert_curve;
end

% %% ---- Average NonPerturbation Biomechanical Factors Together ---- %%
% %CoP and Weight have NOT had the BioMechanical factor found
% for i = 1:size(mean_pert_pos_profiles, 1)
%     
%     %EMG signals have
%     total_mean_EMG_TA(i)  = mean(mean_pert_TA_profiles(MSE_NP_trials_vec(i, :)));
%     total_mean_EMG_PL(i)  = mean(mean_pert_PL_profiles(MSE_NP_trials_vec(i, :)));
%     total_mean_EMG_SOL(i) = mean(mean_p0_SOL_profiles(MSE_NP_trials_vec(i, :)));  
%     total_mean_EMG_GCA(i) = mean(mean_pert_GCA_profiles(MSE_NP_trials_vec(i, :)));
% 
% end

bio_factors_struct.CoP = mean_pert_cop_profiles;
bio_factors_struct.Weight = mean_pert_weight_vals;

EMG_DATA_OUT.TA  = mean_pert_TA_profiles;
EMG_DATA_OUT.PL  = mean_pert_PL_profiles;
EMG_DATA_OUT.SOL = mean_pert_SOL_profiles;
EMG_DATA_OUT.GCA = mean_pert_GCA_profiles;

bio_factors_struct.EMG = EMG_DATA_OUT;

bio_factors_struct.ankle_ang = mean_pert_ang;


% %----MSE NonPerturb Average Reduction Method----%                                                    
% mean_p0_pos_profiles = [];
% mean_p0_torque_profiles = [];
%                                                   
% RESAMPLE_COUNT = 40;
% bs_selection_count = round(0.6*size(p0_pos_profiles, 1))
% for i = 1:RESAMPLE_COUNT
%    %Random selection of Foot Position Curves 
%     [pos_data_sample, sample_ind] = datasample(p0_pos_profiles, bs_selection_count);
%     
%     pos_mean = mean(pos_data_sample);
%     torque_mean = mean(p0_torque_profiles(sample_ind, :));
%     
%     mean_p0_pos_profiles(i, :) = pos_mean;
%     mean_p0_torque_profiles(i, :) = torque_mean;
% end
% 
% p0_pre_perturbation_segments = mean_p0_pos_profiles(:, 1:100) - mean_p0_pos_profiles(:, 1);
% pre_perturbation_segments = pert_pos_profiles(:, 1:100) - pert_pos_profiles(:, 1);
% 
% MSE_NP_trials_vec = zeros(size(pert_pos_profiles, 1),1);
% 
% for pert_trial = 1:size(pert_pos_profiles, 1)
%     err_vecs = (p0_pre_perturbation_segments - pre_perturbation_segments(pert_trial, :));
%     MSE_squared = dot(err_vecs, err_vecs, 2);
%     
%     [min_err, ideal_non_pert_index] = min(MSE_squared);
%     
%     MSE_NP_trials_vec(pert_trial) = ideal_non_pert_index;
% end
% 
% %Differential Position
% for i = 1:size(pert_pos_profiles, 1)
%     pert_curve = pert_pos_profiles(i, :) - pert_pos_profiles(i, 100);
%     non_pert_curve = mean_p0_pos_profiles(MSE_NP_trials_vec(i), :) - mean_p0_pos_profiles(MSE_NP_trials_vec(i), 100);
%     diff_pos(i, :) = pert_curve - non_pert_curve;
% end
% 
% %Differential Torque
% for i = 1:size(pert_pos_profiles, 1)
%     pert_curve = pert_torque_profiles(i, :) - pert_torque_profiles(i, 100);
%     non_pert_curve = mean_p0_torque_profiles(MSE_NP_trials_vec(i), :) - mean_p0_torque_profiles(MSE_NP_trials_vec(i), 100);
%     diff_torque(i, :) = pert_curve - non_pert_curve;
% end


%%----MSE Reduction Method----%
% p0_pre_perturbation_segments = p0_pos_profiles(:, 1:100) - p0_pos_profiles(:, 1);
% pre_perturbation_segments = pert_pos_profiles(:, 1:100) - pert_pos_profiles(:, 1);
% 
% MSE_NP_trials_vec = zeros(size(pert_pos_profiles, 1),1);
% 
% for pert_trial = 1:size(pert_pos_profiles, 1)
%     err_vecs = (p0_pre_perturbation_segments - pre_perturbation_segments(pert_trial, :));
%     MSE_squared = dot(err_vecs, err_vecs, 2);
%     
%     [min_err, ideal_non_pert_index] = min(MSE_squared);
%     
%     MSE_NP_trials_vec(pert_trial) = ideal_non_pert_index;
% end
% 
% %Differential Position
% for i = 1:size(pert_pos_profiles, 1)
%     pert_curve = pert_pos_profiles(i, :) - pert_pos_profiles(i, 100);
%     non_pert_curve = p0_pos_profiles(MSE_NP_trials_vec(i), :) - p0_pos_profiles(MSE_NP_trials_vec(i), 100);
%     diff_pos(i, :) = pert_curve - non_pert_curve;
% end
% 
% %Differential Torque
% for i = 1:size(pert_pos_profiles, 1)
%     pert_curve = pert_torque_profiles(i, :) - pert_torque_profiles(i, 100);
%     non_pert_curve = p0_torque_profiles(MSE_NP_trials_vec(i), :) - p0_torque_profiles(MSE_NP_trials_vec(i), 100);
%     diff_torque(i, :) = pert_curve - non_pert_curve;
% end

end

