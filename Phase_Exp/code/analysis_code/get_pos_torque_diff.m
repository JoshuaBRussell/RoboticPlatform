function [diff_pos,diff_torque, mean_plat_pos_profiles] = get_pos_torque_diff(p0_pos_profiles, pert_pos_profiles, ...
                                                      p0_torque_profiles, pert_torque_profiles, ...
                                                      diff_plat_pos_profiles)



                                                  
%% New Section Here
                                                  
mean_p0_pos_profiles = [];
mean_p0_torque_profiles = [];
                                                  
RESAMPLE_COUNT = 100;
bs_selection_count = round(0.6*size(p0_pos_profiles, 1));
for i = 1:RESAMPLE_COUNT
   %Random selection of Foot Position Curves 
    [pos_data_sample, sample_ind] = datasample(p0_pos_profiles, bs_selection_count);
    
    pos_mean = mean(pos_data_sample);
    torque_mean = mean(p0_torque_profiles(sample_ind, :));
    
    mean_p0_pos_profiles(i, :) = pos_mean;
    mean_p0_torque_profiles(i, :) = torque_mean;
end

for i = 1:RESAMPLE_COUNT
   %Random selection of Foot Position Curves 
    [pos_data_sample, sample_ind] = datasample(pert_pos_profiles, bs_selection_count);
    
    pos_mean = mean(pos_data_sample);
    torque_mean = mean(pert_torque_profiles(sample_ind, :));
    diff_plat_pos_mean = mean(diff_plat_pos_profiles(sample_ind, :));
    
    mean_pert_pos_profiles(i, :) = pos_mean;
    mean_pert_torque_profiles(i, :) = torque_mean;
    mean_plat_pos_profiles(i, :) = diff_plat_pos_mean;
end

p0_pre_perturbation_segments = mean_p0_pos_profiles(:, 1:100) - mean_p0_pos_profiles(:, 1);
pre_perturbation_segments = mean_pert_pos_profiles(:, 1:100) - mean_pert_pos_profiles(:, 1);

MSE_NP_trials_vec = zeros(size(mean_pert_pos_profiles, 1),1);

for pert_trial = 1:size(mean_pert_pos_profiles, 1)
    err_vecs = (p0_pre_perturbation_segments - pre_perturbation_segments(pert_trial, :));
    MSE_squared = dot(err_vecs, err_vecs, 2);
    
    [min_err, ideal_non_pert_index] = min(MSE_squared);
    
    MSE_NP_trials_vec(pert_trial) = ideal_non_pert_index;
end

%Differential Position
for i = 1:size(mean_pert_pos_profiles, 1)
    pert_curve = mean_pert_pos_profiles(i, :) - mean_pert_pos_profiles(i, 100);
    non_pert_curve = mean_p0_pos_profiles(MSE_NP_trials_vec(i), :) - mean_p0_pos_profiles(MSE_NP_trials_vec(i), 100);
    diff_pos(i, :) = pert_curve - non_pert_curve;
end

%Differential Torque
for i = 1:size(mean_pert_pos_profiles, 1)
    pert_curve = mean_pert_torque_profiles(i, :) - mean_pert_torque_profiles(i, 100);
    non_pert_curve = mean_p0_torque_profiles(MSE_NP_trials_vec(i), :) - mean_p0_torque_profiles(MSE_NP_trials_vec(i), 100);
    diff_torque(i, :) = pert_curve - non_pert_curve;
end

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

