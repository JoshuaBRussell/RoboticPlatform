function [mean_coeff, CI_coeff, goodness_of_fit, remaining_after_VAF] = regress_after_bootstrap(diff_foot_pos, diff_foot_vel, ...
                                                           diff_foot_acc, diff_plat_torque, ...
                                                           A, B)

                                                       
regress_coeffs = zeros(size(diff_foot_pos, 1), 3);
bootstrap_sample_VAF_vec = zeros(size(diff_foot_pos, 1), 1);
for i = 1:size(diff_foot_pos, 1)
    data_matrix = [diff_foot_pos(i, :)', ...
                   diff_foot_vel(i, :)', ...
                   diff_foot_acc(i, :)'];
               
    torque_vec =  diff_plat_torque(i, :)';
    regression_coeff = lsqlin(data_matrix, torque_vec, A, B);
    regress_coeffs(i, :) = regression_coeff;
    
    
    %Get a goodness of fit as well
    vartor=var((diff_plat_torque(i, :)));
    varimp=var((diff_plat_torque(i, :))-((diff_foot_pos(i, :))*regression_coeff(1) ... 
                                       + (diff_foot_vel(i, :))*regression_coeff(2) ...
                                       + (diff_foot_acc(i, :))*regression_coeff(3)));
    goodness_of_fit = 100*(1-(varimp/vartor));
    bootstrap_sample_VAF_vec(i) = goodness_of_fit;
    
    
end
%Find which bootstrap sample gave poor VAF
bootstrap_sample_VAF_vec = bootstrap_sample_VAF_vec(bootstrap_sample_VAF_vec > 70);
regress_coeffs = regress_coeffs(bootstrap_sample_VAF_vec > 70, :);
remaining_after_VAF = length(bootstrap_sample_VAF_vec);

%Return mean of regression coefficients
mean_coeff = mean(regress_coeffs);

%Estimate and return 95% CI
ci=prctile(regress_coeffs,[2.5 97.5]);
CI_coeff = ci(2, :) - mean_coeff;


%Get a goodness of fit as well
% vartor=var(mean(diff_plat_torque));
% varimp=var(mean(diff_plat_torque)-(mean(diff_foot_pos)*mean_coeff(1)+mean(diff_foot_vel)*mean_coeff(2)+mean(diff_foot_acc)*mean_coeff(3)));
% goodness_of_fit = 100*(1-(varimp/vartor));
goodness_of_fit = mean(bootstrap_sample_VAF_vec);


%Very_Temp_Code
hist_fig = figure();
hist(bootstrap_sample_VAF_vec);

if ~(isfile(strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST1.png")))
    saveas(hist_fig, strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST1.png"));
elseif ~(isfile(strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST2.png")))
    saveas(hist_fig, strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST2.png"));
elseif ~(isfile(strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST3.png")))
    saveas(hist_fig, strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST3.png"));
elseif ~(isfile(strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST4.png")))
    saveas(hist_fig, strcat(evalin('caller', 'RESULTS_DIR'), evalin('caller', 'sub_name'), "_VAF_HIST4.png"));
end    

    

                                                       
end

