function [mean_coeff, CI_coeff, goodness_of_fit] = regress_after_bootstrap(diff_foot_pos, diff_foot_vel, ...
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
                                                       
end

