function [mean_coeff, CI_coeff, goodness_of_fit] = regress_after_bootstrap(diff_foot_pos, diff_foot_vel, ...
                                                           diff_foot_acc, diff_plat_torque, ...
                                                           A, B)

                                                       
regress_coeffs = zeros(size(diff_foot_pos, 1), 3);
for i = 1:size(diff_foot_pos, 1)
    data_matrix = [diff_foot_pos(i, :)', ...
                   diff_foot_vel(i, :)', ...
                   diff_foot_acc(i, :)'];
               
    torque_vec =  diff_plat_torque(i, :)';
    regression_coeff = lsqlin(data_matrix, torque_vec, A, B);
    regress_coeffs(i, :) = regression_coeff;
end        

%Return mean of regression coefficients
mean_coeff = mean(regress_coeffs);

%Estimate and return 95% CI
ci=prctile(regress_coeffs,[2.5 97.5]);
CI_coeff = ci(2, :) - mean_coeff;


%Get a goodness of fit as well
vartor=var(mean(diff_plat_torque));
varimp=var(mean(diff_plat_torque)-(mean(diff_foot_pos)*mean_coeff(1)+mean(diff_foot_vel)*mean_coeff(2)+mean(diff_foot_acc)*mean_coeff(3)));
goodness_of_fit = 100*(1-(varimp/vartor));

                                                       
end

