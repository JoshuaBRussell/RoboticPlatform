function [mean_coeff, CI_coeff, goodness_of_fit] = bootstrapping_regression(pos_data, vel_data, acc_data, torque_data, perc, A, B)

BOOTSTRAP_COUNT = 100;

regression_coeffs = zeros(BOOTSTRAP_COUNT, 3);
for bootstrapping_counter = 1:BOOTSTRAP_COUNT
    total_num_of_ex = size(pos_data, 1);

    bs_selection_count = round(perc*total_num_of_ex); %Bootstrapping Selection #

    %Randomly sample with replacement
    [pos_data_sample, sample_ind] = datasample(pos_data, bs_selection_count);
    pos_data_sample = pos_data_sample';
    %Also grab the corresponding trial samples from the remaining data
    vel_data_sample = vel_data(sample_ind, :)';
    acc_data_sample = acc_data(sample_ind, :)';
    torque_data_sample = torque_data(sample_ind, :)';
    
    data_matrix = [pos_data_sample(:), vel_data_sample(:), acc_data_sample(:)];
%     data_matrix = [mean(pos_data_sample')', mean(vel_data_sample')', mean(acc_data_sample')']; 
     regression_coeff = lsqlin(data_matrix, torque_data_sample(:), A, B);
%     regression_coeff = lsqlin(data_matrix, mean(torque_data_sample'), A, B);
    regression_coeffs(bootstrapping_counter, :) = regression_coeff';
end

%Return mean of regression coefficients
mean_coeff = mean(regression_coeffs);

%Estimate and return 95% CI
ci=prctile(regression_coeffs,[2.5 97.5]);
CI_coeff = ci(2, :) - mean_coeff;


%Get a goodness of fit as well
vartor=var(mean(torque_data));
varimp=var(mean(torque_data)-(mean(pos_data)*mean_coeff(1)+mean(vel_data)*mean_coeff(2)+mean(acc_data)*mean_coeff(3)));
goodness_of_fit = 100*(1-(varimp/vartor));

end

