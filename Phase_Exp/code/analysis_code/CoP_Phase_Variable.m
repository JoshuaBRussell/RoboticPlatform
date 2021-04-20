% ---- Look at CoP as Phase Variable ---- %

p0_sample_length = p0_peakend-p0_peakst;
stance_phase_duration = p0_sample_length * (1/2000);

figure();
for i = 1:size(copr15, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   cop_trq_i = p0_cop_torque(i, p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   weight_i = weight4(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   cop_i = cop_trq_i./weight_i;
   ankle_angle_i = p0_foot_pos(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   [ankle_vel_i, ankle_acc_i] = get_derivatives(ankle_angle_i, 1/2000);
   start_cop_vec(i) = cop_i(p0_peakst(p0_raw_data_ind(i)));
   end_cop_vec(i) = cop_i(p0_peakend(p0_raw_data_ind(i)) - p0_peakst(p0_raw_data_ind(i)));
   plot(cop_i,ankle_angle_i); hold on;
end

relative_TO_cop_mean = mean(end_cop_vec - start_cop_vec);



%% Mean Relative CoP
ref_phaseline = [0:1:999]*(relative_TO_cop_mean/999);

data_total = zeros(size(copr15, 1), length(ref_phaseline));
figure();
for i = 1:size(copr15, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/1000;
   
   foot_pos_i = p0_foot_pos(i, p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   data_i = interp1(1:length(foot_pos_i), foot_pos_i, ref_phaseline');
   

   data_total(i, :) = data_i;
end

plot(data_total')

%% Individual Relative CoP
ref_phaseline = [0:1:999]*(relative_TO_cop_mean/999);

data_total = zeros(size(copr15, 1), length(ref_phaseline));
figure();
for i = 1:size(copr15, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/1000;
   
   foot_pos_i = p0_foot_pos(i, p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   data_i = interp1(1:length(foot_pos_i), foot_pos_i, [0:1:999]');
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end

plot(data_total')
