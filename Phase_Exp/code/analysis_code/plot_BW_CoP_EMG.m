%% Meant to be run after walking_study_analysis.m so certian variable can 
%% carry over -- basically they are globals.

%% Normalized (Time) CoP Plot

p0_sample_length = p0_peakend-p0_peakst;
stance_phase_duration = p0_sample_length * (1/2000);


%Preallocates memory
data_total = zeros(size(copr15, 1), max(p0_sample_length)+1);
figure();
for i = 1:size(copr15, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   cop_i = p0_cop_torque(i, p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   weight_i = weight4(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i)));
   data_i = interp1(normalized_time_i', cop_i./weight_i, 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end

plot(0:1/max(p0_sample_length):1,100*mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

formatSpec = 'CoP Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('CoP(cm)');
ylim([-5, 20])
saveas(gcf,strcat(RESULTS_DIR, 'copf_plot.jpg'));


%% Normalized (Time) EMG Plot

%TA

p0_sample_length = p0_peakend-p0_peakst;
stance_phase_duration = p0_sample_length * (1/2000);


%Preallocates memory
data_total = zeros(length(ta15), max(p0_sample_length)+1);
for i = 1:size(data_total, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   data_i = interp1(normalized_time_i', ta_emg(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i))), 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end
figure();
subplot(4,1,1);
plot(0:1/max(p0_sample_length):1, mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

formatSpec = 'EMG Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('TA EMG');

%PL
for i = 1:size(data_total, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   data_i = interp1(normalized_time_i', pl_emg(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i))), 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end
subplot(4,1,2);
plot(0:1/max(p0_sample_length):1, mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('PL EMG');

%SOL
for i = 1:size(data_total, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   data_i = interp1(normalized_time_i', sol_emg(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i))), 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end
subplot(4,1,3);
plot(0:1/max(p0_sample_length):1, mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('SOL EMG');


%GCA
for i = 1:size(data_total, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   data_i = interp1(normalized_time_i', gca_emg(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i))), 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end
subplot(4,1,4);
plot(0:1/max(p0_sample_length):1, mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('GCA EMG');

saveas(gcf, strcat(RESULTS_DIR,'emg_plot.jpg'));


%% Normalized (Time) Weight Plot
for i = 1:size(data_total, 1)
   time_i = 0:(1/2000):stance_phase_duration(p0_raw_data_ind(i));
   normalized_time_i = time_i/stance_phase_duration(p0_raw_data_ind(i));
   
   data_i = interp1(normalized_time_i', weight4(p0_raw_data_ind(i), p0_peakst(p0_raw_data_ind(i)):p0_peakend(p0_raw_data_ind(i))), 0:1/max(p0_sample_length):1);
   
   %plot(normalized_time_i, (p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)));
   %plot(0:1/1441:1, data_i);
   %hold on;
   data_total(i, :) = data_i;
end
figure();
plot(0:1/max(p0_sample_length):1, mean(data_total));
xline(0.18,'-','18%')
xline(0.31,'-','31%')
xline(0.44,'-','44%')
xline(0.57,'-','57%')

formatSpec = 'CoP Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('Weight');

saveas(gcf,strcat(RESULTS_DIR,'weight_plot.jpg'));
