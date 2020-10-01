
%% ---- Weight Data Points ---- %%
for i=1:size(copf4,1) %number of pertubrations at that percentage
weightf_15(i)=mean(weightr15(i,99:101));
weightf_30(i)=mean(weightr30(i,99:101));
weightf_45(i)=mean(weightr45(i,99:101));
weightf_60(i)=mean(weightr60(i,99:101));
end

weightf_15m=mean(weightf_15);
weightf_30m=mean(weightf_30);
weightf_45m=mean(weightf_45);
weightf_60m=mean(weightf_60);

t=tinv([0.025 0.975],(size(copf4,1)))
weightf_15s=t(2)*std(weightf_15,'omitNaN')/sqrt(size(weightr15,1));
weightf_30s=t(2)*std(weightf_30,'omitNaN')/sqrt(size(weightr15,1));
weightf_45s=t(2)*std(weightf_45,'omitNaN')/sqrt(size(weightr15,1));
weightf_60s=t(2)*std(weightf_60,'omitNaN')/sqrt(size(weightr15,1));

weight_m = [weightf_15m, weightf_30m,weightf_45m,weightf_60m]';
weight_s = [weightf_15s, weightf_30s,weightf_45s,weightf_60s]';

%% ---- CoP Data Points ---- %%

for i=1:size(cop4,1) %number of pertubrations at that percentage
copf_15(i)=mean(copr15(i,99:101));
copf_30(i)=mean(copr30(i,99:101));
copf_45(i)=mean(copr45(i,99:101));
copf_60(i)=mean(copr60(i,99:101));
end

copf_15m=mean(copf_15);
copf_30m=mean(copf_30);
copf_45m=mean(copf_45);
copf_60m=mean(copf_60);

t=tinv([0.025 0.975],(size(copf4,1)));
copf_15s=t(2)*std(copf_15,'omitNaN')/sqrt(size(copf4,1));
copf_30s=t(2)*std(copf_30,'omitNaN')/sqrt(size(copf4,1));
copf_45s=t(2)*std(copf_45,'omitNaN')/sqrt(size(copf4,1));
copf_60s=t(2)*std(copf_60,'omitNaN')/sqrt(size(copf4,1));

cop_m = [copf_15m, copf_30m,copf_45m,copf_60m]';
cop_s = [copf_15s, copf_30s,copf_45s,copf_60s]';

%% Normalized (Time) CoP Plot

p0_sample_length = p0_peakend-p0_peakst;
stance_phase_duration = p0_sample_length * (1/2000);


%Preallocates memory
data_total = zeros(40, max(p0_sample_length)+1);
figure();
for i = 1:length(stance_phase_duration)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', p0_cop_torque(i, p0_peakst(i):p0_peakend(i))./weight4(i, p0_peakst(i):p0_peakend(i)), 0:1/max(p0_sample_length):1);
   
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
hold on;
scatter(0.18, 100*copf_15m); hold on;
scatter(0.31, 100*copf_30m); hold on;
scatter(0.44, 100*copf_45m); hold on;
scatter(0.57, 100*copf_60m); hold on;

formatSpec = 'CoP Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('CoP(cm)');
saveas(gcf,'copf_plot.jpg');


%% Normalized (Time) EMG Plot

%TA

p0_sample_length = p0_peakend-p0_peakst;
stance_phase_duration = p0_sample_length * (1/2000);


%Preallocates memory
data_total = zeros(40, max(p0_sample_length)+1);
for i = 1:length(stance_phase_duration)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', ta_emg(i, p0_peakst(i):p0_peakend(i)), 0:1/max(p0_sample_length):1);
   
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
hold on;
scatter(0.18, ta15m); hold on;
scatter(0.31, ta30m); hold on;
scatter(0.44, ta45m); hold on;
scatter(0.57, ta60m); hold on;

formatSpec = 'EMG Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('TA EMG');

%PL
for i = 1:length(stance_phase_duration)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', pl_emg(i, p0_peakst(i):p0_peakend(i)), 0:1/max(p0_sample_length):1);
   
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
hold on;
scatter(0.18, pl15m); hold on;
scatter(0.31, pl30m); hold on;
scatter(0.44, pl45m); hold on;
scatter(0.57, pl60m); hold on;


str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('PL EMG');

%SOL
for i = 1:length(stance_phase_duration)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', sol_emg(i, p0_peakst(i):p0_peakend(i)), 0:1/max(p0_sample_length):1);
   
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
hold on;
scatter(0.18, sol15m); hold on;
scatter(0.31, sol30m); hold on;
scatter(0.44, sol45m); hold on;
scatter(0.57, sol60m); hold on;


str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('SOL EMG');


%GCA
for i = 1:length(stance_phase_duration)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', gca_emg(i, p0_peakst(i):p0_peakend(i)), 0:1/max(p0_sample_length):1);
   
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
hold on;
scatter(0.18, gca15m); hold on;
scatter(0.31, gca30m); hold on;
scatter(0.44, gca45m); hold on;
scatter(0.57, gca60m); hold on;


str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('GCA EMG');

saveas(gcf,'emg_plot.jpg');


%% Normalized (Time) Weight Plot
for i = 1:length(stance_phase_duration)
   time_i = 0:(1/2000):stance_phase_duration(i);
   normalized_time_i = time_i/stance_phase_duration(i);
   
   data_i = interp1(normalized_time_i', weight4(i, p0_peakst(i):p0_peakend(i)), 0:1/max(p0_sample_length):1);
   
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
hold on;
scatter(0.18, weightf_15m); hold on;
scatter(0.31, weightf_30m); hold on;
scatter(0.44, weightf_45m); hold on;
scatter(0.57, weightf_60m); hold on;

formatSpec = 'CoP Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('Gait Cycle (%)');
ylabel('Weight');

saveas(gcf,'weight_plot.jpg');
