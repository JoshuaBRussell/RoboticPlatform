
%----Weight Data Points

%The data points start 100 samples (50ms) before the perturbation. This
%gives +/- 50 samples (25ms) around the perturbation point.
%mean() averages all rows together. The transposes allow us to average all points about a
%single perturbation in a single times series together first, THEN average
%over the averages.
weightr15_mean_across_time = mean(weightr15(:, 50:150)')';
weightr30_mean_across_time = mean(weightr30(:, 50:150)')';
weightr45_mean_across_time = mean(weightr45(:, 50:150)')';
weightr60_mean_across_time = mean(weightr60(:, 50:150)')';


weightr15_mean_across_samples = mean(weightr15_mean_across_time);
weightr30_mean_across_samples = mean(weightr30_mean_across_time);
weightr45_mean_across_samples = mean(weightr45_mean_across_time);
weightr60_mean_across_samples = mean(weightr60_mean_across_time);

weightm = [weightr15_mean_across_samples, weightr30_mean_across_samples, ... 
           weightr45_mean_across_samples, weightr60_mean_across_samples]';

t=tinv([0.025 0.975],(size(weightr15_mean_across_time,1)));
weightr15_s=t(2)*std(weightr15_mean_across_time,'omitNaN')/sqrt(size(weight4,1));
weightr30_s=t(2)*std(weightr30_mean_across_time,'omitNaN')/sqrt(size(weight4,1));
weightr45_s=t(2)*std(weightr45_mean_across_time,'omitNaN')/sqrt(size(weight4,1));
weightr60_s=t(2)*std(weightr60_mean_across_time,'omitNaN')/sqrt(size(weight4,1));

weight_s = [weightr15_s, weightr30_s, weightr45_s, weightr60_s]';

%----CoP Data Points
%The data points start 100 samples (50ms) before the perturbation. This
%gives +/- 50 samples (25ms) around the perturbation point.
%mean() averages all rows together. The transposes allow us to average all points about a
%single perturbation in a single times series together first, THEN average
%over the averages.
copr15_mean_across_time = mean(copr15(:, 50:150)')';
copr30_mean_across_time = mean(copr30(:, 50:150)')';
copr45_mean_across_time = mean(copr45(:, 50:150)')';
copr60_mean_across_time = mean(copr60(:, 50:150)')';


copr15_mean_across_samples = mean(copr15_mean_across_time);
copr30_mean_across_samples = mean(copr30_mean_across_time);
copr45_mean_across_samples = mean(copr45_mean_across_time);
copr60_mean_across_samples = mean(copr60_mean_across_time);

copm = [copr15_mean_across_samples, copr30_mean_across_samples, ... 
           copr45_mean_across_samples, copr60_mean_across_samples]';

t=tinv([0.025 0.975],(size(copr15_mean_across_time,1)));
copr15_s=t(2)*std(copr15_mean_across_time,'omitNaN')/sqrt(size(cop4,1));
copr30_s=t(2)*std(copr30_mean_across_time,'omitNaN')/sqrt(size(cop4,1));
copr45_s=t(2)*std(copr45_mean_across_time,'omitNaN')/sqrt(size(cop4,1));
copr60_s=t(2)*std(copr60_mean_across_time,'omitNaN')/sqrt(size(cop4,1));

cop_s = [copr15_s, copr30_s, copr45_s, copr60_s]';

%% ----CoP Plot

for i=1:size(cop4,1)
   for j=1:length(cop4)
copf4(i,j)=p0_plat_torque(i,j)*100/weight4(i,j);
   end
end
copf4m=trimmean(copf4,30);

figure
plot(time,copf4m,'b');
axis([0 700 -5 25])
xline(135,'-','18%')
xline(232,'-','31%')
xline(330,'-','44%')
xline(428,'-','57%')
formatSpec = 'CoP Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid');
xlabel('time(ms)');
ylabel('CoP(cm)');
saveas(gcf,'copf_plot.jpg');


