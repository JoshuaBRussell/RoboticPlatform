function [] = plot_regression_results(pos_data,vel_data, acc_data, observed_data, regress_coeff, title_str)
%plot_regression_results 
%Plots the position, velocity, and acceleration time series data used
%for regression. Also, plots the regression model output along side the
%original output data.


LOWER_SAMPLE_PLOT_LIMIT = 60;
UPPER_SAMPLE_PLOT_LIMIT = 360;
figure();

%% Position Time Series Data
ax1=subplot(4,1,1)

s_time=-20:0.5:130;

for i=1:28 
    plot(s_time, pos_data(i,LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'Color',[0.7 0.7 0.7]);
    hold on
end
pos_data_mean=mean(pos_data);
plot(s_time,pos_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'k');
ylim([-0.1 0.1])
title(title_str);

%% Velocity Time Series Data
ax2=subplot(4,1,2)
for i=1:28 

    plot(s_time, vel_data(i,LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'Color',[0.7 0.7 0.7]);
    hold on
end
vel_data_mean = mean(vel_data);
plot(s_time,vel_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'k');
ylim([-1, 1]);

%% Acceleration Time Series Data
ax3=subplot(4,1,3)
for i=1:28
    plot(s_time, acc_data(i,LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'Color',[0.7 0.7 0.7]);
    hold on
end
acc_data_mean = mean(acc_data);
plot(s_time, acc_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'k');
ylim([-50 50]);

%% Regresion Model Time Series Data
ax4=subplot(4,1,4)
hold on
observed_data_mean = mean(observed_data);
plot(s_time, observed_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT),'k');

plot(s_time, pos_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT)*regress_coeff(1),'r');
plot(s_time, vel_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT)*regress_coeff(2),'g');
plot(s_time, acc_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT)*regress_coeff(3),'b');

plot(s_time,pos_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT)*regress_coeff(1)+vel_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT)*regress_coeff(2)+acc_data_mean(LOWER_SAMPLE_PLOT_LIMIT:UPPER_SAMPLE_PLOT_LIMIT)*regress_coeff(3),'m');
ylim([-5 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,strcat(title_str, '_rigid.jpg'));


end

