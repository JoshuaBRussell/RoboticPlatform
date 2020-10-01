emg_summary = [emg_rigid;emg_haptic];
data_summary = [Final_value_constrained,boot_imp,goodnessb,emg_summary];
filename = 'datasummary.xlsx';
writematrix(data_summary,filename,'Range','B3:O8');
