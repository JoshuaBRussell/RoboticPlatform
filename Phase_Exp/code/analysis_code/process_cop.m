function [cop_m, cop_s] = process_cop(cop_p1, cop_p2, cop_p3, cop_p4)

%% ---- CoP Data Points ---- %%

copf_p1 = mean(cop_p1(:, 99:101), 2);
copf_p2 = mean(cop_p2(:, 99:101), 2);
copf_p3 = mean(cop_p3(:, 99:101), 2);
copf_p4 = mean(cop_p4(:, 99:101), 2);

copf_1m=mean(copf_p1);
copf_2m=mean(copf_p2);
copf_3m=mean(copf_p3);
copf_4m=mean(copf_p4);

t_p1=tinv([0.025 0.975],(size(copf_p1,1)));
t_p2=tinv([0.025 0.975],(size(copf_p2,1)));
t_p3=tinv([0.025 0.975],(size(copf_p3,1)));
t_p4=tinv([0.025 0.975],(size(copf_p4,1)));


copf_1s=t_p1(2)*std(copf_p1,'omitNaN')/sqrt(size(copf_p1,1));
copf_2s=t_p2(2)*std(copf_p2,'omitNaN')/sqrt(size(copf_p2,1));
copf_3s=t_p3(2)*std(copf_p3,'omitNaN')/sqrt(size(copf_p3,1));
copf_4s=t_p4(2)*std(copf_p4,'omitNaN')/sqrt(size(copf_p4,1));

cop_m = [copf_1m, copf_2m,copf_3m,copf_4m]';
cop_s = [copf_1s,copf_2s,copf_3s, copf_4s]';