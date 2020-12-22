function [weight_m, weight_s] = process_weight(weight_p1, weight_p2, weight_p3, weight_p4)

%% ---- Weight Data Points ---- %%

weightf_p1 = mean(weight_p1(:, 99:101), 2);
weightf_p2 = mean(weight_p2(:, 99:101), 2);
weightf_p3 = mean(weight_p3(:, 99:101), 2);
weightf_p4 = mean(weight_p4(:, 99:101), 2);

weightf_1m=mean(weightf_p1);
weightf_2m=mean(weightf_p2);
weightf_3m=mean(weightf_p3);
weightf_4m=mean(weightf_p4);


t_p1=tinv([0.025 0.975],(size(weight_p1,1)));
t_p2=tinv([0.025 0.975],(size(weight_p2,1)));
t_p3=tinv([0.025 0.975],(size(weight_p3,1)));
t_p4=tinv([0.025 0.975],(size(weight_p4,1)));

weightf_1s=t_p1(2)*std(weightf_p1,'omitNaN')/sqrt(size(weight_p1,1));
weightf_2s=t_p2(2)*std(weightf_p2,'omitNaN')/sqrt(size(weight_p2,1));
weightf_3s=t_p3(2)*std(weightf_p3,'omitNaN')/sqrt(size(weight_p3,1));
weightf_4s=t_p4(2)*std(weightf_p4,'omitNaN')/sqrt(size(weight_p4,1));

weight_m = [weightf_1m, weightf_2m,weightf_3m,weightf_4m]';
weight_s = [weightf_1s,weightf_2s,weightf_3s, weightf_4s]';


end