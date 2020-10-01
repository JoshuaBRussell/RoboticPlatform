figure
subplot(4,1,1)

plot(time,ta_emgm,'b');
hold on
plot(time,ta_emghm,'r');
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
title('TA EMG');
ylabel('TA (%mvc)');
legend('Rigid','Haptic')
subplot(4,1,2)

plot(time,sol_emgm,'b');
hold on
plot(time,sol_emghm,'r');
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
title('SOL EMG');
ylabel('SOL (%mvc)');
subplot(4,1,3)

plot(time,pl_emgm,'b');
hold on
plot(time,pl_emghm,'r');
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
title('PL EMG');
ylabel('PL (%mvc)');
subplot(4,1,4)

plot(time,gca_emgm,'b');
hold on
plot(time,gca_emghm,'r');
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
title('GCA EMG');
ylabel('GCA (%mvc)');
saveas(gcf,'emg_plot.jpg');

formatSpec = 'Muscle Activation %s';
str = sprintf(formatSpec,sub_name);
sgtitle(str)

figure

plot(time,weight4m,'b');
hold on
plot(time,weight7m,'r');
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
formatSpec = 'Weight Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid','Haptic');
xlabel('time(ms)');
ylabel('Weight(N)');
saveas(gcf,'weight_plot.jpg');
%%
for i=1:size(ta_emg)
ta_20(i)=mean(ta_emg(i,650:750));
ta_40(i)=mean(ta_emg(i,950:1050));
ta_60(i)=mean(ta_emg(i,1230:1330));
end
ta_20m=mean(ta_20);
ta_40m=mean(ta_40);
ta_60m=mean(ta_60);

for i=1:size(sol_emg)
sol_20(i)=mean(sol_emg(i,650:750));
sol_40(i)=mean(sol_emg(i,950:1050));
sol_60(i)=mean(sol_emg(i,1230:1330));
end
sol_20m=mean(sol_20);
sol_40m=mean(sol_40);
sol_60m=mean(sol_60);

for i=1:size(pl_emg)
pl_20(i)=mean(pl_emg(i,650:750));
pl_40(i)=mean(pl_emg(i,950:1050));
pl_60(i)=mean(pl_emg(i,1230:1330));
end
pl_20m=mean(pl_20);
pl_40m=mean(pl_40);
pl_60m=mean(pl_60);

for i=1:size(gca_emg)
gca_20(i)=mean(gca_emg(i,650:750));
gca_40(i)=mean(gca_emg(i,950:1050));
gca_60(i)=mean(gca_emg(i,1230:1330));
end
gca_20m=mean(gca_20);
gca_40m=mean(gca_40);
gca_60m=mean(gca_60);

for i=1:size(ta_emgh)
ta_20h(i)=mean(ta_emgh(i,650:750));
ta_40h(i)=mean(ta_emgh(i,950:1050));
ta_60h(i)=mean(ta_emgh(i,1230:1330));
end
ta_20hm=mean(ta_20h);
ta_40hm=mean(ta_40h);
ta_60hm=mean(ta_60h);

for i=1:size(sol_emgh)
sol_20h(i)=mean(sol_emgh(i,650:750));
sol_40h(i)=mean(sol_emgh(i,950:1050));
sol_60h(i)=mean(sol_emgh(i,1230:1330));
end
sol_20hm=mean(sol_20h);
sol_40hm=mean(sol_40h);
sol_60hm=mean(sol_60h);

for i=1:size(pl_emgh)
pl_20h(i)=mean(pl_emgh(i,650:750));
pl_40h(i)=mean(pl_emgh(i,950:1050));
pl_60h(i)=mean(pl_emgh(i,1230:1330));
end
pl_20hm=mean(pl_20h);
pl_40hm=mean(pl_40h);
pl_60hm=mean(pl_60h);

for i=1:size(gca_emgh)
gca_20h(i)=mean(gca_emgh(i,650:750));
gca_40h(i)=mean(gca_emgh(i,950:1050));
gca_60h(i)=mean(gca_emgh(i,1230:1330));
end
gca_20hm=mean(gca_20h);
gca_40hm=mean(gca_40h);
gca_60hm=mean(gca_60h);

for i=1:size(weight4)
weight_20(i)=mean(weight4(i,650:750));
weight_40(i)=mean(weight4(i,950:1050));
weight_60(i)=mean(weight4(i,1230:1330));
end
weight_20m=mean(weight_20);
weight_40m=mean(weight_40);
weight_60m=mean(weight_60);
t=tinv([0.025 0.975],(size(weight4,1)))
weight_20s=t(2)*std(weight_20)/sqrt(size(weight4,1));
weight_40s=t(2)*std(weight_40)/sqrt(size(weight4,1));
weight_60s=t(2)*std(weight_60)/sqrt(size(weight4,1));

for i=1:size(weight7)
weighth_20(i)=mean(weight7(i,650:750));
weighth_40(i)=mean(weight7(i,950:1050));
weighth_60(i)=mean(weight7(i,1230:1330));
end
weighth_20m=mean(weighth_20);
weighth_40m=mean(weighth_40);
weighth_60m=mean(weighth_60);
t=tinv([0.025 0.975],(size(weight7,1)))
weighth_20s=t(2)*std(weighth_20)/sqrt(size(weight7,1));
weighth_40s=t(2)*std(weighth_40)/sqrt(size(weight7,1));
weighth_60s=t(2)*std(weighth_60)/sqrt(size(weight7,1));



for i=1:size(cop7,1)
   for j=1:length(cop7)
copf7(i,j)=cop7(i,j)*100/weight7(i,j);
copf4(i,j)=cop4(i,j)*100/weight4(i,j);
   end
end
copf4m=trimmean(copf4,30);
copf7m=trimmean(copf7,30);
for i=1:size(copf4)
copf_20(i)=trimmean(copf4(i,650:750),30);
copf_40(i)=trimmean(copf4(i,950:1050),30);
copf_60(i)=trimmean(copf4(i,1230:1330),30);
end
copf_20m=trimmean(copf_20,30);
copf_40m=trimmean(copf_40,30);
copf_60m=trimmean(copf_60,30);
t=tinv([0.025 0.975],(size(copf4,1)))
copf_20s=t(2)*std(copf_20)/sqrt(size(copf4,1));
copf_40s=t(2)*std(copf_40)/sqrt(size(copf4,1));
copf_60s=t(2)*std(copf_60)/sqrt(size(copf4,1));

for i=1:size(copf7)
copfh_20(i)=trimmean(copf7(i,650:750),30);
copfh_40(i)=trimmean(copf7(i,950:1050),30);
copfh_60(i)=trimmean(copf7(i,1230:1330),30);
end
copfh_20m=trimmean(copfh_20,30);
copfh_40m=trimmean(copfh_40,30);
copfh_60m=trimmean(copfh_60,30);
t=tinv([0.025 0.975],(size(copf7,1)))
copfh_20s=t(2)*std(copfh_20)/sqrt(size(copf7,1));
copfh_40s=t(2)*std(copfh_40)/sqrt(size(copf7,1));
copfh_60s=t(2)*std(copfh_60)/sqrt(size(copf7,1));
figure
plot(time,copf4m,'b');
hold on
plot(time,copf7m,'r');
axis([0 700 -5 25])
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
formatSpec = 'CoP Profile %s';
str = sprintf(formatSpec,sub_name);
title(str)
legend('Rigid','Haptic');
xlabel('time(ms)');
ylabel('CoP(cm)');
saveas(gcf,'copf_plot.jpg');


emg_rigid=[ta_20m, sol_20m, pl_20m, gca_20m,weight_20m,copf_20m;ta_40m,sol_40m,pl_40m,gca_40m,weight_40m,copf_40m;ta_60m,sol_60m,pl_60m,gca_60m,weight_60m,copf_60m;]
emg_haptic=[ta_20hm, sol_20hm, pl_20hm, gca_20hm,weighth_20m,copfh_20m;ta_40hm,sol_40hm,pl_40hm,gca_40hm,weighth_40m,copfh_40m;ta_60hm,sol_60hm,pl_60hm,gca_60hm,weighth_60m,copfh_60m;]
weight_rigid=[weight_20m, weight_20s,copf_20m,copf_20s;weight_40m, weight_40s,copf_40m,copf_40s;weight_60m, weight_60s,copf_60m,copf_60s]
weighth_rigid=[weighth_20m, weighth_20s,copfh_20m,copfh_20s;weighth_40m, weighth_40s,copfh_40m,copfh_40s;weighth_60m, weighth_60s,copfh_60m,copfh_60s]
save('emg.mat','emg_rigid');
save('emg_haptic.mat','emg_haptic');
w=[weight_rigid;weighth_rigid];

