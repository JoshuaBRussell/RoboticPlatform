figure
plot(time,p0_plat_torquem,'k','LineWidth',1)
hold on
plot(time,p1_plat_torquem,'r','LineWidth',1)
plot(time,p2_plat_torquem,'g','LineWidth',1)
plot(time,p3_plat_torquem,'b','LineWidth',1)
plot(time,p7_plat_torquem,'k:','LineWidth',1)
hold on
plot(time,p4_plat_torquem,'r:','LineWidth',1)
plot(time,p5_plat_torquem,'g:','LineWidth',1)
plot(time,p6_plat_torquem,'b:','LineWidth',1)
xline(150,'-','20%')
xline(300,'-','40%')
xline(450,'-','60%')
legend('No pert rigid','20% rigid','40% rigid','60% rigid','No pert haptix','20% haptic','40% haptic','60% haptic')
xlabel('time (ms)')
ylabel('Torque (Nm)')
saveas(gcf,'torquecomparison.jpg');
figure