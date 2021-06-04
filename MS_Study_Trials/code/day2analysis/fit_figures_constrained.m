if(plot_figs_constrained==1)

figure
ax1=subplot(4,1,1)
s_time=-50:0.5:150;
for i=1:p1-1 
    plot(s_time,diff_p1_foot_pos(i,1200:1600),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_posm(1200:1600),'k');
ylim([-0.1 0.1])
title('rigid constrained');
ax2=subplot(4,1,2)
for i=1:p1-1 

    plot(s_time,diff_p1_foot_vel(i,1200:1600),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_velm(1200:1600),'k');
ax3=subplot(4,1,3)
for i=1:p1-1 

    plot(s_time,diff_p1_foot_acc(i,1200:1600),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_accm(1200:1600),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p1_plat_torqueimpm(1200:1600),'k');
plot(s_time,diff_p1_foot_velm(1200:1600)*p1impm2(2),'g');
plot(s_time,diff_p1_foot_posm(1200:1600)*p1impm2(1),'r');
plot(s_time,diff_p1_foot_accm(1200:1600)*(p1impm2(3)),'b');
plot(s_time,diff_p1_foot_posm(1200:1600)*p1impm2(1)+diff_p1_foot_velm(1200:1600)*p1impm2(2)+diff_p1_foot_accm(1200:1600)*p1impm2(3),'m');
ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'20% rigid_const.jpg');

end





