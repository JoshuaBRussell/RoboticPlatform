if(plot_figs_constrained==1)

figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:analysis_value-1 
    plot(s_time,diff_p1_foot_pos(i,660:960),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_posm(660:960),'k');
ylim([-0.1 0.1])
title('20% rigid constrained');
ax2=subplot(4,1,2)
for i=1:analysis_value-1 

    plot(s_time,diff_p1_foot_vel(i,660:960),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_velm(660:960),'k');
ax3=subplot(4,1,3)
for i=1:analysis_value-1 

    plot(s_time,diff_p1_foot_acc(i,660:960),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_accm(660:960),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p1_plat_torqueimpm(660:960),'k');
plot(s_time,diff_p1_foot_velm(660:960)*p1impm2(2),'g');
plot(s_time,diff_p1_foot_posm(660:960)*p1impm2(1),'r');
plot(s_time,diff_p1_foot_accm(660:960)*(p1impm2(3)),'b');
plot(s_time,diff_p1_foot_posm(660:960)*p1impm2(1)+diff_p1_foot_velm(660:960)*p1impm2(2)+diff_p1_foot_accm(660:960)*p1impm2(3),'m');
ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'20% rigid_const.jpg');
figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:analysis_value-1 
    plot(s_time,diff_p2_foot_pos(i,960:1260),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p2_foot_posm(960:1260),'k');
ylim([-0.1 0.1])
title('40% rigid constrained');
ax2=subplot(4,1,2)
for i=1:analysis_value-1 

    plot(s_time,diff_p2_foot_vel(i,960:1260),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p2_foot_velm(960:1260),'k');
ax3=subplot(4,1,3)
for i=1:analysis_value-1 

    plot(s_time,diff_p2_foot_acc(i,960:1260),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p2_foot_accm(960:1260),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p2_plat_torqueimpm(960:1260),'k');
plot(s_time,diff_p2_foot_velm(960:1260)*p2impm2(2),'g');
plot(s_time,diff_p2_foot_posm(960:1260)*p2impm2(1),'r');
plot(s_time,diff_p2_foot_accm(960:1260)*(p2impm2(3)),'b');
plot(s_time,diff_p2_foot_posm(960:1260)*p2impm2(1)+(diff_p2_foot_velm(960:1260)*p2impm2(2))+diff_p2_foot_accm(960:1260)*p2impm2(3),'m');
ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');


saveas(gcf,'40% rigid_const.jpg');figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:analysis_value-1 
    plot(s_time,diff_p3_foot_pos(i,1240:1540),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p3_foot_posm(1240:1540),'k');
ylim([-0.1 0.1])
title('60% rigid constrained');
ax2=subplot(4,1,2)
for i=1:analysis_value-1 

    plot(s_time,diff_p3_foot_vel(i,1240:1540),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p3_foot_velm(1240:1540),'k');
ax3=subplot(4,1,3)
for i=1:analysis_value-1 

    plot(s_time,diff_p3_foot_acc(i,1240:1540),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p3_foot_accm(1240:1540),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p3_plat_torqueimpm(1260:1560),'k');
plot(s_time,diff_p3_foot_velm(1260:1560)*p3impm2(2),'g');
plot(s_time,diff_p3_foot_posm(1260:1560)*p3impm2(1),'r');
plot(s_time,diff_p3_foot_accm(1260:1560)*(p3impm2(3)),'b');
plot(s_time,diff_p3_foot_posm(1260:1560)*p3impm2(1)+(diff_p3_foot_velm(1260:1560)*p3impm2(2))+diff_p3_foot_accm(1260:1560)*p3impm2(3),'m');
ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'60% rigid_const.jpg');figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:analysis_value-1 
    plot(s_time,diff_p4_foot_pos(i,680:980),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p4_foot_posm(680:980),'k');
ylim([-0.1 0.1])
title('20% haptic constrained');
ax2=subplot(4,1,2)
for i=1:analysis_value-1 

    plot(s_time,diff_p4_foot_vel(i,680:980),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p4_foot_velm(680:980),'k');
ax3=subplot(4,1,3)
for i=1:analysis_value-1 

    plot(s_time,diff_p4_foot_acc(i,680:980),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p4_foot_accm(680:980),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p4_plat_torqueimpm(680:980),'k');
plot(s_time,diff_p4_foot_velm(680:980)*p4impm2(2),'g');
plot(s_time,diff_p4_foot_posm(680:980)*p4impm2(1),'r');
plot(s_time,diff_p4_foot_accm(680:980)*(p4impm2(3)),'b');
plot(s_time,diff_p4_foot_posm(680:980)*p4impm2(1)+diff_p4_foot_velm(680:980)*p4impm2(2)+diff_p4_foot_accm(680:980)*p4impm2(3),'m');
ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'20% Haptic_const.jpg');figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:analysis_value-1 
    plot(s_time,diff_p5_foot_pos(i,980:1280),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p5_foot_posm(980:1280),'k');
ylim([-0.1 0.1])
title('40% haptic constrained');
ax2=subplot(4,1,2)
for i=1:analysis_value-1 

    plot(s_time,diff_p5_foot_vel(i,980:1280),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p5_foot_velm(980:1280),'k');
ax3=subplot(4,1,3)
for i=1:analysis_value-1 

    plot(s_time,diff_p5_foot_acc(i,980:1280),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p5_foot_accm(980:1280),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p5_plat_torqueimpm(980:1280),'k');
plot(s_time,diff_p5_foot_velm(980:1280)*p5impm2(2),'g');
plot(s_time,diff_p5_foot_posm(980:1280)*p5impm2(1),'r');
plot(s_time,diff_p5_foot_accm(980:1280)*(p5impm2(3)),'b');
plot(s_time,diff_p5_foot_posm(980:1280)*p5impm2(1)+(diff_p5_foot_velm(980:1280)*p5impm2(2))+diff_p5_foot_accm(980:1280)*p5impm2(3),'m');
ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');


saveas(gcf,'40% Haptic_const.jpg');
figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:analysis_value-1 
    plot(s_time,diff_p6_foot_pos(i,1260:1560),'Color',[0.7 0.7 0.7]);
    hold on
end
ylim([-0.1 0.1])
plot(s_time,diff_p6_foot_posm(1260:1560),'k');
title('60% haptic constrained');
ax2=subplot(4,1,2)
for i=1:analysis_value-1 

    plot(s_time,diff_p6_foot_vel(i,1260:1560),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p6_foot_velm(1260:1560),'k');
ax3=subplot(4,1,3)
for i=1:analysis_value-1 

    plot(s_time,diff_p6_foot_acc(i,1260:1560),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p6_foot_accm(1260:1560),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p6_plat_torqueimpm(1260:1560),'k');
plot(s_time,diff_p6_foot_velm(1260:1560)*p6impm2(2),'g');
plot(s_time,diff_p6_foot_posm(1260:1560)*p6impm2(1),'r');
plot(s_time,diff_p6_foot_accm(1260:1560)*(p6impm2(3)),'b');
plot(s_time,diff_p6_foot_posm(1260:1560)*p6impm2(1)+(diff_p6_foot_velm(1260:1560)*p6impm2(2))+diff_p6_foot_accm(1260:1560)*p6impm2(3),'m');
ylim([-10 30])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'60% haptic_const.jpg')
end





