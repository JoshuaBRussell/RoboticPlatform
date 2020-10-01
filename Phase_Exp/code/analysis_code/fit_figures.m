% if(plot_figs==1)
figure

ax1=subplot(4,1,1)

s_time=-20:0.5:130;

for i=1:28 
    plot(s_time,diff_p1_foot_posv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_posvm(60:360),'k');
%ylim([-0.1 0.1])
title('30% rigid');
ax2=subplot(4,1,2)
for i=1:28 

    plot(s_time,diff_p1_foot_velv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_velvm(60:360),'k');
ax3=subplot(4,1,3)
for i=1:28
    plot(s_time,diff_p1_foot_accv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p1_foot_accvm(60:360),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p1_plat_torqueimpm(60:360),'k');
plot(s_time,diff_p1_foot_velvm(60:360)*p1impm(2),'g');
plot(s_time,diff_p1_foot_posvm(60:360)*p1impm(1),'r');
plot(s_time,diff_p1_foot_accvm(60:360)*(p1impm(3)),'b');
plot(s_time,diff_p1_foot_posvm(60:360)*p1impm(1)+diff_p1_foot_velvm(60:360)*p1impm(2)+diff_p1_foot_accvm(60:360)*p1impm(3),'m');
%ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'30% rigid.jpg');
figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:28 
    plot(s_time,diff_p2_foot_posv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p2_foot_posvm(60:360),'k');
%ylim([-0.1 0.1])
title('45% rigid');
ax2=subplot(4,1,2)
for i=1:28 

    plot(s_time,diff_p2_foot_velv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p2_foot_velvm(60:360),'k');
ax3=subplot(4,1,3)
for i=1:28 

    plot(s_time,diff_p2_foot_accv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p2_foot_accvm(60:360),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p2_plat_torqueimpm(60:360),'k');
plot(s_time,diff_p2_foot_velvm(60:360)*p2impm(2),'g');
plot(s_time,diff_p2_foot_posvm(60:360)*p2impm(1),'r');
plot(s_time,diff_p2_foot_accvm(60:360)*(p2impm(3)),'b');
plot(s_time,diff_p2_foot_posvm(60:360)*p2impm(1)+(diff_p2_foot_velvm(60:360)*p2impm(2))+diff_p2_foot_accvm(60:360)*p2impm(3),'m');
%ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');


saveas(gcf,'45% rigid.jpg');figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:28 
    plot(s_time,diff_p3_foot_posv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p3_foot_posvm(60:360),'k');
%ylim([-0.1 0.1])
title('60% rigid');
ax2=subplot(4,1,2)
for i=1:28 

    plot(s_time,diff_p3_foot_velv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p3_foot_velvm(60:360),'k');
ax3=subplot(4,1,3)
for i=1:28 

    plot(s_time,diff_p3_foot_accv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p3_foot_accvm(60:360),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p3_plat_torqueimpm(60:360),'k');
plot(s_time,diff_p3_foot_velvm(60:360)*p3impm(2),'g');
plot(s_time,diff_p3_foot_posvm(60:360)*p3impm(1),'r');
plot(s_time,diff_p3_foot_accvm(60:360)*(p3impm(3)),'b');
plot(s_time,diff_p3_foot_posvm(60:360)*p3impm(1)+(diff_p3_foot_velvm(60:360)*p3impm(2))+diff_p3_foot_accvm(60:360)*p3impm(3),'m');
%ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'60% rigid.jpg');figure
ax1=subplot(4,1,1)
s_time=-20:0.5:130;

for i=1:28 
    plot(s_time,diff_p4_foot_posv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p4_foot_posvm(60:360),'k');
%ylim([-0.1 0.1])
title('15% rigid');
ax2=subplot(4,1,2)
for i=1:28 

    plot(s_time,diff_p4_foot_velv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p4_foot_velvm(60:360),'k');
ax3=subplot(4,1,3)
for i=1:28 

    plot(s_time,diff_p4_foot_accv(i,60:360),'Color',[0.7 0.7 0.7]);
    hold on
end
plot(s_time,diff_p4_foot_accvm(60:360),'k');
ax4=subplot(4,1,4)
hold on
plot(s_time,diff_p4_plat_torqueimpm(60:360),'k');
plot(s_time,diff_p4_foot_velvm(60:360)*p4impm(2),'g');
plot(s_time,diff_p4_foot_posvm(60:360)*p4impm(1),'r');
plot(s_time,diff_p4_foot_accvm(60:360)*(p4impm(3)),'b');
plot(s_time,diff_p4_foot_posvm(60:360)*p4impm(1)+diff_p4_foot_velvm(60:360)*p4impm(2)+diff_p4_foot_accvm(60:360)*p4impm(3),'m');
%ylim([-10 20])
linkaxes([ax1,ax2,ax3,ax4],'x');

saveas(gcf,'15% Rigid.jpg');figure
% end




