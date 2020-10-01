ideal=floor(num_pert*0.7);

for i=1:analysis_value-1
if(isnan(p1impg(i))==1)
    p1impg(i)=-1000;
end
if(isnan(p2impg(i))==1)
    p2impg(i)=-1000;
end
if(isnan(p3impg(i))==1)
    p3impg(i)=-1000;
end
if(isnan(p4impg(i))==1)
    p4impg(i)=-1000;
end

end
%%


for i=1:analysis_value-1
   p1imp(i,4)=p1impg(i);
   p2imp(i,4)=p2impg(i);
   p3imp(i,4)=p3impg(i);
   p4imp(i,4)=p4impg(i);
  
end

[~,idx] = sort(p1imp(:,4),'descend'); % sort just the first column
p1imps = p1imp(idx,:);
diff_p1_plat_torqueimps=diff_p1_plat_torqueimp(idx,:);
diff_p1_foot_poss=diff_p1_foot_pos(idx,:);
diff_p1_foot_vels=diff_p1_foot_vel(idx,:);
diff_p1_foot_accs=diff_p1_foot_acc(idx,:);

diff_p1_plat_torqueimpv=diff_p1_plat_torqueimps(1:ideal,:);
diff_p1_foot_posv=diff_p1_foot_poss(1:ideal,:);
diff_p1_foot_velv=diff_p1_foot_vels(1:ideal,:);
diff_p1_foot_accv=diff_p1_foot_accs(1:ideal,:);



[~,idx] = sort(p2imp(:,4),'descend'); % sort just the first column
p2imps = p2imp(idx,:);
diff_p2_plat_torqueimps=diff_p2_plat_torqueimp(idx,:);
diff_p2_foot_poss=diff_p2_foot_pos(idx,:);
diff_p2_foot_vels=diff_p2_foot_vel(idx,:);
diff_p2_foot_accs=diff_p2_foot_acc(idx,:);

diff_p2_plat_torqueimpv=diff_p2_plat_torqueimps(1:ideal,:);
diff_p2_foot_posv=diff_p2_foot_poss(1:ideal,:);
diff_p2_foot_velv=diff_p2_foot_vels(1:ideal,:);
diff_p2_foot_accv=diff_p2_foot_accs(1:ideal,:);

[~,idx] = sort(p3imp(:,4),'descend'); % sort just the first column
p3imps = p3imp(idx,:);
diff_p3_plat_torqueimps=diff_p3_plat_torqueimp(idx,:);
diff_p3_foot_poss=diff_p3_foot_pos(idx,:);
diff_p3_foot_vels=diff_p3_foot_vel(idx,:);
diff_p3_foot_accs=diff_p3_foot_acc(idx,:);

diff_p3_plat_torqueimpv=diff_p3_plat_torqueimps(1:ideal,:);
diff_p3_foot_posv=diff_p3_foot_poss(1:ideal,:);
diff_p3_foot_velv=diff_p3_foot_vels(1:ideal,:);
diff_p3_foot_accv=diff_p3_foot_accs(1:ideal,:);

[~,idx] = sort(p4imp(:,4),'descend'); % sort just the first column
p4imps = p4imp(idx,:);
diff_p4_plat_torqueimps=diff_p4_plat_torqueimp(idx,:);
diff_p4_foot_poss=diff_p4_foot_pos(idx,:);
diff_p4_foot_vels=diff_p4_foot_vel(idx,:);
diff_p4_foot_accs=diff_p4_foot_acc(idx,:);

diff_p4_plat_torqueimpv=diff_p4_plat_torqueimps(1:ideal,:);
diff_p4_foot_posv=diff_p4_foot_poss(1:ideal,:);
diff_p4_foot_velv=diff_p4_foot_vels(1:ideal,:);
diff_p4_foot_accv=diff_p4_foot_accs(1:ideal,:);



for i=1:analysis_value-1
   if(p1imps(i,4)>80)
   p1_fin=i;
   end
   if(p2imps(i,4)>80)
   p2_fin=i;
   end
   if(p3imps(i,4)>80)
   p3_fin=i;
   end
   if(p4imps(i,4)>80)
   p4_fin=i;
   end
   
end

% p1_fin2=max(p1_fin,ideal);
% p2_fin2=max(p2_fin,ideal);
% p3_fin2=max(p3_fin,ideal);
% p4_fin2=max(p4_fin,ideal);
% p5_fin2=max(p5_fin,ideal);
% p6_fin2=max(p6_fin,ideal);
%  p_fin=[p1_fin p2_fin p3_fin p4_fin p5_fin p6_fin]';
%  p_fin2=[p1_fin2 p2_fin2 p3_fin2 p4_fin2 p5_fin2 p6_fin2]';
%%

% indp2(1,:)=mean(p1imps(1:p1_fin2,:));
% indp2(2,:)=mean(p2imps(1:p2_fin2,:));
% indp2(3,:)=mean(p3imps(1:p3_fin2,:));
% indp2(4,:)=mean(p4imps(1:p4_fin2,:));
% indp2(5,:)=mean(p5imps(1:p5_fin2,:));
% indp2(6,:)=mean(p6imps(1:p6_fin2,:));



% indp(1,:)=mean(p1imps(1:p1_fin,:));
% indp(2,:)=mean(p2imps(1:p2_fin,:));
% indp(3,:)=mean(p3imps(1:p3_fin,:));
% indp(4,:)=mean(p4imps(1:p4_fin,:));
% indp(5,:)=mean(p5imps(1:p5_fin,:));
% indp(6,:)=mean(p6imps(1:p6_fin,:));

% sindp2(1,:)=std(p1imps(1:p1_fin2,:))/sqrt(p1_fin2);
% sindp2(2,:)=std(p2imps(1:p2_fin2,:))/sqrt(p2_fin2);
% sindp2(3,:)=std(p3imps(1:p3_fin2,:))/sqrt(p3_fin2);
% sindp2(4,:)=std(p4imps(1:p4_fin2,:))/sqrt(p4_fin2);
% sindp2(5,:)=std(p5imps(1:p5_fin2,:))/sqrt(p5_fin2);
% sindp2(6,:)=std(p6imps(1:p6_fin2,:))/sqrt(p6_fin2);
% 
% sindp(1,:)=std(p1imps(1:p1_fin,:))/sqrt(p1_fin);
% sindp(2,:)=std(p2imps(1:p2_fin,:))/sqrt(p2_fin);
% sindp(3,:)=std(p3imps(1:p3_fin,:))/sqrt(p3_fin);
% sindp(4,:)=std(p4imps(1:p4_fin,:))/sqrt(p4_fin);
% sindp(5,:)=std(p5imps(1:p5_fin,:))/sqrt(p5_fin);
% sindp(6,:)=std(p6imps(1:p6_fin,:))/sqrt(p6_fin);
% 
% 
% tfin2(1,:)=tinv([0.025 0.975],(p1_fin2-1));
% tfin2(2,:)=tinv([0.025 0.975],(p2_fin2-1));
% tfin2(3,:)=tinv([0.025 0.975],(p3_fin2-1));
% tfin2(4,:)=tinv([0.025 0.975],(p4_fin2-1));

% tfin(1,:)=tinv([0.025 0.975],(p1_fin-1));
% tfin(2,:)=tinv([0.025 0.975],(p2_fin-1));
% tfin(3,:)=tinv([0.025 0.975],(p3_fin-1));
% tfin(4,:)=tinv([0.025 0.975],(p4_fin-1));


% for i=1:4
%     for j=1:4
%         
%         cindp(i,j)=(tfin(i,2)*sindp(i,j));
% %         cindp2(i,j)=(tfin2(i,2)*sindp2(i,j));
%         
%     end
%     
%     
% end


% fin_ind=[indp cindp(:,1) p_fin indp2 cindp2(:,1) p_fin2]