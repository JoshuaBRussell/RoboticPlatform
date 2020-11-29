
ideal=floor(NUM_PERT*0.7);

diff_p1_plat_torqueimpvm=mean(diff_p1_plat_torqueimpv);
diff_p1_foot_posvm=mean(diff_p1_foot_posv);
diff_p1_foot_velvm=mean(diff_p1_foot_velv);
diff_p1_foot_accvm=mean(diff_p1_foot_accv);

diff_p2_plat_torqueimpvm=mean(diff_p2_plat_torqueimpv);
diff_p2_foot_posvm=mean(diff_p2_foot_posv);
diff_p2_foot_velvm=mean(diff_p2_foot_velv);
diff_p2_foot_accvm=mean(diff_p2_foot_accv);

diff_p3_plat_torqueimpvm=mean(diff_p3_plat_torqueimpv);
diff_p3_foot_posvm=mean(diff_p3_foot_posv);
diff_p3_foot_velvm=mean(diff_p3_foot_velv);
diff_p3_foot_accvm=mean(diff_p3_foot_accv);

diff_p4_plat_torqueimpvm=mean(diff_p4_plat_torqueimpv);
diff_p4_foot_posvm=mean(diff_p4_foot_posv);
diff_p4_foot_velvm=mean(diff_p4_foot_velv);
diff_p4_foot_accvm=mean(diff_p4_foot_accv);


for i=1:BOOTSTRAPPING_LOOPS
    unique1=0;
    unique2=0;
    unique3=0;
    unique4=0;
   
    
    while(unique1<60)
        [bunch1,idx1]=datasample(diff_p1_plat_torqueimpv,ideal);
        numb=length(unique(idx1));
        unique1=100*(numb/ideal);
    end
    unique1=0;
    while(unique2<60)
        [bunch2,idx2]=datasample(diff_p2_plat_torqueimpv,ideal);
        numb=length(unique(idx2));
        unique2=100*(numb/ideal);
    end
    unique2=0;
    while(unique3<60)
        [bunch3,idx3]=datasample(diff_p3_plat_torqueimpv,ideal);
        numb=length(unique(idx3));
        unique3=100*(numb/ideal);
    end
    unique3=0;
    while(unique4<60)
        [bunch4,idx4]=datasample(diff_p4_plat_torqueimpv,ideal);
        numb=length(unique(idx4));
        unique4=100*(numb/ideal);
    end
    unique4=0;
   
    
    
    for j=1:ideal
        p_bunch1(j,:)=diff_p1_foot_posv(idx1(j),:);
        v_bunch1(j,:)=diff_p1_foot_velv(idx1(j),:);
        a_bunch1(j,:)=diff_p1_foot_accv(idx1(j),:);
        p_bunch2(j,:)=diff_p2_foot_posv(idx2(j),:);
        v_bunch2(j,:)=diff_p2_foot_velv(idx2(j),:);
        a_bunch2(j,:)=diff_p2_foot_accv(idx2(j),:);
        p_bunch3(j,:)=diff_p3_foot_posv(idx3(j),:);
        v_bunch3(j,:)=diff_p3_foot_velv(idx3(j),:);
        a_bunch3(j,:)=diff_p3_foot_accv(idx3(j),:);
        p_bunch4(j,:)=diff_p4_foot_posv(idx4(j),:);
        v_bunch4(j,:)=diff_p4_foot_velv(idx4(j),:);
        a_bunch4(j,:)=diff_p4_foot_accv(idx4(j),:);
       
    end
    
    bunch1m=nanmean(bunch1);
    p_bunch1m=nanmean(p_bunch1);
    v_bunch1m=nanmean(v_bunch1);
    a_bunch1m=nanmean(a_bunch1);
    
    % boot_impv1{i}=regress(bunch1m(100:300)',[p_bunch1m(100:300)' v_bunch1m(100:300)' a_bunch1m(100:300)']);
    C=[p_bunch1m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' v_bunch1m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' a_bunch1m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d=bunch1m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    boot_impv1(i,:)=lsqlin(C,d,A,B);
    
    bunch2m=nanmean(bunch2);
    p_bunch2m=nanmean(p_bunch2);
    v_bunch2m=nanmean(v_bunch2);
    a_bunch2m=nanmean(a_bunch2);
    
    % boot_impv1{i}=regress(bunch2m(100:300)',[p_bunch2m(100:300)' v_bunch2m(100:300)' a_bunch2m(100:300)']);
    C=[p_bunch2m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' v_bunch2m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' a_bunch2m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d=bunch2m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    boot_impv2(i,:)=lsqlin(C,d,A,B);
    
    bunch3m=nanmean(bunch3);
    p_bunch3m=nanmean(p_bunch3);
    v_bunch3m=nanmean(v_bunch3);
    a_bunch3m=nanmean(a_bunch3);
    
    % boot_impv1{i}=regress(bunch3m(100:300)',[p_bunch3m(100:300)' v_bunch3m(100:300)' a_bunch3m(100:300)']);
    C=[p_bunch3m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' v_bunch3m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' a_bunch3m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d=bunch3m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    boot_impv3(i,:)=lsqlin(C,d,A,B);
    
    
    bunch4m=nanmean(bunch4);
    p_bunch4m=nanmean(p_bunch4);
    v_bunch4m=nanmean(v_bunch4);
    a_bunch4m=nanmean(a_bunch4);
    
    % boot_impv1{i}=regress(bunch4m(100:300)',[p_bunch4m(100:300)' v_bunch4m(100:300)' a_bunch4m(100:300)']);
    C=[p_bunch4m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' v_bunch4m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)' a_bunch4m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)'];
    d=bunch4m(REGRESSION_START_INDEX:REGRESSION_END_INDEX)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*INERTIAL_LOWER_LIM;INERTIAL_UPPER_LIMIT];
    boot_impv4(i,:)=lsqlin(C,d,A,B);
   
end


boot_impv(1,:)=nanmean(boot_impv1);
boot_impv(2,:)=nanmean(boot_impv2);
boot_impv(3,:)=nanmean(boot_impv3);
boot_impv(4,:)=nanmean(boot_impv4);


ci=prctile(boot_impv1,[2.5 97.5])
boot_impvs(1,:)=ci(2,:)-boot_impv(1,1:3);
ci=prctile(boot_impv2,[2.5 97.5])
boot_impvs(2,:)=ci(2,:)-boot_impv(2,1:3);
ci=prctile(boot_impv3,[2.5 97.5])
boot_impvs(3,:)=ci(2,:)-boot_impv(3,1:3);
ci=prctile(boot_impv4,[2.5 97.5])
boot_impvs(4,:)=ci(2,:)-boot_impv(4,1:3);



save('bootstrap_impv.mat','boot_impv');

%%


vartor1=var(diff_p1_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX));
varimp1=var(diff_p1_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p1_foot_posvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(1,1)+diff_p1_foot_velvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(1,2)+diff_p1_foot_accvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(1,3)));
goodnessb(1,1)=100*(1-(varimp1/vartor1));




vartor2=var(diff_p2_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX));
varimp2=var(diff_p2_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p2_foot_posvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(2,1)+diff_p2_foot_velvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(2,2)+diff_p2_foot_accvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(2,3)));
goodnessb(2,1)=100*(1-(varimp2/vartor2))

vartor3=var(diff_p3_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX));
varimp3=var(diff_p3_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p3_foot_posvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(3,1)+diff_p3_foot_velvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(3,2)+diff_p3_foot_accvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(3,3)));
goodnessb(3,1)=100*(1-(varimp3/vartor3))

vartor4=var(diff_p4_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX));
varimp4=var(diff_p4_plat_torqueimpvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)-(diff_p4_foot_posvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(4,1)+diff_p4_foot_velvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(4,2)+diff_p4_foot_accvm(REGRESSION_START_INDEX:REGRESSION_END_INDEX)*boot_impv(4,3)));
goodnessb(4,1)=100*(1-(varimp4/vartor1));



boot_impv=[boot_impv,goodnessb]