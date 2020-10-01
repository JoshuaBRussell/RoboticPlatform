
ideal=floor((analysis_value-1)*1);


for i=1:loops
    unique1=0;
    unique2=0;
    unique3=0;
    unique4=0;
    unique5=0;
    unique6=0;
    
    while(unique1<60)
        [bunch1,idx1]=datasample(diff_p1_plat_torqueimp,ideal);
        numb=length(unique(idx1));
        unique1=100*(numb/ideal);
    end
    unique1=0;
    while(unique2<60)
        [bunch2,idx2]=datasample(diff_p2_plat_torqueimp,ideal);
        numb=length(unique(idx2));
        unique2=100*(numb/ideal);
    end
    unique2=0;
    while(unique3<60)
        [bunch3,idx3]=datasample(diff_p3_plat_torqueimp,ideal);
        numb=length(unique(idx3));
        unique3=100*(numb/ideal);
    end
    unique3=0;
    while(unique4<60)
        [bunch4,idx4]=datasample(diff_p4_plat_torqueimp,ideal);
        numb=length(unique(idx4));
        unique4=100*(numb/ideal);
    end
    unique4=0;
    while(unique5<60)
        [bunch5,idx5]=datasample(diff_p5_plat_torqueimp,ideal);
        numb=length(unique(idx5));
        unique5=100*(numb/ideal);
    end
    unique5=0;
    while(unique6<60)
        [bunch6,idx6]=datasample(diff_p6_plat_torqueimp,ideal);
        numb=length(unique(idx6));
        unique6=100*(numb/ideal);
    end
    unique6=0;
    
    for j=1:ideal
        p_bunch1(j,:)=diff_p1_foot_pos(idx1(j),:);
        v_bunch1(j,:)=diff_p1_foot_vel(idx1(j),:);
        a_bunch1(j,:)=diff_p1_foot_acc(idx1(j),:);
        p_bunch2(j,:)=diff_p2_foot_pos(idx2(j),:);
        v_bunch2(j,:)=diff_p2_foot_vel(idx2(j),:);
        a_bunch2(j,:)=diff_p2_foot_acc(idx2(j),:);
        p_bunch3(j,:)=diff_p3_foot_pos(idx3(j),:);
        v_bunch3(j,:)=diff_p3_foot_vel(idx3(j),:);
        a_bunch3(j,:)=diff_p3_foot_acc(idx3(j),:);
        p_bunch4(j,:)=diff_p4_foot_pos(idx4(j),:);
        v_bunch4(j,:)=diff_p4_foot_vel(idx4(j),:);
        a_bunch4(j,:)=diff_p4_foot_acc(idx4(j),:);
        p_bunch5(j,:)=diff_p5_foot_pos(idx5(j),:);
        v_bunch5(j,:)=diff_p5_foot_vel(idx5(j),:);
        a_bunch5(j,:)=diff_p5_foot_acc(idx5(j),:);
        p_bunch6(j,:)=diff_p6_foot_pos(idx6(j),:);
        v_bunch6(j,:)=diff_p6_foot_vel(idx6(j),:);
        a_bunch6(j,:)=diff_p6_foot_acc(idx6(j),:);
    end
    
    bunch1m=nanmean(bunch1);
    p_bunch1m=nanmean(p_bunch1);
    v_bunch1m=nanmean(v_bunch1);
    a_bunch1m=nanmean(a_bunch1);
    
    % boot_imp1{i}=regress(bunch1m(700:900)',[p_bunch1m(700:900)' v_bunch1m(700:900)' a_bunch1m(700:900)']);
    C=[p_bunch1m(700:900)' v_bunch1m(700:900)' a_bunch1m(700:900)'];
    d=bunch1m(700:900)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    boot_imp1(i,:)=lsqlin(C,d,A,B);
    
    bunch2m=nanmean(bunch2);
    p_bunch2m=nanmean(p_bunch2);
    v_bunch2m=nanmean(v_bunch2);
    a_bunch2m=nanmean(a_bunch2);
    
    % boot_imp1{i}=regress(bunch2m(1000:1200)',[p_bunch2m(1000:1200)' v_bunch2m(1000:1200)' a_bunch2m(1000:1200)']);
    C=[p_bunch2m(1000:1200)' v_bunch2m(1000:1200)' a_bunch2m(1000:1200)'];
    d=bunch2m(1000:1200)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    boot_imp2(i,:)=lsqlin(C,d,A,B);
    
    bunch3m=nanmean(bunch3);
    p_bunch3m=nanmean(p_bunch3);
    v_bunch3m=nanmean(v_bunch3);
    a_bunch3m=nanmean(a_bunch3);
    
    % boot_imp1{i}=regress(bunch3m(1280:1480)',[p_bunch3m(1280:1480)' v_bunch3m(1280:1480)' a_bunch3m(1280:1480)']);
    C=[p_bunch3m(1280:1480)' v_bunch3m(1280:1480)' a_bunch3m(1280:1480)'];
    d=bunch3m(1280:1480)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    boot_imp3(i,:)=lsqlin(C,d,A,B);
    
    
    bunch4m=nanmean(bunch4);
    p_bunch4m=nanmean(p_bunch4);
    v_bunch4m=nanmean(v_bunch4);
    a_bunch4m=nanmean(a_bunch4);
    
    % boot_imp1{i}=regress(bunch4m(720:920)',[p_bunch4m(720:920)' v_bunch4m(720:920)' a_bunch4m(720:920)']);
    C=[p_bunch4m(720:920)' v_bunch4m(720:920)' a_bunch4m(720:920)'];
    d=bunch4m(720:920)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    boot_imp4(i,:)=lsqlin(C,d,A,B);
    
    bunch5m=nanmean(bunch5);
    p_bunch5m=nanmean(p_bunch5);
    v_bunch5m=nanmean(v_bunch5);
    a_bunch5m=nanmean(a_bunch5);
    
    % boot_imp1{i}=regress(bunch5m(1020:1220)',[p_bunch5m(1020:1220)' v_bunch5m(1020:1220)' a_bunch5m(1020:1220)']);
    C=[p_bunch5m(1020:1220)' v_bunch5m(1020:1220)' a_bunch5m(1020:1220)'];
    d=bunch5m(1020:1220)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    boot_imp5(i,:)=lsqlin(C,d,A,B);
    
    bunch6m=nanmean(bunch6);
    p_bunch6m=nanmean(p_bunch6);
    v_bunch6m=nanmean(v_bunch6);
    a_bunch6m=nanmean(a_bunch6);
    
    % boot_imp1{i}=regress(bunch6m(1300:1500)',[p_bunch6m(1300:1500)' v_bunch6m(1300:1500)' a_bunch6m(1300:1500)']);
    C=[p_bunch6m(1300:1500)' v_bunch6m(1300:1500)' a_bunch6m(1300:1500)'];
    d=bunch6m(1300:1500)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;0.07];
    boot_imp6(i,:)=lsqlin(C,d,A,B);
    
    
end


boot_imp(1,:)=nanmean(boot_imp1);
boot_imp(2,:)=nanmean(boot_imp2);
boot_imp(3,:)=nanmean(boot_imp3);
boot_imp(4,:)=nanmean(boot_imp4);
boot_imp(5,:)=nanmean(boot_imp5);
boot_imp(6,:)=nanmean(boot_imp6);

save('bootstrap_imp.mat','boot_imp');
