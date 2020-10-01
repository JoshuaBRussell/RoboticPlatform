clear boot_impvaf;
diff_p1_plat_torqueimpfin=diff_p1_plat_torqueimps(1:p1_fin,:);
diff_p1_foot_posfin=diff_p1_foot_poss(1:p1_fin,:);
diff_p1_foot_velfin=diff_p1_foot_vels(1:p1_fin,:);
diff_p1_foot_accfin=diff_p1_foot_accs(1:p1_fin,:);

diff_p2_plat_torqueimpfin=diff_p2_plat_torqueimps(1:p2_fin,:);
diff_p2_foot_posfin=diff_p2_foot_poss(1:p2_fin,:);
diff_p2_foot_velfin=diff_p2_foot_vels(1:p2_fin,:);
diff_p2_foot_accfin=diff_p2_foot_accs(1:p2_fin,:);

diff_p3_plat_torqueimpfin=diff_p3_plat_torqueimps(1:p3_fin,:);
diff_p3_foot_posfin=diff_p3_foot_poss(1:p3_fin,:);
diff_p3_foot_velfin=diff_p3_foot_vels(1:p3_fin,:);
diff_p3_foot_accfin=diff_p3_foot_accs(1:p3_fin,:);


diff_p4_plat_torqueimpfin=diff_p4_plat_torqueimps(1:p4_fin,:);
diff_p4_foot_posfin=diff_p4_foot_poss(1:p4_fin,:);
diff_p4_foot_velfin=diff_p4_foot_vels(1:p4_fin,:);
diff_p4_foot_accfin=diff_p4_foot_accs(1:p4_fin,:);

diff_p5_plat_torqueimpfin=diff_p5_plat_torqueimps(1:p5_fin,:);
diff_p5_foot_posfin=diff_p5_foot_poss(1:p5_fin,:);
diff_p5_foot_velfin=diff_p5_foot_vels(1:p5_fin,:);
diff_p5_foot_accfin=diff_p5_foot_accs(1:p5_fin,:);

diff_p6_plat_torqueimpfin=diff_p6_plat_torqueimps(1:p6_fin,:);
diff_p6_foot_posfin=diff_p6_foot_poss(1:p6_fin,:);
diff_p6_foot_velfin=diff_p6_foot_vels(1:p6_fin,:);
diff_p6_foot_accfin=diff_p6_foot_accs(1:p6_fin,:);






%%




for i=1:loops
    unique1=0;
    unique2=0;
    unique3=0;
    unique4=0;
    unique5=0;
    unique6=0;
    
    while(unique1<60)
        [vafbunch1,idx1]=datasample(diff_p1_plat_torqueimpv,p1_fin);
        numb=length(unique(idx1));
        unique1=100*(numb/p1_fin);
    end
    unique1=0;
    while(unique2<60)
        [vafbunch2,idx2]=datasample(diff_p2_plat_torqueimpv,p2_fin);
        numb=length(unique(idx2));
        unique2=100*(numb/p2_fin);
    end
    unique2=0;
    while(unique3<60)
        [vafbunch3,idx3]=datasample(diff_p3_plat_torqueimpv,p3_fin);
        numb=length(unique(idx3));
        unique3=100*(numb/p3_fin);
    end
    unique3=0;
    while(unique4<60)
        [vafbunch4,idx4]=datasample(diff_p4_plat_torqueimpv,p4_fin);
        numb=length(unique(idx4));
        unique4=100*(numb/p4_fin);
    end
    unique4=0;
    while(unique5<60)
        [vafbunch5,idx5]=datasample(diff_p5_plat_torqueimpv,p5_fin);
        numb=length(unique(idx5));
        unique5=100*(numb/p5_fin);
    end
    unique5=0;
    while(unique6<60)
        [vafbunch6,idx6]=datasample(diff_p6_plat_torqueimpv,p6_fin);
        numb=length(unique(idx6));
        unique6=100*(numb/p6_fin);
    end
    unique6=0;
    
    for j=1:p1_fin
        p_vafbunch1(j,:)=diff_p1_foot_posv(idx1(j),:);
        v_vafbunch1(j,:)=diff_p1_foot_velv(idx1(j),:);
        a_vafbunch1(j,:)=diff_p1_foot_accv(idx1(j),:);
    end
    for j=1:p2_fin
        p_vafbunch2(j,:)=diff_p2_foot_posv(idx2(j),:);
        v_vafbunch2(j,:)=diff_p2_foot_velv(idx2(j),:);
        a_vafbunch2(j,:)=diff_p2_foot_accv(idx2(j),:);
    end
    for j=1:p3_fin
        p_vafbunch3(j,:)=diff_p3_foot_posv(idx3(j),:);
        v_vafbunch3(j,:)=diff_p3_foot_velv(idx3(j),:);
        a_vafbunch3(j,:)=diff_p3_foot_accv(idx3(j),:);
    end
    for j=1:p4_fin
        p_vafbunch4(j,:)=diff_p4_foot_posv(idx4(j),:);
        v_vafbunch4(j,:)=diff_p4_foot_velv(idx4(j),:);
        a_vafbunch4(j,:)=diff_p4_foot_accv(idx4(j),:);
    end
    for j=1:p5_fin
        p_vafbunch5(j,:)=diff_p5_foot_posv(idx5(j),:);
        v_vafbunch5(j,:)=diff_p5_foot_velv(idx5(j),:);
        a_vafbunch5(j,:)=diff_p5_foot_accv(idx5(j),:);
    end
    for j=1:p6_fin
        p_vafbunch6(j,:)=diff_p6_foot_posv(idx6(j),:);
        v_vafbunch6(j,:)=diff_p6_foot_velv(idx6(j),:);
        a_vafbunch6(j,:)=diff_p6_foot_accv(idx6(j),:);
    end
    
    vafbunch1m=nanmean(vafbunch1);
    p_vafbunch1m=nanmean(p_vafbunch1);
    v_vafbunch1m=nanmean(v_vafbunch1);
    a_vafbunch1m=nanmean(a_vafbunch1);
    
    % boot_impvaf1{i}=regress(vafbunch1m(700:900)',[p_vafbunch1m(700:900)' v_vafbunch1m(700:900)' a_vafbunch1m(700:900)']);
    C=[p_vafbunch1m(700:900)' v_vafbunch1m(700:900)' a_vafbunch1m(700:900)'];
    d=vafbunch1m(700:900)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;l_lim];
    boot_impvaf1(i,:)=lsqlin(C,d,A,B);
    
    vafbunch2m=nanmean(vafbunch2);
    p_vafbunch2m=nanmean(p_vafbunch2);
    v_vafbunch2m=nanmean(v_vafbunch2);
    a_vafbunch2m=nanmean(a_vafbunch2);
    
    % boot_impvaf1{i}=regress(vafbunch2m(1000:1200)',[p_vafbunch2m(1000:1200)' v_vafbunch2m(1000:1200)' a_vafbunch2m(1000:1200)']);
    C=[p_vafbunch2m(1000:1200)' v_vafbunch2m(1000:1200)' a_vafbunch2m(1000:1200)'];
    d=vafbunch2m(1000:1200)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;l_lim];
    boot_impvaf2(i,:)=lsqlin(C,d,A,B);
    
    vafbunch3m=nanmean(vafbunch3);
    p_vafbunch3m=nanmean(p_vafbunch3);
    v_vafbunch3m=nanmean(v_vafbunch3);
    a_vafbunch3m=nanmean(a_vafbunch3);
    
    % boot_impvaf1{i}=regress(vafbunch3m(1280:1480)',[p_vafbunch3m(1280:1480)' v_vafbunch3m(1280:1480)' a_vafbunch3m(1280:1480)']);
    C=[p_vafbunch3m(1280:1480)' v_vafbunch3m(1280:1480)' a_vafbunch3m(1280:1480)'];
    d=vafbunch3m(1280:1480)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;l_lim];
    boot_impvaf3(i,:)=lsqlin(C,d,A,B);
    
    
    vafbunch4m=nanmean(vafbunch4);
    p_vafbunch4m=nanmean(p_vafbunch4);
    v_vafbunch4m=nanmean(v_vafbunch4);
    a_vafbunch4m=nanmean(a_vafbunch4);
    
    % boot_impvaf1{i}=regress(vafbunch4m(720:920)',[p_vafbunch4m(720:920)' v_vafbunch4m(720:920)' a_vafbunch4m(720:920)']);
    C=[p_vafbunch4m(720:920)' v_vafbunch4m(720:920)' a_vafbunch4m(720:920)'];
    d=vafbunch4m(720:920)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;l_lim];
    boot_impvaf4(i,:)=lsqlin(C,d,A,B);
    
    vafbunch5m=nanmean(vafbunch5);
    p_vafbunch5m=nanmean(p_vafbunch5);
    v_vafbunch5m=nanmean(v_vafbunch5);
    a_vafbunch5m=nanmean(a_vafbunch5);
    
    % boot_impvaf1{i}=regress(vafbunch5m(1020:1220)',[p_vafbunch5m(1020:1220)' v_vafbunch5m(1020:1220)' a_vafbunch5m(1020:1220)']);
    C=[p_vafbunch5m(1020:1220)' v_vafbunch5m(1020:1220)' a_vafbunch5m(1020:1220)'];
    d=vafbunch5m(1020:1220)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;l_lim];
    boot_impvaf5(i,:)=lsqlin(C,d,A,B);
    
    vafbunch6m=nanmean(vafbunch6);
    p_vafbunch6m=nanmean(p_vafbunch6);
    v_vafbunch6m=nanmean(v_vafbunch6);
    a_vafbunch6m=nanmean(a_vafbunch6);
    
    % boot_impvaf1{i}=regress(vafbunch6m(1300:1500)',[p_vafbunch6m(1300:1500)' v_vafbunch6m(1300:1500)' a_vafbunch6m(1300:1500)']);
    C=[p_vafbunch6m(1300:1500)' v_vafbunch6m(1300:1500)' a_vafbunch6m(1300:1500)'];
    d=vafbunch6m(1300:1500)';
    A=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
    B=[0 ;0 ;1000;1000;-1*lim;l_lim];
    boot_impvaf6(i,:)=lsqlin(C,d,A,B);
    
    
end


boot_impvaf(1,:)=nanmean(boot_impvaf1);
boot_impvaf(2,:)=nanmean(boot_impvaf2);
boot_impvaf(3,:)=nanmean(boot_impvaf3);
boot_impvaf(4,:)=nanmean(boot_impvaf4);
boot_impvaf(5,:)=nanmean(boot_impvaf5);
boot_impvaf(6,:)=nanmean(boot_impvaf6);
diff_p1_plat_torqueimpfinm=mean(diff_p1_plat_torqueimpfin);
diff_p1_foot_posfinm=mean(diff_p1_foot_posfin);
diff_p1_foot_velfinm=mean(diff_p1_foot_velfin);
diff_p1_foot_accfinm=mean(diff_p1_foot_accfin);

diff_p2_plat_torqueimpfinm=mean(diff_p2_plat_torqueimpfin);
diff_p2_foot_posfinm=mean(diff_p2_foot_posfin);
diff_p2_foot_velfinm=mean(diff_p2_foot_velfin);
diff_p2_foot_accfinm=mean(diff_p2_foot_accfin);

diff_p3_plat_torqueimpfinm=mean(diff_p3_plat_torqueimpfin);
diff_p3_foot_posfinm=mean(diff_p3_foot_posfin);
diff_p3_foot_velfinm=mean(diff_p3_foot_velfin);
diff_p3_foot_accfinm=mean(diff_p3_foot_accfin);

diff_p4_plat_torqueimpfinm=mean(diff_p4_plat_torqueimpfin);
diff_p4_foot_posfinm=mean(diff_p4_foot_posfin);
diff_p4_foot_velfinm=mean(diff_p4_foot_velfin);
diff_p4_foot_accfinm=mean(diff_p4_foot_accfin);

diff_p5_plat_torqueimpfinm=mean(diff_p5_plat_torqueimpfin);
diff_p5_foot_posfinm=mean(diff_p5_foot_posfin);
diff_p5_foot_velfinm=mean(diff_p5_foot_velfin);
diff_p5_foot_accfinm=mean(diff_p5_foot_accfin);

diff_p6_plat_torqueimpfinm=mean(diff_p6_plat_torqueimpfin);
diff_p6_foot_posfinm=mean(diff_p6_foot_posfin);
diff_p6_foot_velfinm=mean(diff_p6_foot_velfin);
diff_p6_foot_accfinm=mean(diff_p6_foot_accfin);


%%
vartor1=var(diff_p1_plat_torqueimpfinm(700:900));
varimp1=var(diff_p1_plat_torqueimpfinm(700:900)-(diff_p1_foot_posfinm(700:900)*boot_impvaf(1,1)+diff_p1_foot_velfinm(700:900)*boot_impvaf(1,2)+diff_p1_foot_accfinm(700:900)*boot_impvaf(1,3)));
goodnessb(1,1)=100*(1-(varimp1/vartor1));




vartor2=var(diff_p2_plat_torqueimpfinm(1000:1200));
varimp2=var(diff_p2_plat_torqueimpfinm(1000:1200)-(diff_p2_foot_posfinm(1000:1200)*boot_impvaf(2,1)+diff_p2_foot_velfinm(1000:1200)*boot_impvaf(2,2)+diff_p2_foot_accfinm(1000:1200)*boot_impvaf(2,3)));
goodnessb(2,1)=100*(1-(varimp2/vartor2))

vartor3=var(diff_p3_plat_torqueimpfinm(1280:1480));
varimp3=var(diff_p3_plat_torqueimpfinm(1280:1480)-(diff_p3_foot_posfinm(1280:1480)*boot_impvaf(3,1)+diff_p3_foot_velfinm(1280:1480)*boot_impvaf(3,2)+diff_p3_foot_accfinm(1280:1480)*boot_impvaf(3,3)));
goodnessb(3,1)=100*(1-(varimp3/vartor3))

vartor4=var(diff_p4_plat_torqueimpfinm(720:920));
varimp4=var(diff_p4_plat_torqueimpfinm(720:920)-(diff_p4_foot_posfinm(720:920)*boot_impvaf(4,1)+diff_p4_foot_velfinm(720:920)*boot_impvaf(4,2)+diff_p4_foot_accfinm(720:920)*boot_impvaf(4,3)));
goodnessb(4,1)=100*(1-(varimp4/vartor1));




vartor5=var(diff_p5_plat_torqueimpfinm(1020:1220));
varimp5=var(diff_p5_plat_torqueimpfinm(1020:1220)-(diff_p5_foot_posfinm(1020:1220)*boot_impvaf(5,1)+diff_p5_foot_velfinm(1020:1220)*boot_impvaf(5,2)+diff_p5_foot_accfinm(1020:1220)*boot_impvaf(5,3)));
goodnessb(5,1)=100*(1-(varimp5/vartor2))

vartor6=var(diff_p6_plat_torqueimpfinm(1300:1500));
varimp6=var(diff_p6_plat_torqueimpfinm(1300:1500)-(diff_p6_foot_posfinm(1300:1500)*boot_impvaf(6,1)+diff_p6_foot_velfinm(1300:1500)*boot_impvaf(6,2)+diff_p6_foot_accfinm(1300:1500)*boot_impvaf(6,3)));
goodnessb(6,1)=100*(1-(varimp6/vartor3))

boot_impvaf=[boot_impvaf,goodnessb]