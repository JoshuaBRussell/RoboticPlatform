

vartor1=var(diff_p1_plat_torqueimpm(700:900));
varimp1=var(diff_p1_plat_torqueimpm(700:900)-(diff_p1_foot_posm(700:900)*p1impm2(1)+diff_p1_foot_velm(700:900)*p1impm2(2)+diff_p1_foot_accm(700:900)*p1impm2(3)));
goodnessc(1,1)=100*(1-(varimp1/vartor1))




vartor2=var(diff_p2_plat_torqueimpm(1000:1200));
varimp2=var(diff_p2_plat_torqueimpm(1000:1200)-(diff_p2_foot_posm(1000:1200)*p2impm2(1)+diff_p2_foot_velm(1000:1200)*p2impm2(2)+diff_p2_foot_accm(1000:1200)*p2impm2(3)));
goodnessc(2,1)=100*(1-(varimp2/vartor2));

vartor3=var(diff_p3_plat_torqueimpm(1280:1480));
varimp3=var(diff_p3_plat_torqueimpm(1280:1480)-(diff_p3_foot_posm(1280:1480)*p3impm2(1)+diff_p3_foot_velm(1280:1480)*p3impm2(2)+diff_p3_foot_accm(1280:1480)*p3impm2(3)));
goodnessc(3,1)=100*(1-(varimp3/vartor3))

vartor4=var(diff_p4_plat_torqueimpm(720:920));
varimp4=var(diff_p4_plat_torqueimpm(720:920)-(diff_p4_foot_posm(720:920)*p4impm2(1)+diff_p4_foot_velm(720:920)*p4impm2(2)+diff_p4_foot_accm(720:920)*p4impm2(3)));
goodnessc(4,1)=100*(1-(varimp4/vartor1))




vartor5=var(diff_p5_plat_torqueimpm(1000:1200));
varimp5=var(diff_p5_plat_torqueimpm(1000:1200)-(diff_p5_foot_posm(1000:1200)*p5impm2(1)+diff_p5_foot_velm(1000:1200)*p5impm2(2)+diff_p5_foot_accm(1000:1200)*p5impm2(3)));
goodnessc(5,1)=100*(1-(varimp5/vartor2));

vartor6=var(diff_p6_plat_torqueimpm(1280:1480));
varimp6=var(diff_p6_plat_torqueimpm(1280:1480)-(diff_p6_foot_posm(1280:1480)*p6impm2(1)+diff_p6_foot_velm(1280:1480)*p6impm2(2)+diff_p6_foot_accm(1280:1480)*p6impm2(3)));
goodnessc(6,1)=100*(1-(varimp6/vartor3))