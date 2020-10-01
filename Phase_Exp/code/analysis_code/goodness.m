

vartor1=var(diff_p1_plat_torqueimpm(700:900));
varimp1=var(diff_p1_plat_torqueimpm(700:900)-(diff_p1_foot_posm(700:900)*p1impm(1)+diff_p1_foot_velm(700:900)*p1impm(2)+diff_p1_foot_accm(700:900)*p1impm(3)));
goodnessn(1,1)=100*(1-(varimp1/vartor1));




vartor2=var(diff_p2_plat_torqueimpm(1000:1200));
varimp2=var(diff_p2_plat_torqueimpm(1000:1200)-(diff_p2_foot_posm(1000:1200)*p2impm(1)+diff_p2_foot_velm(1000:1200)*p2impm(2)+diff_p2_foot_accm(1000:1200)*p2impm(3)));
goodnessn(2,1)=100*(1-(varimp2/vartor2))

vartor3=var(diff_p3_plat_torqueimpm(1280:1480));
varimp3=var(diff_p3_plat_torqueimpm(1280:1480)-(diff_p3_foot_posm(1280:1480)*p3impm(1)+diff_p3_foot_velm(1280:1480)*p3impm(2)+diff_p3_foot_accm(1280:1480)*p3impm(3)));
goodnessn(3,1)=100*(1-(varimp3/vartor3))

vartor4=var(diff_p4_plat_torqueimpm(720:920));
varimp4=var(diff_p4_plat_torqueimpm(720:920)-(diff_p4_foot_posm(720:920)*p4impm(1)+diff_p4_foot_velm(720:920)*p4impm(2)+diff_p4_foot_accm(720:920)*p4impm(3)));
goodnessn(4,1)=100*(1-(varimp4/vartor1));




vartor5=var(diff_p5_plat_torqueimpm(1020:1220));
varimp5=var(diff_p5_plat_torqueimpm(1020:1220)-(diff_p5_foot_posm(1020:1220)*p5impm(1)+diff_p5_foot_velm(1020:1220)*p5impm(2)+diff_p5_foot_accm(1020:1220)*p5impm(3)));
goodnessn(5,1)=100*(1-(varimp5/vartor2))

vartor6=var(diff_p6_plat_torqueimpm(1300:1500));
varimp6=var(diff_p6_plat_torqueimpm(1300:1500)-(diff_p6_foot_posm(1300:1500)*p6impm(1)+diff_p6_foot_velm(1300:1500)*p6impm(2)+diff_p6_foot_accm(1300:1500)*p6impm(3)));
goodnessn(6,1)=100*(1-(varimp6/vartor3))