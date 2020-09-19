

vartor1=var(diff_p1_plat_torqueimpm(1120:1320));
varimp1=var(diff_p1_plat_torqueimpm(1120:1320)-(diff_p1_foot_posm(1120:1320)*p1impm(1)+diff_p1_foot_velm(1120:1320)*p1impm(2)+diff_p1_foot_accm(1120:1320)*p1impm(3)));
goodnessn(1,1)=100*(1-(varimp1/vartor1));


