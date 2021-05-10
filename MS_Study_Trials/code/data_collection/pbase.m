h = fopen('PBASE.DAT')
live_data=fread(h);
Input3= SimulinkRealTime.utils.getFileScopeData(live_data);

f1_off=nanmean(tsmovavg(Input3.data(1:40000,1),'s',1000,1));
f2_off=nanmean(tsmovavg(Input3.data(1:40000,2),'s',1000,1));
f3_off=nanmean(tsmovavg(Input3.data(1:40000,3),'s',1000,1));
f4_off=nanmean(tsmovavg(Input3.data(1:40000,4),'s',1000,1));
f5_off=nanmean(tsmovavg(Input3.data(1:40000,5),'s',1000,1));
f6_off=nanmean(tsmovavg(Input3.data(1:40000,6),'s',1000,1));
f7_off=nanmean(tsmovavg(Input3.data(1:40000,7),'s',1000,1));
f8_off=nanmean(tsmovavg(Input3.data(1:40000,8),'s',1000,1));

fclose all