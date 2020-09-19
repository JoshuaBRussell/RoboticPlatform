h = fopen('WEIGHT.DAT')
live_data=fread(h);
Input3= SimulinkRealTime.utils.getFileScopeData(live_data);

weight=nanmean(tsmovavg(Input3.data(1:40000,1),'s',1000,1));
tweight=weight/2;
weight=weight/9.8;