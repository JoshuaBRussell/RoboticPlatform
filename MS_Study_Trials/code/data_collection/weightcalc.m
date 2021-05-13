weight_file_handle = fopen('WEIGHT.DAT')
live_data=fread(weight_file_handle);
Input3= SimulinkRealTime.utils.getFileScopeData(live_data);

weight=nanmean(tsmovavg(Input3.data(1:40000,1),'s',1000,1));
tweight=weight/2;

fclose(weight_file_handle);