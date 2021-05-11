function [loadcell_bias] = pbase(loadcell_bias_file);
    %Finds the bias for each loadcell.
    %Input:  loadcell_bias_file: "*.DAT" filename as a string.
    %Output: struct where each field is the corresponding loadcell's DC bias.
    h = fopen(loadcell_bias_file)
    live_data=fread(h);
    Input3= SimulinkRealTime.utils.getFileScopeData(live_data);

    loadcell_bias.f1_off=nanmean(tsmovavg(Input3.data(1:40000,1),'s',1000,1));
    loadcell_bias.f2_off=nanmean(tsmovavg(Input3.data(1:40000,2),'s',1000,1));
    loadcell_bias.f3_off=nanmean(tsmovavg(Input3.data(1:40000,3),'s',1000,1));
    loadcell_bias.f4_off=nanmean(tsmovavg(Input3.data(1:40000,4),'s',1000,1));
    loadcell_bias.f5_off=nanmean(tsmovavg(Input3.data(1:40000,5),'s',1000,1));
    loadcell_bias.f6_off=nanmean(tsmovavg(Input3.data(1:40000,6),'s',1000,1));
    loadcell_bias.f7_off=nanmean(tsmovavg(Input3.data(1:40000,7),'s',1000,1));
    loadcell_bias.f8_off=nanmean(tsmovavg(Input3.data(1:40000,8),'s',1000,1));

    fclose all
    
end