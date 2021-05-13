function [mvc_data] = mv_reflex(mvc_filename);
    %This assumes only a single foot is on the platform.
    %Gets the nominal base CoP (in both IE and DP directions). 
    %This --should-- be independent of which platform is used (left or right).
    
    %In order to start making this more modular, this is now a function and 
    %will return the needed data in a struct.
    
    
    TA_EMG_SIG  = 1;
    SOL_EMG_SIG = 2;
    PL_EMG_SIG  = 3;
    GCA_EMG_SIG = 4;

    WEIGHT_SIG = 17; 

    DP_COP_TORQUE_SIG = 18;
    IE_COP_TORQUE_SIG = 8;
    DP_TORQUE_SIG = 7;

    h = fopen(mvc_filename)
    live_data2=fread(h);
    Input3= SimulinkRealTime.utils.getFileScopeData(live_data2);
    siz=size(Input3.data)
    siz=size(Input3.data)
    for i=1:siz(1,1)
    if(or(Input3.data(i,2)>10,Input3.data(i,2)<-10))
    Input3.data(i,2)=0;
    end
    end

    for i=1:siz(1,1)
    if(or(Input3.data(i,1)>10,Input3.data(i,1)<-10))
    Input3.data(i,1)=0;
    end
    end

    %----Find EMG DC Bias----%
    mvc_data.emg.off_TA=nanmean(tsmovavg(Input3.data(1:20000,TA_EMG_SIG),'s',1000,1));
    mvc_data.emg.off_SOL=nanmean(tsmovavg(Input3.data(1:20000,SOL_EMG_SIG),'s',1000,1));
    mvc_data.emg.off_PL=nanmean(tsmovavg(Input3.data(1:20000,PL_EMG_SIG),'s',1000,1));
    mvc_data.emg.off_GCA=nanmean(tsmovavg(Input3.data(1:20000,GCA_EMG_SIG),'s',1000,1));

    

    %----Find EMG MVC----%
    mvc_data.emg.TA_MVC=max((tsmovavg(abs(Input3.data(:,TA_EMG_SIG)-mvc_data.emg.off_TA),'s',1000,1)));
    mvc_data.emg.SOL_MVC=max((tsmovavg(abs(Input3.data(:,SOL_EMG_SIG)-mvc_data.emg.off_SOL),'s',1000,1)));
    mvc_data.emg.PL_MVC=max((tsmovavg(abs(Input3.data(:,PL_EMG_SIG)-mvc_data.emg.off_PL),'s',1000,1)));
    mvc_data.emg.GCA_MVC=max((tsmovavg(abs(Input3.data(:,GCA_EMG_SIG)-mvc_data.emg.off_GCA),'s',1000,1)));

    for i=2:siz(1,1)
    if(or(Input3.data(i,7)-Input3.data(i-1,7)>100,Input3.data(i,7)-Input3.data(i-1,7)<-100))
        Input3.data(i,7)=Input3.data(i-1,7);
    end
    end
    for i=2:siz(1,1)
    if(abs(Input3.data(i,8)-Input3.data(i-1,8))>100)
        Input3.data(i,8)=Input3.data(i-1,8);
    end
    end
    mvc_data.DPbase=nanmean(tsmovavg(Input3.data(1:20000,DP_TORQUE_SIG),'s',1000,1));
    mvc_data.IEbase=nanmean(tsmovavg(Input3.data(1:20000,IE_COP_TORQUE_SIG),'s',1000,1));
    mvc_data.CoPbase=nanmean(tsmovavg(Input3.data(1:20000,DP_COP_TORQUE_SIG),'s',1000,1));
    mvc_data.single_weight = mean(Input3.data(1:20000, WEIGHT_SIG)); %Gets the base weight while a single foot is on the platform. 

    


    %This uses the torque/(actual weight) on the platform for a single foot (instead of
    %just assuming the weight is evenly distributed between the both legs (i.e. torque/(weight/2))).
    mvc_data.DP_cop=mean(Input3.data(1:20000, DP_COP_TORQUE_SIG)./Input3.data(1:20000, WEIGHT_SIG)) * 100; %in cm - x100 converts to cm
    mvc_data.IE_cop=mean(Input3.data(1:20000, IE_COP_TORQUE_SIG)./Input3.data(1:20000, WEIGHT_SIG)) * 100; %in cm - x100 converts to cm






    %----EMG Plot----%
    figure
    subplot(4,1,1)
    plot((tsmovavg(abs(Input3.data(:,1)-mvc_data.emg.off_TA),'s',1000,1)));
    title("TA");
    subplot(4,1,2)
    plot((tsmovavg(abs(Input3.data(:,2)-mvc_data.emg.off_SOL),'s',1000,1)));
    title("SOL");
    subplot(4,1,3)
    plot((tsmovavg(abs(Input3.data(:,3)-mvc_data.emg.off_PL),'s',1000,1)));
    title("PL");
    subplot(4,1,4)
    plot((tsmovavg(abs(Input3.data(:,4)-mvc_data.emg.off_GCA),'s',1000,1)));
    title("GCA");
    fclose(h);


end