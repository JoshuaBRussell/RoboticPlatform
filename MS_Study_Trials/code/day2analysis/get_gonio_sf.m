function [foot_gonio_sf, plat_gonio_sf] = get_gonio_sf(rel_data_path, perturb_type)

%----Needed Constants----%
DP_FOOT_GON_POS_SIG = 3;
DP_PLAT_GON_POS_SIG = 4;
DP_ENC_SIG = 1;

IE_FOOT_GON_POS_SIG = 3;
IE_PLAT_GON_POS_SIG = 4;
IE_ENC_SIG = 2;

GON_DELAY_IN_SAMPLES = 221;

if(perturb_type=='D' || perturb_type=='P')
    h = fopen(strcat(rel_data_path, 'DPFOOT','.dat'));
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);

    gonio=filtfilt(d1,Input1.data(:,DP_FOOT_GON_POS_SIG));
    gonio=detrend(gonio,'Linear');
    plot(gonio)
    [a,b]=findpeaks(gonio);
    encoder=filtfilt(d1,Input1.data(:,DP_ENC_SIG));

    [c,d]=findpeaks(encoder);
    for i=1:length(gonio)-GON_DELAY_IN_SAMPLES
    c_gonio(i)=gonio(i+GON_DELAY_IN_SAMPLES);
    end
    c_gonio=c_gonio-mean(c_gonio);
    DP_foot_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
    figure
    plot(DP_foot_gonio*c_gonio)
    hold on
    plot(encoder)

    % Platform goniometer
    gonio=filtfilt(d1,Input1.data(:,DP_PLAT_GON_POS_SIG));
    plot(gonio)
    [a,b]=findpeaks(gonio);
    encoder=filtfilt(d1,Input1.data(:,DP_ENC_SIG));

    [c,d]=findpeaks(encoder);
    for i=1:length(gonio)-GON_DELAY_IN_SAMPLES
    c_gonio(i)=gonio(i+GON_DELAY_IN_SAMPLES);
    end
    c_gonio=c_gonio-mean(c_gonio);
    DP_plat_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
    figure
    plot(DP_plat_gonio*c_gonio)
    hold on
    plot(encoder)

    foot_gonio_sf = DP_foot_gonio;
    plat_gonio_sf = DP_plat_gonio;

else
    h = fopen(strcat(rel_data_path, 'IEFOOT','.dat'));
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);

    gonio=filtfilt(d1,Input1.data(:,IE_FOOT_GON_POS_SIG));
    gonio=detrend(gonio,'Linear');
    plot(gonio)
    [a,b]=findpeaks(gonio);
    encoder=filtfilt(d1,Input1.data(:,IE_ENC_SIG));

    [c,d]=findpeaks(encoder);
    for i=1:length(gonio)-GON_DELAY_IN_SAMPLES
    c_gonio(i)=gonio(i+GON_DELAY_IN_SAMPLES);
    end
    c_gonio=c_gonio-mean(c_gonio);
    IE_foot_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
    figure
    plot(IE_foot_gonio*c_gonio)
    hold on
    plot(encoder)

    %% Platform goniometer
    gonio=filtfilt(d1,Input1.data(:,IE_PLAT_GON_POS_SIG));
    plot(gonio)
    [a,b]=findpeaks(gonio);
    encoder=filtfilt(d1,Input1.data(:,IE_ENC_SIG));

    [c,d]=findpeaks(encoder);
    for i=1:length(gonio)-GON_DELAY_IN_SAMPLES
    c_gonio(i)=gonio(i+GON_DELAY_IN_SAMPLES);
    end
    c_gonio=c_gonio-mean(c_gonio);
    IE_plat_gonio=regress(encoder(1:39700),c_gonio(1:39700)')
    figure
    plot(IE_plat_gonio*c_gonio)
    hold on
    plot(encoder)

    foot_gonio_sf = IE_foot_gonio;
    plat_gonio_sf = IE_plat_gonio;


    end
    pause(3);
    close all


end