function [TA_NORM, SOL_NORM, PL_NORM, GCA_NORM] = find_emg_normalization_walking(DATA_DIR, DATA_FILE_NAME)

%% Used for when certain subjects do not step on the platform during the EMG normalization section 
%% of the data collection.

%% In this case we just take the maximum EMG during a brief walking section.


TA_EMG_SIG  = 1;
SOL_EMG_SIG = 2;
PL_EMG_SIG  = 3;
GCA_EMG_SIG = 4;


d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);



h = fopen(strcat(DATA_DIR, DATA_FILE_NAME));

live_data=fread(h);
Input1 = SimulinkRealTime.utils.getFileScopeData(live_data);


ta=Input1.data(:,TA_EMG_SIG);
sol=Input1.data(:,SOL_EMG_SIG);
pl=Input1.data(:,PL_EMG_SIG);
gca=Input1.data(:,GCA_EMG_SIG); 

off_TA = mean(ta);
off_SOL = mean(sol);
off_PL = mean(pl);
off_GCA = mean(gca);

ta = ta - off_TA;
sol = sol - off_SOL;
pl = pl - off_PL;
gca = gca - off_GCA;

ta = filtfilt(d1,ta);
sol = filtfilt(d1,sol);
pl = filtfilt(d1,pl);
gca = filtfilt(d1,gca);

TA_NORM  = max(ta);
SOL_NORM = max(sol);
PL_NORM  = max(pl);
GCA_NORM = max(gca);


end

