function [emg_means] = process_EMG(EMG_p1, EMG_p2, EMG_p3 ,EMG_p4)


emg_means = [mean(EMG_p1.TA), mean(EMG_p1.PL), mean(EMG_p1.SOL), mean(EMG_p1.GCA);
             mean(EMG_p2.TA), mean(EMG_p2.PL), mean(EMG_p2.SOL), mean(EMG_p2.GCA);
             mean(EMG_p3.TA), mean(EMG_p3.PL), mean(EMG_p3.SOL), mean(EMG_p3.GCA);
             mean(EMG_p4.TA), mean(EMG_p4.PL), mean(EMG_p4.SOL), mean(EMG_p4.GCA);];






end

