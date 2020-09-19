clear all

%visual=0 for eyes closed, 1 for eyes opened
visual=1;
%Enter Subjects first initial
sub_initial='V';


if visual==0
    str='COPN';
elseif visual==1
    str='COP';
end

%%
for section = 1:2
    head={'Compliant',' Rigid'};
    for trial=1:3
    exp={'Trial 1','Trial 2','Trial 3'};
    h = fopen(strcat(sub_initial,str,num2str(section),num2str(trial),'.dat'));
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    emg_data{1} = Input1;
    d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
    save(strcat(sub_initial,str,num2str(section),num2str(trial),'.mat'),'Input1');

    coptorque=filtfilt(d1,Input1.data(:,18));
    weight=filtfilt(d1,Input1.data(:,17));
    ietorque=filtfilt(d1,Input1.data(:,8));
    for i=1:length(ietorque)
    dp_cop(i)=coptorque(i)*100/weight(i);
    ie_cop(i)=ietorque(i)*100/weight(i);

    end
    dp_pos=filtfilt(d1,Input1.data(:,5));
    ie_pos=filtfilt(d1,Input1.data(:,6));
    % dp_cop=dp_cop-mean(dp_cop);
    % ie_cop=ie_cop-mean(ie_cop);
    
    
  %plotting and recording 
    
  figure;
   
    subplot(1,2,1)
    scatter(ie_cop(1:120000),dp_cop(1:120000),2);
    
    xlabel('IE variance (cm)')
    ylabel('DP variance (cm)')
    title(strcat('CoP distribution ',head(section),' ', exp(trial)));

    axis([-2 2 3 13])
    subplot(1,2,2)

    scatter(ie_pos(1:120000),dp_pos(1:120000),2);
    title(strcat('Position distribution',head(section),' ',exp(trial)));
    axis([-10 10 -10 15])
     xlabel('IE variance (deg)')
    ylabel('DP variance (deg)')

    vardp(trial,section)=var(dp_cop(120000:1));
    varie(trial,section)=var(ie_cop(120000:1));
    
    
    end
  
end
