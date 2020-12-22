close all
clear all
mvc_evaluation;
% Insert subject initial and name. Make sure it matches the format for naming
sub_initial='C';
NUM_TRAINING_TRIALS = 40;
[num,den] = butter(2,10/1000);
[filt_num, filt_den] = butter(2,10/1000);
% Add the trials you want to exclude in here
exclude=[0];%exclude=[1,2];
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','Samplerate',2000);
%% Section to calculate goniometer gains
t=gonio_values_func;
DP_foot_gonio=t(1);
DP_plat_gonio=t(2);
%close all    %uncomment to close all figures
%%
i=0;
p1=1;
p2=1;
p3=1;
p0=1;
p4=1;
p5=1;
p6=1;
p7=1;

for trials=1:2
    
  if(ismember(trials,exclude)==0)  
    if(trials<10)
        h = fopen(strcat(sub_initial,'PW1',num2str(trials),'.dat'));
    end
    
    live_data=fread(h);
    Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
    siz=size(Input1.data);
    
    %             %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','CWamplerate',2000);
    %d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 20, 'CWampleRate', 2000, 'DesignMethod', 'window');
    pert_torque=filtfilt(d1,Input1.data(:,7));
    f1=Input1.data(:,9)*53.4;
    
    f2=Input1.data(:,10)*53.4;

    f3=Input1.data(:,11)*53.4;

    f4=Input1.data(:,12)*53.4;
 
    f5=Input1.data(:,20)*53.4/2;
  
    f6=Input1.data(:,21)*53.4/2;
 
  
    w1=(Input1.data(:,18));
  
    flag=Input1.data(:,17);
    foot_pos_data=Input1.data(:,13);
    foot_const=mean(foot_pos_data);
    foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
    plat_pos_data=Input1.data(:,14);
    plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
    ramp=Input1.data(:,22);
    
    filtered_weight = Input1.data(:,23);
    filtered_fy = Input1.data(:,24);
    phase_variable = Input1.data(:,25);
    
    [test,peaks]=findpeaks(Input1.data(:,17));
    
    for i=1:length(peaks)
        
        time=[-200:0.5:1000];
       
        if test(i)==4
            force1r_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f2(peaks(i)-360);
            force1r_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f3(peaks(i)-360);
            force1r_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f4(peaks(i)-360);
            force1r_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f5(peaks(i)-360);
            force1r_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f6(peaks(i)-360);
            filtered_weightr(p0,:)= filtered_weight(peaks(i)-400:peaks(i)+2000)-filtered_weight(peaks(i)-360);
            filtered_fyr(p0,:)    = filtered_fy(peaks(i)-400:peaks(i)+2000)-filtered_fy(peaks(i)-360);
            phase_variabler(p0,:) = phase_variable(peaks(i)-400:peaks(i)+2000)-phase_variable(peaks(i)-360);
            
            weight1r(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
            
            p1r_plat_pos(p0,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
            p1r_foot_pos(p0,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
            phase1r(p0,:)=ramp(peaks(i)-400:peaks(i)+2000);
            p0=p0+1;
        end
        
       
       
        if test(i)==8
            force1h_1(p7,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1h_2(p7,:)=f2(peaks(i)-400:peaks(i)+2000)-f2(peaks(i)-360);
            force1h_3(p7,:)=f3(peaks(i)-400:peaks(i)+2000)-f3(peaks(i)-360);
            force1h_4(p7,:)=f4(peaks(i)-400:peaks(i)+2000)-f4(peaks(i)-360);
            force1h_5(p7,:)=f5(peaks(i)-400:peaks(i)+2000)-f5(peaks(i)-360);
            force1h_6(p7,:)=f6(peaks(i)-400:peaks(i)+2000)-f6(peaks(i)-360);
            weight1h(p7,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
            p1h_plat_pos(p7,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
            p1h_foot_pos(p7,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
            phase1h(p7,:)=ramp(peaks(i)-400:peaks(i)+2000);
            
            filtered_weighth(p7,:)=filtered_weight(peaks(i)-400:peaks(i)+2000)-filtered_weight(peaks(i)-360);
            filtered_fyh(p7,:)    = filtered_fy(peaks(i)-400:peaks(i)+2000)-filtered_fy(peaks(i)-360);
            phase_variableh(p7,:) = phase_variable(peaks(i)-400:peaks(i)+2000)-phase_variable(peaks(i)-360);
            
            p7=p7+1;
        end
        
        
    end
    
    fclose all;
    clear live_data Input1 peaks Img;
    
  end
    
end

%True Gait Phase
max_weight = max(max([weight1r, weight1h]));
start_index_vec_r = zeros(1, NUM_TRAINING_TRIALS/2); 
end_index_vec_r = zeros(1, NUM_TRAINING_TRIALS/2);

start_index_vec_h = zeros(1, NUM_TRAINING_TRIALS/2);
end_index_vec_h = zeros(1, NUM_TRAINING_TRIALS/2);


for i = 1:20
    start_index = min(find(weight1r(i,:) > 0.015*max_weight))
    end_index = max(find(weight1r(i,:) > 0.03*max_weight))
    start_index_vec_r(i) = start_index;
    end_index_vec_r(i) = end_index;
end
for i = 1:20
    start_index = min(find(weight1h(i,:) > 0.015*max_weight))
    end_index = max(find(weight1h(i,:) > 0.03*max_weight))
    start_index_vec_h(i) = start_index;
    end_index_vec_h(i) = end_index;
end



figure();
%Since the trigger point isn't the true start of stance, we need to find
%the time difference (in samples) between the true start and the trigger
%point.
stance_start_to_trigger_r = 400 - start_index_vec_r;
stance_start_to_trigger_h = 400 - start_index_vec_h;

stance_phase_length_r = end_index_vec_r-start_index_vec_r;
stance_phase_length_h = end_index_vec_h-start_index_vec_h;
stance_phase_duration_r = stance_phase_length_r * (1/2000);
stance_phase_duration_h = stance_phase_length_h * (1/2000);
output1r = zeros(size(phase1r));
output1h = zeros(size(phase1h));
for i=1:20
%     output1r(i, start_index_vec_r(i):(start_index_vec_r(i) + stance_phase_length_r(i))) = linspace(0, 100.0,  stance_phase_length_r(i)+1);
%     output1h(i, start_index_vec_h(i):(start_index_vec_h(i) + stance_phase_length_h(i))) = linspace(0, 100.0,  stance_phase_length_h(i)+1);

   value1r=max(phase1r(i,:));
   output1r(i,:)=(phase1r(i,:))*100/value1r;
   value1h=max(phase1h(i,:));
   output1h(i,:)=(phase1h(i,:))*100/value1h;
%     value2r=max(phase2r(i,:));
%    output2r(i,:)=(phase2r(i,:))*100/value2r;
%    value2h=max(phase2h(i,:));
%    output2h(i,:)=(phase2h(i,:))*100/value2h;
%    value3r=max(phase3r(i,:));
%    output3r(i,:)=(phase3r(i,:))*100/value3r;
%    value3h=max(phase3h(i,:));
%    output3h(i,:)=(phase3h(i,:))*100/value3h;
end

figure()
plot(output1r')
hold on;
plot(weight1r'/8)
stem(start_index_vec_r, 100*ones(size(start_index_vec_r)))

figure()
plot(output1h')
hold on;
plot(weight1h'/8)
stem(start_index_vec_h, 100*ones(size(start_index_vec_h)))

%%
input1r_f1=0;
input1r_f2=0;
input1r_f3=0;
input1r_f4=0;
input1r_f5=0;
input1r_f6=0;
input1r_w=0;
input1r_p=0;
input1r_foot=0;

input1r_filt_w = 0;
input1r_filt_fy = 0;
input1r_phase_var = 0;
fin_out1r=0;

for i=1:20
    
%     input1r_f1=[input1r_f1;force1r_1(i,:)'];
%     input1r_f2=[input1r_f2;force1r_2(i,:)'];
%     input1r_f3=[input1r_f3;force1r_3(i,:)'];
%     input1r_f4=[input1r_f4;force1r_4(i,:)'];
%     input1r_f5=[input1r_f5;force1r_5(i,:)'];
%     input1r_f6=[input1r_f6;force1r_6(i,:)'];
%     input1r_w=[input1r_w;weight1r(i,:)'];
     input1r_p=[input1r_p;phase1r(i,:)'];
     input1r_foot=[input1r_foot;p1r_foot_pos(i,:)'];
    input1r_filt_w = [input1r_filt_w; filtered_weightr(i,:)'];
    %input1r_filt_fy = [input1r_filt_fy; filtered_fyr(i,:)'];
    %input1r_phase_var = [input1r_phase_var; phase_variabler(i,:)'];

    
    fin_out1r=[fin_out1r;output1r(i,:)'];
end 



input1h_f1=0;
input1h_f2=0;
input1h_f3=0;
input1h_f4=0;
input1h_f5=0;
input1h_f6=0;
input1h_w=0;
input1h_p=0;
input1h_foot=0;

input1h_filt_w = 0;
input1h_filt_fy = 0;
input1h_phase_var = 0;
fin_out1h=0;
 for i=1:20
    
%     input1h_f1=[input1h_f1;force1h_1(i,:)'];
%     input1h_f2=[input1h_f2;force1h_2(i,:)'];
%     input1h_f3=[input1h_f3;force1h_3(i,:)'];
%     input1h_f4=[input1h_f4;force1h_4(i,:)'];
%     input1h_f5=[input1h_f5;force1h_5(i,:)'];
%     input1h_f6=[input1h_f6;force1h_6(i,:)'];
%     input1h_w=[input1h_w;weight1h(i,:)'];
     input1h_p=[input1h_p;phase1h(i,:)'];
     input1h_foot=[input1h_foot;p1h_foot_pos(i,:)'];
    
    input1h_filt_w = [input1h_filt_w; filtered_weighth(i,:)'];
    %input1h_filt_fy = [input1h_filt_fy; filtered_fyh(i,:)'];
    %input1h_phase_var = [input1h_phase_var; phase_variableh(i,:)'];
    
    fin_out1h=[fin_out1h;output1h(i,:)'];
 end 

%fin_in1r=[input1r_f1,input1r_f2,input1r_f3,input1r_f4,input1r_f5,input1r_f6,input1r_w,input1r_p];

% fin_in1h=[input1h_f1,input1h_f2,input1h_f3,input1h_f4,input1h_f5,input1h_f6,input1h_w,input1h_p]; 
 




%% Total Training %%
[min_grf,max_grf] = bounds([input1r_filt_w; input1h_filt_w],  'all');
[min_pos,max_pos] = bounds([input1r_foot; input1h_foot],  'all');
[min_phas,max_phas] = bounds([input1r_p; input1h_p], 'all');
input1r_foot = filter(num,den,input1r_foot);
input1h_foot = filter(num,den,input1h_foot);
fin_in1r = [(input1r_filt_w-min(input1r_filt_w))/(max(input1r_filt_w)-min(input1r_filt_w)), (input1r_foot -min(input1r_foot))/(max(input1r_foot)-min(input1r_foot)), (input1r_p -min(input1r_p))/(max(input1r_p)-min(input1r_p))];
fin_in1h = [(input1h_filt_w-min(input1h_filt_w))/(max(input1h_filt_w)-min(input1h_filt_w)), (input1h_foot -min(input1h_foot))/(max(input1h_foot)-min(input1h_foot)), (input1h_p -min(input1h_p))/(max(input1h_p)-min(input1h_p))];


total_training_data = [ fin_in1r; fin_in1h];
total_training_output = [fin_out1r; fin_out1h];