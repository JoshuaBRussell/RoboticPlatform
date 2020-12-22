close all
clear all
mvc_evaluation;
% Insert subject initial and name. Make sure it matches the format for naming
sub_initial='T';

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
    foot_pos_data=Input1.data(:,14);
    foot_const=mean(foot_pos_data);
    foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
    plat_pos_data=Input1.data(:,16);
    plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
    ramp=Input1.data(:,22);
    [test,peaks]=findpeaks(Input1.data(:,17));
    
    for i=1:length(peaks)
        
        time=[-200:0.5:1000];
       
        if test(i)==4
            force1r_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1r_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            
            weight1r(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
            
            p1r_plat_pos(p0,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
            p1r_foot_pos(p0,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
            phase1r(p0,:)=ramp(peaks(i)-400:peaks(i)+2000);
            p0=p0+1;
        end
        
       
       
        if test(i)==8
            force1h_1(p7,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1h_2(p7,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1h_3(p7,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1h_4(p7,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1h_5(p7,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            force1h_6(p7,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
            weight1h(p7,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
            p1h_plat_pos(p7,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
            p1h_foot_pos(p7,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
            phase1h(p7,:)=ramp(peaks(i)-400:peaks(i)+2000);
            p7=p7+1;
        end
        
        
    end
    
    fclose all;
    clear live_data Input1 peaks Img;
    
  end
    
end

% p0=1;
% p7=1;
% for trials=1:4
%     
%   if(ismember(trials,exclude)==0)  
%     if(trials<10)
%         h = fopen(strcat(sub_initial,'PW2',num2str(trials),'.dat'));
%     end
%     
%     live_data=fread(h);
%     Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
%     siz=size(Input1.data);
%     
%     %             %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','CWamplerate',2000);
%     %d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 20, 'CWampleRate', 2000, 'DesignMethod', 'window');
%     pert_torque=filtfilt(d1,Input1.data(:,7));
%     f1=Input1.data(:,9)*53.4;
%     
%     f2=Input1.data(:,10)*53.4;
% 
%     f3=Input1.data(:,11)*53.4;
% 
%     f4=Input1.data(:,12)*53.4;
%  
%     f5=Input1.data(:,20)*53.4/2;
%   
%     f6=Input1.data(:,21)*53.4/2;
%  
%   
%     w1=filtfilt(d1,Input1.data(:,18));
%   
%     flag=Input1.data(:,17);
%     foot_pos_data=Input1.data(:,14);
%     foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
%     plat_pos_data=Input1.data(:,16);
%     plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
%     ramp=Input1.data(:,22);
%     [test,peaks]=findpeaks(Input1.data(:,17));
%     
%     for i=1:length(peaks)
%         
%         time=[-200:0.5:1000];
%        
%         if test(i)==4
%             force2r_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2r_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2r_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2r_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2r_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2r_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             
%             weight2r(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
%             
%             p2r_plat_pos(p0,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
%             p2r_foot_pos(p0,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
%             phase2r(p0,:)=ramp(peaks(i)-400:peaks(i)+2000);
%             p0=p0+1;
%         end
%         
%        
%        
%         if test(i)==8
%             force2h_1(p7,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2h_2(p7,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2h_3(p7,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2h_4(p7,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2h_5(p7,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force2h_6(p7,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             weight2h(p7,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
%             p2h_plat_pos(p7,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
%             p2h_foot_pos(p7,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
%             phase2h(p7,:)=ramp(peaks(i)-400:peaks(i)+2000);
%             p7=p7+1;
%         end
%         
%         
%     end
%     
%     fclose all;
%     clear live_data Input1 peaks Img;
%     
%   end
%     
% end
% 
% p0=1;
% p7=1;
% for trials=1:4
%     
%   if(ismember(trials,exclude)==0)  
%     if(trials<10)
%         h = fopen(strcat(sub_initial,'PW1',num2str(trials),'.dat'));
%     end
%     
%     live_data=fread(h);
%     Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
%     siz=size(Input1.data);
%     
%     %             %d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',20,'DesignMethod','butter','CWamplerate',2000);
%     %d1 = designfilt('lowpassfir', 'FilterOrder', 50, 'CutoffFrequency', 20, 'CWampleRate', 2000, 'DesignMethod', 'window');
%     pert_torque=filtfilt(d1,Input1.data(:,7));
%     f1=Input1.data(:,9)*53.4;
%     
%     f2=Input1.data(:,10)*53.4;
% 
%     f3=Input1.data(:,11)*53.4;
% 
%     f4=Input1.data(:,12)*53.4;
%  
%     f5=Input1.data(:,20)*53.4/2;
%   
%     f6=Input1.data(:,21)*53.4/2;
%  
%   
%     w1=filtfilt(d1,Input1.data(:,18));
%   
%     flag=Input1.data(:,17);
%     foot_pos_data=Input1.data(:,14);
%     foot_pos_data=((foot_pos_data-mean(foot_pos_data))*DP_foot_gonio*pi/180);
%     plat_pos_data=Input1.data(:,16);
%     plat_pos_data=((plat_pos_data-mean(plat_pos_data))*DP_plat_gonio*pi/180);
%     ramp=Input1.data(:,22);
%     [test,peaks]=findpeaks(Input1.data(:,17));
%     
%     for i=1:length(peaks)
%         
%         time=[-200:0.5:1000];
%        
%         if test(i)==4
%             force3r_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3r_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3r_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3r_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3r_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3r_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             
%             weight3r(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
%             
%             p3r_plat_pos(p0,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
%             p3r_foot_pos(p0,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
%             phase3r(p0,:)=ramp(peaks(i)-400:peaks(i)+2000);
%             p0=p0+1;
%         end
%         
%        
%        
%         if test(i)==8
%             force3h_1(p7,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3h_2(p7,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3h_3(p7,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3h_4(p7,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3h_5(p7,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             force3h_6(p7,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
%             weight3h(p7,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
%             p3h_plat_pos(p7,:)=plat_pos_data(peaks(i)-400:peaks(i)+2000);
%             p3h_foot_pos(p7,:)=foot_pos_data(peaks(i)-400:peaks(i)+2000);
%             phase3h(p7,:)=ramp(peaks(i)-400:peaks(i)+2000);
%             p7=p7+1;
%         end
%         
%         
%     end
%     
%     fclose all;
%     clear live_data Input1 peaks Img;
%     
%   end
%     
% end

for i=1:20
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
fin_out1r=0;

for i=1:20
    
    input1r_f1=[input1r_f1;force1r_1(i,:)'];
    input1r_f2=[input1r_f2;force1r_2(i,:)'];
    input1r_f3=[input1r_f3;force1r_3(i,:)'];
    input1r_f4=[input1r_f4;force1r_4(i,:)'];
    input1r_f5=[input1r_f5;force1r_5(i,:)'];
    input1r_f6=[input1r_f6;force1r_6(i,:)'];
    input1r_w=[input1r_w;weight1r(i,:)'];
    input1r_p=[input1r_p;phase1r(i,:)'];
    input1r_foot=[input1r_foot;p1r_foot_pos(i,:)'];
    fin_out1r=[fin_out1r;output1r(i,:)'];
end 
 fin_in1r=[input1r_f1,input1r_f2,input1r_f3,input1r_f4,input1r_f5,input1r_f6,input1r_w,input1r_p];     

 
%  input2r_f1=0;
% input2r_f2=0;
% input2r_f3=0;
% input2r_f4=0;
% input2r_f5=0;
% input2r_f6=0;
% input2r_w=0;
% input2r_p=0;
% input2r_foot=0;
% fin_out2r=0;
%  for i=1:40
%     
%     input2r_f1=[input2r_f1;force2r_1(i,:)'];
%     input2r_f2=[input2r_f2;force2r_2(i,:)'];
%     input2r_f3=[input2r_f3;force2r_3(i,:)'];
%     input2r_f4=[input2r_f4;force2r_4(i,:)'];
%     input2r_f5=[input2r_f5;force2r_5(i,:)'];
%     input2r_f6=[input2r_f6;force2r_6(i,:)'];
%     input2r_w=[input2r_w;weight2r(i,:)'];
%     input2r_p=[input2r_p;phase2r(i,:)'];
%     input2r_foot=[input2r_foot;p3h_foot_pos(i,:)'];
%     fin_out2r=[fin_out2r;output2r(i,:)'];
% end 
%  fin_in2r=[input2r_f1,input2r_f2,input2r_f3,input2r_f4,input2r_f5,input2r_f6,input2r_w,input2r_p,input2r_foot];     
% 
%  input3r_f1=0;
% input3r_f2=0;
% input3r_f3=0;
% input3r_f4=0;
% input3r_f5=0;
% input3r_f6=0;
% input3r_w=0;
% input3r_p=0;
% input3r_foot=0;
% fin_out3r=0;
%  for i=1:40
%     
%     input3r_f1=[input3r_f1;force3r_1(i,:)'];
%     input3r_f2=[input3r_f2;force3r_2(i,:)'];
%     input3r_f3=[input3r_f3;force3r_3(i,:)'];
%     input3r_f4=[input3r_f4;force3r_4(i,:)'];
%     input3r_f5=[input3r_f5;force3r_5(i,:)'];
%     input3r_f6=[input3r_f6;force3r_6(i,:)'];
%     input3r_w=[input3r_w;weight3r(i,:)'];
%     input3r_p=[input3r_p;phase3r(i,:)'];
%     input3r_foot=[input3r_foot;p3h_foot_pos(i,:)'];
%     fin_out3r=[fin_out3r;output3r(i,:)'];
% end 
%  fin_in3r=[input3r_f1,input3r_f2,input3r_f3,input3r_f4,input3r_f5,input3r_f6,input3r_w,input3r_p,input3r_foot];     
input1h_f1=0;
input1h_f2=0;
input1h_f3=0;
input1h_f4=0;
input1h_f5=0;
input1h_f6=0;
input1h_w=0;
input1h_p=0;
input1h_foot=0;
fin_out1h=0;
 for i=1:20
    
    input1h_f1=[input1h_f1;force1h_1(i,:)'];
    input1h_f2=[input1h_f2;force1h_2(i,:)'];
    input1h_f3=[input1h_f3;force1h_3(i,:)'];
    input1h_f4=[input1h_f4;force1h_4(i,:)'];
    input1h_f5=[input1h_f5;force1h_5(i,:)'];
    input1h_f6=[input1h_f6;force1h_6(i,:)'];
    input1h_w=[input1h_w;weight1h(i,:)'];
    input1h_p=[input1h_p;phase1h(i,:)'];
    input1h_foot=[input1h_foot;p1h_foot_pos(i,:)'];
    fin_out1h=[fin_out1h;output1h(i,:)'];
end 
 fin_in1h=[input1h_f1,input1h_f2,input1h_f3,input1h_f4,input1h_f5,input1h_f6,input1h_w,input1h_p];     
% input2h_f1=0;
% input2h_f2=0;
% input2h_f3=0;
% input2h_f4=0;
% input2h_f5=0;
% input2h_f6=0;
% input2h_w=0;
% input2h_p=0;
% input2h_foot=0;
% fin_out2h=0;
%  for i=1:40
%     
%     input2h_f1=[input2h_f1;force2h_1(i,:)'];
%     input2h_f2=[input2h_f2;force2h_2(i,:)'];
%     input2h_f3=[input2h_f3;force2h_3(i,:)'];
%     input2h_f4=[input2h_f4;force2h_4(i,:)'];
%     input2h_f5=[input2h_f5;force2h_5(i,:)'];
%     input2h_f6=[input2h_f6;force2h_6(i,:)'];
%     input2h_w=[input2h_w;weight2h(i,:)'];
%     input2h_p=[input2h_p;phase2h(i,:)'];
%     input2h_foot=[input2h_foot;p3h_foot_pos(i,:)'];
%     fin_out2h=[fin_out2h;output2h(i,:)'];
% end 
%  fin_in2h=[input2h_f1,input2h_f2,input2h_f3,input2h_f4,input2h_f5,input2h_f6,input2h_w,input2h_p,input2h_foot];     
% input3h_f1=0;
% input3h_f2=0;
% input3h_f3=0;
% input3h_f4=0;
% input3h_f5=0;
% input3h_f6=0;
% input3h_w=0;
% input3h_p=0;
% input3h_foot=0;
% fin_out3h=0;
%  for i=1:40
%     
%     input3h_f1=[input3h_f1;force3h_1(i,:)'];
%     input3h_f2=[input3h_f2;force3h_2(i,:)'];
%     input3h_f3=[input3h_f3;force3h_3(i,:)'];
%     input3h_f4=[input3h_f4;force3h_4(i,:)'];
%     input3h_f5=[input3h_f5;force3h_5(i,:)'];
%     input3h_f6=[input3h_f6;force3h_6(i,:)'];
%     input3h_w=[input3h_w;weight3h(i,:)'];
%     input3h_p=[input3h_p;phase3h(i,:)'];
%     input3h_foot=[input3h_foot;p3h_foot_pos(i,:)'];
%     fin_out3h=[fin_out3h;output3h(i,:)'];
% end 
%  fin_in3h=[input3h_f1,input3h_f2,input3h_f3,input3h_f4,input3h_f5,input3h_f6,input3h_w,input3h_p,input3h_foot];     
% 
%  fin_inr=[fin_in1r;fin_in2r;fin_in3r];
%  fin_outr=[fin_out1r;fin_out2r;fin_out3r];
%  
%   fin_inh=[fin_in1h;fin_in2h;fin_in3h];
%  fin_outh=[fin_out1h;fin_out2h;fin_out3h];