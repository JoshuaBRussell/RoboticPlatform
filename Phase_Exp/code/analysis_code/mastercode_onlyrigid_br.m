
%----mastercode_onlyrigid_br.m----%
% Same as mastercode_onlyrigid.m except no bootstrapping. Just basic
% regression with outlier removal

%% Setup, Filtering and Data Visualization settings
close all
clear all
% Getting the MVC values
mvc_evaluation;
% Insert subject initial and name.
%Make sure it matches the format for naming
sub_initial='M';
sub_name='Morgan';
%Add number of perturbation you actually ran
num_pert=40;
% insert lower limit of inertia of foot in the fit
% u_lim is the upper limit of the inertia and lim
% is the lower limit
lim=0.007;
u_lim=0.2;
% change these flags to 1 for figures (normal fit and constrained fit)
plot_figs=1;
% plot_figs_constrained=0;
% Change to 1 to get ankle position histogram
plot_hist=1;
%  Change to get torque plot comparison figure
plot_torque=1;
shift=0;
%add the number of loops you want to run bootstrapping here
loops=100;
% Add the trials you want to exclude in here
exclude=[];%exclude=[1,2];
d3 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',...
    5,'DesignMethod','butter','Samplerate',2000);
d1 = designfilt('lowpassiir','FilterOrder',4,'HalfPowerFrequency',15,...
    'DesignMethod','butter','Samplerate',2000);
%% Section to calculate goniometer gains
t=gonio_values_func;
DP_foot_gonio=t(1);
DP_plat_gonio=t(2);

%close all    %uncomment to close all figures
%% Acquisition of trial data from different files
i=0;
p1=1;
p2=1;
p3=1;
p0=1;
p4=1;
p5=1;
p6=1;
p7=1;
% Make sure you only add the total number of
% trials you actually run
for trials=1:10
    
    if(ismember(trials,exclude)==0)
        if(trials<10)
            h = fopen(strcat(sub_initial,'W0',num2str(trials),'.dat'));
        else
            h = fopen(strcat(sub_initial,'W',num2str(trials),'.dat'));
        end
        
        live_data=fread(h);
        Input1= SimulinkRealTime.utils.getFileScopeData(live_data);
        siz=size(Input1.data);
        Img_flag=Input1.data(:,6);
        [x,img_st]=findpeaks(diff(Img_flag));
        img_st=round(img_st(1)/20);
        if(trials<10)
            Img=csvread(strcat(sub_name,'_00',num2str(trials),'.csv'));
        else
            Img=csvread(strcat(sub_name,'_0',num2str(trials),'.csv'));
        end
        % obtaining data from channels
        pert_torque=filtfilt(d1,Input1.data(:,7));
        f1=Input1.data(:,9)*53.4;
        f1=filtfilt(d1,f1);
        f2=Input1.data(:,10)*53.4;
        f2=filtfilt(d1,f2);
        f3=Input1.data(:,11)*53.4;
        f3=filtfilt(d1,f3);
        f4=Input1.data(:,12)*53.4;
        f4=filtfilt(d1,f4);
        f5=Input1.data(:,20)*53.4/2;
        f5=filtfilt(d1,f5);
        f6=Input1.data(:,21)*53.4/2;
        f6=filtfilt(d1,f6);
        ta=Input1.data(:,1);
        ta=abs(ta-off_TA)*100/mvc_ta;
        sol=Input1.data(:,2);
        sol=abs(sol-off_SOL)*100/mvc_sol;
        pl=Input1.data(:,3);
        pl=abs(pl-off_PL)*100/mvc_pl;
        gca=Input1.data(:,4);
        gca=abs(gca-off_GCA)*100/mvc_gca;
        w1=filtfilt(d1,Input1.data(:,18));
        cop=filtfilt(d1,Input1.data(:,19));
        flag=Input1.data(:,17);
        rigid_phase_tot=Input1.data(:,15);
        haptic_phase_tot=Input1.data(:,16);
        perturb_start=Input1.data(:,22);
        foot_pos_data=filtfilt(d1,Input1.data(:,13));
        foot_pos_data=((foot_pos_data-mean(foot_pos_data))...
            *DP_foot_gonio*pi/180);
        plat_pos_data=filtfilt(d1,Input1.data(:,14));
        plat_pos_data=((plat_pos_data-mean(plat_pos_data))...
            *DP_plat_gonio*pi/180);
        % 17 records an impulse everytime a
        % perturbation occurs with diff amplitudes
        [test,peaks]=findpeaks(Input1.data(:,17));
        for i=1:length(peaks)
     % Each test value corresponds to a different pert       
            time=[-200:0.5:1000];
            if test(i)==1
                force1_1(p1,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_2(p1,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_3(p1,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_4(p1,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_5(p1,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force1_6(p1,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                weight1(p1,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p1_plat_torque(p1,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p1_plat_pos(p1,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p1_foot_pos(p1,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p1_phase(p1,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p1_pert(p1,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img1_pos(p1)=getmin(peaks(i),img_st,Img);
                cop1(p1,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p1_pert(p1,:)));
                p1_peakst(p1)=b;
                p1=p1+1;
                
            end
            if test(i)==2
                force2_1(p2,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_2(p2,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_3(p2,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_4(p2,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_5(p2,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force2_6(p2,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                
                weight2(p2,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p2_plat_torque(p2,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p2_plat_pos(p2,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p2_foot_pos(p2,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p2_phase(p2,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p2_pert(p2,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img2_pos(p2)=getmin(peaks(i),img_st,Img);
                cop2(p2,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p2_pert(p2,:)));

                p2_peakst(p2)=b;
                p2=p2+1;
            end
            if test(i)==3
                force3_1(p3,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_2(p3,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_3(p3,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_4(p3,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_5(p3,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force3_6(p3,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                
                weight3(p3,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p3_plat_torque(p3,:)=1*(pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50));
                p3_plat_pos(p3,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p3_foot_pos(p3,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p3_phase(p3,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p3_pert(p3,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img3_pos(p3)=getmin(peaks(i),img_st,Img);
                cop3(p3,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p3_pert(p3,:)));
                p3_peakst(p3)=b;
                p3=p3+1;
            end
            if test(i)==4
                force0_1(p0,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_2(p0,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_3(p0,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_4(p0,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_5(p0,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force0_6(p0,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                
                ta_emg(p0,:)=ta(peaks(i)-400:peaks(i)+2000);
                sol_emg(p0,:)=sol(peaks(i)-400:peaks(i)+2000);
                pl_emg(p0,:)=pl(peaks(i)-400:peaks(i)+2000);
                gca_emg(p0,:)=gca(peaks(i)-400:peaks(i)+2000);
                weight4(p0,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p0_plat_torque(p0,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p0_plat_pos(p0,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p0_foot_pos(p0,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p0_phase(p0,:)=rigid_phase_tot(peaks(i)-400:peaks(i)+2000);
                p0_pert(p0,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img0_pos(p0)=getmin(peaks(i),img_st,Img);
                cop4(p0,:)=cop(peaks(i)-400:peaks(i)+2000);
                [a,b]=findpeaks(diff(p0_pert(p0,:)));

                p0_peakst(p0)=400;%b;

                [a,b]=min(diff(p0_phase(p0,:)));

                p0_peakend(p0)=b;
                p0=p0+1;
            end
            if test(i)==5
                force4_1(p4,:)=f1(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_2(p4,:)=f2(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_3(p4,:)=f3(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_4(p4,:)=f4(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_5(p4,:)=f5(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                force4_6(p4,:)=f6(peaks(i)-400:peaks(i)+2000)-f1(peaks(i)-360);
                weight3(p4,:)=w1(peaks(i)-400:peaks(i)+2000)-w1(peaks(i)-360);%+w2(peaks(i)-400:peaks(i)+2000)-w2(peaks(i)-360)+w3(peaks(i)-400:peaks(i)+2000)-w3(peaks(i)-360)+w4(peaks(i)-400:peaks(i)+2000)-w4(peaks(i)-20);
                p4_plat_torque(p4,:)=pert_torque(peaks(i)-400:peaks(i)+2000)-pert_torque(peaks(i)+50);
                p4_plat_pos(p4,:)=plat_pos_data(peaks(i)-179:peaks(i)+2221);
                p4_foot_pos(p4,:)=foot_pos_data(peaks(i)-179+shift:peaks(i)+2221+shift);
                p4_phase(p4,:)=haptic_phase_tot(peaks(i)-400:peaks(i)+2000);
                p4_pert(p4,:)=perturb_start(peaks(i)-400:peaks(i)+2000);
                img4_pos(p4)=getmin(peaks(i),img_st,Img);
                [a,b]=findpeaks(diff(p4_pert(p4,:)));

                p4_peakst(p4)=b;
        
                p4=p4+1;
            end
            
        end
        
%         fclose all;
%         clear live_data Input1 peaks Img;
        
    end
    
end
p0_er=0;
p1_er=0;
p2_er=0;
p3_er=0;
p4_er=0;


p0_er=0;
p1_er=0;
p2_er=0;
p3_er=0;
p4_er=0;

%% removing outliers using motion capture data: Its redundant but I left it as it wont really change anything
for i=1:p1-1
    if(img1_pos(i)>500)
        img1_pos(i)=NaN;
    end
end
for i=1:p2-1
    if(img2_pos(i)>500)
        img2_pos(i)=NaN;
    end
end
for i=1:p3-1
    if(img3_pos(i)>500)
        img3_pos(i)=NaN;
    end
end
for i=1:p0-1
    if(img0_pos(i)>500)
        img0_pos(i)=NaN;
    end
end


if plot_hist==1
    im=[img0_pos img1_pos img2_pos img3_pos img4_pos];
end

p0_removed_ind = abs(img0_pos) > 2.0; 

%% Obtaining perturbation torque and CoP data from different trials the act_torque is calculated by using the foot position obtained form motion capture
for i=1:p1-1
    if (isnan(p1_plat_torque(i,1))==1)
        p1_act_torque(i,1:2201)=NaN;
        ssss='Found one';
        disp(ssss);
    else
        p1_act_torque(i,:)=force1_1(i,:)*(0.315-(img1_pos(i)/100))+force1_4(i,:)*(0.315-(img1_pos(i)/100))-force1_2(i,:)*(0.105+(img1_pos(i)/100))-force1_3(i,:)*(0.105+(img1_pos(i)/100))-force1_5(i,:)*(0.095)-force1_6(i,:)*(0.095);
    end
end
for i=1:p2-1
    if (isnan(p2_plat_torque(i,1))==1)
        p2_act_torque(i,1:2201)=NaN;
    else
        p2_act_torque(i,:)=force2_1(i,:)*(0.315-(img2_pos(i)/100))+force2_4(i,:)*(0.315-(img2_pos(i)/100))-force2_2(i,:)*(0.105+(img2_pos(i)/100))-force2_3(i,:)*(0.105+(img2_pos(i)/100))-force2_5(i,:)*(0.095)-force2_6(i,:)*(0.095);
    end
end


for i=1:p3-1
    if (isnan(p3_plat_torque(i,1))==1)
        p3_act_torque(i,1:2201)=NaN;
    else
        p3_act_torque(i,:)=force3_1(i,:)*(0.315-(img3_pos(i)/100))+force3_4(i,:)*(0.315-(img3_pos(i)/100))-force3_2(i,:)*(0.105+(img3_pos(i)/100))-force3_3(i,:)*(0.105+(img3_pos(i)/100))-force3_5(i,:)*(0.095)-force3_6(i,:)*(0.095);
    end
end
for i=1:p0-1
    if (isnan(p0_plat_torque(i,1))==1)
        p0_act_torque(i,1:2201)=NaN;
        p0_cop_torque(i,1:2201)=NaN;
    else
        p0_act_torque(i,:)=force0_1(i,:)*(0.315-(img0_pos(i)/100))+force0_4(i,:)*(0.315-(img0_pos(i)/100))-force0_2(i,:)*(0.105+(img0_pos(i)/100))-force0_3(i,:)*(0.105+(img0_pos(i)/100))-force0_5(i,:)*(0.095)-force0_6(i,:)*(0.095);
        p0_cop_torque(i,:)=force0_1(i,:)*(0.315-(img0_pos(i)/100))+force0_4(i,:)*(0.315-(img0_pos(i)/100))-force0_2(i,:)*(0.105+(img0_pos(i)/100))-force0_3(i,:)*(0.105+(img0_pos(i)/100))-force0_5(i,:)*(0.025)-force0_6(i,:)*(0.025);
        
    end
end
for i=1:p4-1
    if (isnan(p4_plat_torque(i,1))==1)
        p4_act_torque(i,1:2201)=NaN;
    else
        p4_act_torque(i,:)=force4_1(i,:)*(0.315-(img4_pos(i)/100))+force4_4(i,:)*(0.315-(img4_pos(i)/100))-force4_2(i,:)*(0.105+(img4_pos(i)/100))-force4_3(i,:)*(0.105+(img4_pos(i)/100))-force4_5(i,:)*(0.095)-force4_6(i,:)*(0.095);
    end
end




%% Calculation of means for plotting...not really used in analysis
% weight1m=nanmean(weight1);
% weight2m=nanmean(weight2);
% weight3m=nanmean(weight3);
% weight4m=trimmean(weight4,30);
% 
% cop4m=trimmean(p0_cop_torque,30);
% 
p0_plat_torquem=trimmean(p0_act_torque,30);
p0_plat_posm=trimmean(p0_plat_pos,30);
p0_foot_posm=nanmean(p0_foot_pos);
p1_plat_torquem=trimmean(p1_act_torque,30);
p1_plat_posm=trimmean(p1_plat_pos,30);
p1_foot_posm=trimmean(p1_foot_pos,30);
p2_plat_torquem=trimmean(p2_act_torque,30);
p2_plat_posm=trimmean(p2_plat_pos,30);
p2_foot_posm=trimmean(p2_foot_pos,30);
p3_plat_torquem=trimmean(p3_act_torque,30);
p3_plat_posm=trimmean(p3_plat_pos,30);
p3_foot_posm=trimmean(p3_foot_pos,30);
p4_plat_torquem=trimmean(p4_act_torque,30);
p4_plat_posm=trimmean(p4_plat_pos,30);
p4_foot_posm=trimmean(p4_foot_pos,30);


for i=1:p0-1
    
    ta_emg(i,:)=filtfilt(d3,ta_emg(i,:));
    pl_emg(i,:)=filtfilt(d3,pl_emg(i,:));
    sol_emg(i,:)=filtfilt(d3,sol_emg(i,:));
    gca_emg(i,:)=filtfilt(d3,gca_emg(i,:));
    
end



ta_emgm=trimmean(ta_emg,30);
sol_emgm=trimmean(sol_emg,30);
pl_emgm=trimmean(pl_emg,30);
gca_emgm=trimmean(gca_emg,30);

%% Rigid and Haptic no perturbation torque;
% the torque is collected from 100 samples prior
% to 300 samples post perturbation (heel strike)
%It is then subtracted from observed torque to get
%differential torque.
%Excluded uses min number of perturbations so that
%analysis is unifrom
excluded=[p1,p2,p3,p4];
analysis_value=min(excluded);
for i=1:p0-1
    point_i=floor((p0_peakend(i)-p0_peakst(i))*0.18)+p0_peakst(i)
    p0_plat_torque15(i,:)=p0_act_torque(i,point_i-100:point_i+300);
    p0_cop_torque15(i,:)=p0_cop_torque(i,point_i-100:point_i+300);
    p0_foot_pos15(i,:)=p0_foot_pos(i,point_i-100:point_i+300);
    
    point_j=floor((p0_peakend(i)-p0_peakst(i))*0.31)+p0_peakst(i)
    p0_plat_torque30(i,:)=p0_act_torque(i,point_j-100:point_j+300);
    p0_cop_torque30(i,:)=p0_cop_torque(i,point_j-100:point_j+300);
    p0_foot_pos30(i,:)=p0_foot_pos(i,point_j-100:point_j+300);
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*0.44)+p0_peakst(i)
    p0_plat_torque45(i,:)=p0_act_torque(i,point_k-100:point_k+300);
    p0_cop_torque45(i,:)=p0_cop_torque(i,point_k-100:point_k+300);
    p0_foot_pos45(i,:)=p0_foot_pos(i,point_k-100:point_k+300);
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*0.57)+p0_peakst(i)
    p0_plat_torque60(i,:)=p0_act_torque(i,point_k-100:point_k+300);
    p0_cop_torque60(i,:)=p0_cop_torque(i,point_k-100:point_k+300);
    p0_foot_pos60(i,:)=p0_foot_pos(i,point_k-100:point_k+300);
end

% Weight during perturbation
for i=1:p0-1
    point_i=floor((p0_peakend(i)-p0_peakst(i))*0.18)+p0_peakst(i)
    weightr15(i,:)=weight4(i,point_i-100:point_i+300);
    point_j=floor((p0_peakend(i)-p0_peakst(i))*0.31)+p0_peakst(i)
    weightr30(i,:)=weight4(i,point_j-100:point_j+300);
    point_k=floor((p0_peakend(i)-p0_peakst(i))*0.44)+p0_peakst(i)
    weightr45(i,:)=weight4(i,point_k-100:point_k+300);
    point_k=floor((p0_peakend(i)-p0_peakst(i))*0.57)+p0_peakst(i)
    weightr60(i,:)=weight4(i,point_k-100:point_k+300);
    
end


%CoP during perturbation
for i=1:p0-1
 for j=1:length(weightr15)
     copr15(i,j)=p0_cop_torque15(i,j)/weightr15(i,j);
 end
for j=1:length(weightr30)
     copr30(i,j)=p0_cop_torque30(i,j)/weightr30(i,j);
end
for j=1:length(weightr45)
     copr45(i,j)=p0_cop_torque45(i,j)/weightr45(i,j);
end
 for j=1:length(weightr60)
     copr60(i,j)=p0_cop_torque60(i,j)/weightr60(i,j);
 end
    
end


for i=1:p0-1
    point_i=floor((p0_peakend(i)-p0_peakst(i))*0.18)+p0_peakst(i)
    ta15(i)=mean(ta_emg(i,point_i-2:point_i+2));
    pl15(i)=mean(pl_emg(i,point_i-25:point_i+25));
    sol15(i)=mean(sol_emg(i,point_i-25:point_i+25));
    gca15(i)=mean(gca_emg(i,point_i-25:point_i+25));
   
    point_j=floor((p0_peakend(i)-p0_peakst(i))*0.31)+p0_peakst(i)
    ta30(i)=mean(ta_emg(i,point_j-25:point_j+25));
    pl30(i)=mean(pl_emg(i,point_j-25:point_j+25));
    sol30(i)=mean(sol_emg(i,point_j-25:point_j+25));
    gca30(i)=mean(gca_emg(i,point_j-25:point_j+25));
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*0.44)+p0_peakst(i)
    ta45(i)=mean(ta_emg(i,point_k-25:point_k+25));
    pl45(i)=mean(pl_emg(i,point_k-25:point_k+25));
    sol45(i)=mean(sol_emg(i,point_k-25:point_k+25));
    gca45(i)=mean(gca_emg(i,point_k-25:point_k+25));
    
    point_k=floor((p0_peakend(i)-p0_peakst(i))*0.57)+p0_peakst(i)
    ta60(i)=mean(ta_emg(i,point_k-25:point_k+25));
    pl60(i)=mean(pl_emg(i,point_k-25:point_k+25));
    sol60(i)=mean(sol_emg(i,point_k-25:point_k+25));
    gca60(i)=mean(gca_emg(i,point_k-25:point_k+25));
end

p0_plat_torque15m=mean(p0_plat_torque15);
p0_plat_torque30m=mean(p0_plat_torque30);
p0_plat_torque45m=mean(p0_plat_torque45);
p0_plat_torque60m=mean(p0_plat_torque60);
p0_foot_pos15m=mean(p0_foot_pos15);
p0_foot_pos30m=mean(p0_foot_pos30);
p0_foot_pos45m=mean(p0_foot_pos45);
p0_foot_pos60m=mean(p0_foot_pos60);


ta15m=mean(ta15);
pl15m=mean(pl15);
sol15m=mean(sol15);
gca15m=mean(gca15);
ta30m=mean(ta30);
pl30m=mean(pl30);
sol30m=mean(sol30);
gca30m=mean(gca30);
ta45m=mean(ta45);
pl45m=mean(pl45);
sol45m=mean(sol45);
gca45m=mean(gca45);
ta60m=mean(ta60);
pl60m=mean(pl60);
sol60m=mean(sol60);
gca60m=mean(gca60);
emg_final=[ta15m,pl15m,sol15m,gca15m;ta30m,pl30m,sol30m,gca30m;ta45m,pl45m,sol45m,gca45m;ta60m,pl60m,sol60m,gca60m;];
%% Obtaining differential data for impedance analysis
% data is broken into chunks of 400 points:100
% before pert and 300 after pert. 

for i=1:analysis_value-1
    p1_foot_pos(i,:)=p1_foot_pos(i,:)-p0_foot_posm;
    p2_foot_pos(i,:)=p2_foot_pos(i,:)-p0_foot_posm;
    p3_foot_pos(i,:)=p3_foot_pos(i,:)-p0_foot_posm;
    p4_foot_pos(i,:)=p4_foot_pos(i,:)-p0_foot_posm;
    
    diff_p1_plat_pos(i,:)=p1_plat_pos(i,p1_peakst(i)-100:p1_peakst(i)+300);
    diff_p2_plat_pos(i,:)=p2_plat_pos(i,p2_peakst(i)-100:p2_peakst(i)+300);
    diff_p3_plat_pos(i,:)=p3_plat_pos(i,p3_peakst(i)-100:p3_peakst(i)+300);
    diff_p1_plat_torque(i,:)=p1_act_torque(i,p1_peakst(i)-100:p1_peakst(i)+300);
    diff_p2_plat_torque(i,:)=p2_act_torque(i,p2_peakst(i)-100:p2_peakst(i)+300);
    diff_p3_plat_torque(i,:)=p3_act_torque(i,p3_peakst(i)-100:p3_peakst(i)+300);
    
    diff_p1_plat_torque(i,:)=diff_p1_plat_torque(i,:)-p0_plat_torque30m;
    diff_p2_plat_torque(i,:)=diff_p2_plat_torque(i,:)-p0_plat_torque45m;
    diff_p3_plat_torque(i,:)=diff_p3_plat_torque(i,:)-p0_plat_torque60m;
    
%     diff_p1_foot_pos(i,:)=p1_foot_pos(i,p1_peakst(i)-100:p1_peakst(i)+300)-p0_foot_pos30m;
%     diff_p2_foot_pos(i,:)=p2_foot_pos(i,p2_peakst(i)-100:p2_peakst(i)+300)-p0_foot_pos45m;
%     diff_p3_foot_pos(i,:)=p3_foot_pos(i,p3_peakst(i)-100:p3_peakst(i)+300)-p0_foot_pos60m;

    diff_p1_foot_pos(i,:)=p1_foot_pos(i,p1_peakst(i)-100:p1_peakst(i)+300);
    diff_p2_foot_pos(i,:)=p2_foot_pos(i,p2_peakst(i)-100:p2_peakst(i)+300);
    diff_p3_foot_pos(i,:)=p3_foot_pos(i,p3_peakst(i)-100:p3_peakst(i)+300);
    
%      diff_p1_foot_pos(i,:)=diff_p1_foot_pos(i,:)-p0_foot_pos30m;
%     diff_p2_foot_pos(i,:)=diff_p2_foot_pos(i,:)-p0_foot_pos45m;
%     diff_p3_foot_pos(i,:)=diff_p3_foot_pos(i,:)-p0_foot_pos60m;
    diff_p1_foot_pos(i,:)=diff_p1_foot_pos(i,:)-diff_p1_foot_pos(i,100);
    diff_p2_foot_pos(i,:)=diff_p2_foot_pos(i,:)-diff_p2_foot_pos(i,100);
    diff_p3_foot_pos(i,:)=diff_p3_foot_pos(i,:)-diff_p3_foot_pos(i,100);
   

    
    diff_p1_plat_pos(i,:)=diff_p1_plat_pos(i,:)-diff_p1_plat_pos(i,100);
    diff_p2_plat_pos(i,:)=diff_p2_plat_pos(i,:)-diff_p2_plat_pos(i,100);
    diff_p3_plat_pos(i,:)=diff_p3_plat_pos(i,:)-diff_p3_plat_pos(i,100);
    
    diff_p1_plat_vel(i,1)=0;
    for l=2:length(diff_p1_plat_pos(1,:))
        diff_p1_plat_vel(i,l)=(diff_p1_plat_pos(i,l)-diff_p1_plat_pos(i,l-1))/0.0005;
    end
    diff_p2_plat_vel(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p2_plat_vel(i,l)=(diff_p2_plat_pos(i,l)-diff_p2_plat_pos(i,l-1))/0.0005;
    end
    diff_p3_plat_vel(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p3_plat_vel(i,l)=(diff_p3_plat_pos(i,l)-diff_p3_plat_pos(i,l-1))/0.0005;
    end
    diff_p1_foot_vel(i,1)=0;
    for l=2:length(diff_p1_foot_pos(1,:))
        diff_p1_foot_vel(i,l)=(diff_p1_foot_pos(i,l)-diff_p1_foot_pos(i,l-1))/0.0005;
    end
    diff_p2_foot_vel(i,1)=0;
    for l=2:length(diff_p2_foot_pos(1,:))
        diff_p2_foot_vel(i,l)=(diff_p2_foot_pos(i,l)-diff_p2_foot_pos(i,l-1))/0.0005;
    end
    diff_p3_foot_vel(i,1)=0;
    for l=2:length(diff_p3_foot_pos(1,:))
        diff_p3_foot_vel(i,l)=(diff_p3_foot_pos(i,l)-diff_p3_foot_pos(i,l-1))/0.0005;
    end
    
      diff_p1_plat_acc(i,1)=0;
    for l=2:length(diff_p1_plat_pos(1,:))
        diff_p1_plat_acc(i,l)=(diff_p1_plat_vel(i,l)-diff_p1_plat_vel(i,l-1))/0.0005;
    end
    diff_p2_plat_acc(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p2_plat_acc(i,l)=(diff_p2_plat_vel(i,l)-diff_p2_plat_vel(i,l-1))/0.0005;
    end
    diff_p3_plat_acc(i,1)=0;
    for l=2:length(diff_p3_plat_pos(1,:))
        diff_p3_plat_acc(i,l)=(diff_p3_plat_vel(i,l)-diff_p3_plat_vel(i,l-1))/0.0005;
    end
    diff_p1_foot_acc(i,1)=0;
    for l=2:length(diff_p1_plat_pos(1,:))
        diff_p1_foot_acc(i,l)=(diff_p1_foot_vel(i,l)-diff_p1_foot_vel(i,l-1))/0.0005;
    end
    diff_p2_foot_acc(i,1)=0;
    for l=2:length(diff_p2_plat_pos(1,:))
        diff_p2_foot_acc(i,l)=(diff_p2_foot_vel(i,l)-diff_p2_foot_vel(i,l-1))/0.0005;
    end
    diff_p3_foot_acc(i,1)=0;
    for l=2:length(diff_p3_plat_pos(1,:))
        diff_p3_foot_acc(i,l)=(diff_p3_foot_vel(i,l)-diff_p3_foot_vel(i,l-1))/0.0005;
    end  
   
    
    diff_p4_plat_pos(i,:)=p4_plat_pos(i,p4_peakst(i)-100:p4_peakst(i)+300);
    diff_p4_plat_torque(i,:)=p4_act_torque(i,p4_peakst(i)-100:p4_peakst(i)+300);
    
    diff_p4_plat_torque(i,:)=diff_p4_plat_torque(i,:)-p0_plat_torque15m;
    
%     diff_p4_foot_pos(i,:)=p4_foot_pos(i,p4_peakst(i)-100:p4_peakst(i)+300)-p0_foot_pos15m;
    diff_p4_foot_pos(i,:)=p4_foot_pos(i,p4_peakst(i)-100:p4_peakst(i)+300);
    diff_p4_foot_pos(i,:)=diff_p4_foot_pos(i,:)-diff_p4_foot_pos(i,100);
    
    
    diff_p4_plat_pos(i,:)=diff_p4_plat_pos(i,:)-diff_p4_plat_pos(i,100);
   
    diff_p4_plat_vel(i,1)=0;
    for l=2:length(diff_p4_plat_pos(1,:))
        diff_p4_plat_vel(i,l)=(diff_p4_plat_pos(i,l)-diff_p4_plat_pos(i,l-1))/0.0005;
    end
    diff_p4_foot_vel(i,1)=0;
    for l=2:length(diff_p4_foot_pos(1,:))
        diff_p4_foot_vel(i,l)=(diff_p4_foot_pos(i,l)-diff_p4_foot_pos(i,l-1))/0.0005;
    end
    
    
      diff_p4_plat_acc(i,1)=0;
    for l=2:length(diff_p4_plat_pos(1,:))
        diff_p4_plat_acc(i,l)=(diff_p4_plat_vel(i,l)-diff_p4_plat_vel(i,l-1))/0.0005;
    end
    
 
    diff_p4_foot_acc(i,1)=0;
    for l=2:length(diff_p4_plat_pos(1,:))
        diff_p4_foot_acc(i,l)=(diff_p4_foot_vel(i,l)-diff_p4_foot_vel(i,l-1))/0.0005;
    end
    
   
    
    
end

diff_p1_plat_torquem=trimmean(diff_p1_plat_torque,30);
diff_p1_plat_posm=trimmean(diff_p1_plat_pos,30);
diff_p1_foot_posm=trimmean(diff_p1_foot_pos,30);
diff_p2_plat_torquem=trimmean(diff_p2_plat_torque,30);
diff_p2_plat_posm=trimmean(diff_p2_plat_pos,30);
diff_p2_foot_posm=trimmean(diff_p2_foot_pos,30);
diff_p3_plat_torquem=trimmean(diff_p3_plat_torque,30);
diff_p3_plat_posm=trimmean(diff_p3_plat_pos,30);
diff_p3_foot_posm=trimmean(diff_p3_foot_pos,30);

diff_p1_plat_accm=trimmean(diff_p1_plat_acc,30);
diff_p1_foot_accm=trimmean(diff_p1_foot_acc,30);
diff_p2_plat_accm=trimmean(diff_p2_plat_acc,30);
diff_p2_foot_accm=trimmean(diff_p2_foot_acc,30);
diff_p3_plat_accm=trimmean(diff_p3_plat_acc,30);
diff_p3_foot_accm=trimmean(diff_p3_foot_acc,30);
diff_p4_plat_torquem=trimmean(diff_p4_plat_torque,30);
diff_p4_plat_posm=trimmean(diff_p4_plat_pos,30);
diff_p4_foot_posm=trimmean(diff_p4_foot_pos,30);

diff_p4_plat_accm=trimmean(diff_p4_plat_acc,30);
diff_p4_foot_accm=trimmean(diff_p4_foot_acc,30);


%% Removing platform dynamics and Impedance estimation

%Estimated platform dynamics are removed and
%individual impances are estimated.
exc_rigid=[p1,p2,p3,p4];
analysis_value=min(exc_rigid);
for i=1:analysis_value-1
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torque(i,:)-0.1945*diff_p1_plat_acc(i,:);%+0.02*diff_p1_foot_accm;%+1*diff_p1_foot_velm;
    diff_p1_plat_torqueimp(i,:)=diff_p1_plat_torqueimp(i,:)-diff_p1_plat_torqueimp(i,100);
    
   
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torque(i,:)-0.1945*diff_p2_plat_acc(i,:);%+0.013*diff_p2_foot_accm;
    diff_p2_plat_torqueimp(i,:)=diff_p2_plat_torqueimp(i,:)-diff_p2_plat_torqueimp(i,100);

    
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torque(i,:)-0.1945*diff_p3_plat_acc(i,:);%-0.09*diff_p3_foot_accm;
    diff_p3_plat_torqueimp(i,:)=diff_p3_plat_torqueimp(i,:)-diff_p3_plat_torqueimp(i,100);

    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torque(i,:)-0.1945*diff_p4_plat_acc(i,:);%+0.02*diff_p4_foot_accm;%+1*diff_p4_foot_velm;
    diff_p4_plat_torqueimp(i,:)=diff_p4_plat_torqueimp(i,:)-diff_p4_plat_torqueimp(i,100);
end





diff_p1_foot_accm=diff_p1_foot_accm;
diff_p2_foot_accm=trimmean(diff_p2_foot_acc,30);
diff_p3_foot_accm=trimmean(diff_p3_foot_acc,30);
diff_p1_plat_torqueimpm=trimmean(diff_p1_plat_torqueimp,30);
diff_p1_plat_torqueimpm=trimmean(diff_p1_plat_torqueimp,30);
diff_p2_plat_torqueimpm=trimmean(diff_p2_plat_torqueimp,30);
diff_p2_plat_torqueimpm=trimmean(diff_p2_plat_torqueimp,30);
diff_p3_plat_torqueimpm=trimmean(diff_p3_plat_torqueimp,30);
diff_p3_plat_torqueimpm=trimmean(diff_p3_plat_torqueimp,30);
diff_p1_foot_posm=trimmean(diff_p1_foot_pos,30);
diff_p2_plat_posm=trimmean(diff_p2_plat_pos,30);
diff_p2_foot_posm=trimmean(diff_p2_foot_pos,30);
diff_p3_plat_posm=trimmean(diff_p3_plat_pos,30);
diff_p3_foot_posm=trimmean(diff_p3_foot_pos,30);
diff_p1_plat_velm=trimmean(diff_p1_plat_vel,30);
diff_p1_foot_velm=trimmean(diff_p1_foot_vel,30);
diff_p2_plat_velm=trimmean(diff_p2_plat_vel,30);
diff_p2_foot_velm=trimmean(diff_p2_foot_vel,30);
diff_p3_plat_velm=trimmean(diff_p3_plat_vel,30);
diff_p3_foot_velm=trimmean(diff_p3_foot_vel,30);
diff_p2_plat_torqueimpm=trimmean(diff_p2_plat_torqueimp,30);
diff_p3_plat_torqueimpm=trimmean(diff_p3_plat_torqueimp,30);

diff_p4_foot_accm=diff_p4_foot_accm;
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);
diff_p4_plat_torqueimpm=trimmean(diff_p4_plat_torqueimp,30);

diff_p4_foot_posm=trimmean(diff_p4_foot_pos,30);
diff_p4_foot_velm=trimmean(diff_p4_foot_vel,30);
diff_p4_plat_velm=trimmean(diff_p4_plat_vel,30);
diff_p4_foot_accm=trimmean(diff_p4_foot_acc,30);
% constained linear regression
%Value is made NaN if foot position is not on
%platform

outlier_criterion_std = 2.5; 

outlier_removal_type = "POS"; % POS | ACC_TOR
%outlier_removal_type = "ACC_TOR";

POS_REJECTION_LIMIT = 2.0;


%%---P1---%%

%Outlier Rejection
if strcmp(outlier_removal_type, "POS")
    [p1_outliers_ind, p1_removed_ind] = rmoutliers(diff_p1_foot_pos(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
elseif strcmp(outlier_removal_type, "ACC_TOR")
    [p1_outliers_ind, p1_removed_ind] = rmoutliers(diff_p1_foot_acc(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
end

p1_removed_ind = p1_removed_ind | (abs(img1_pos') > POS_REJECTION_LIMIT);


%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
diff_p1_pos_data_vec = diff_p1_foot_pos(~p1_removed_ind, 100:300)';
diff_p1_vel_data_vec = diff_p1_foot_vel(~p1_removed_ind, 100:300)';
diff_p1_acc_data_vec = diff_p1_foot_acc(~p1_removed_ind, 100:300)';
p1_data_matrix = [diff_p1_pos_data_vec(:), diff_p1_vel_data_vec(:), diff_p1_acc_data_vec(:)];

p1_torque_data_matrix = diff_p1_plat_torqueimp(~p1_removed_ind, 100:300)';
p1_output_data_vec = p1_torque_data_matrix(:);

A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

p1_regressors = lsqlin(p1_data_matrix, p1_output_data_vec, A1, B1);

%%---P2---%%

%Outlier Rejection
if strcmp(outlier_removal_type, "POS")
    [p2_outliers_ind, p2_removed_ind] = rmoutliers(diff_p2_foot_pos(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
elseif strcmp(outlier_removal_type, "ACC_TOR")
    [p2_outliers_ind, p2_removed_ind] = rmoutliers(diff_p2_foot_acc(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
end

p2_removed_ind = p2_removed_ind | (abs(img2_pos') > POS_REJECTION_LIMIT);

%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
diff_p2_pos_data_vec = diff_p2_foot_pos(~p2_removed_ind, 100:300)';
diff_p2_vel_data_vec = diff_p2_foot_vel(~p2_removed_ind, 100:300)';
diff_p2_acc_data_vec = diff_p2_foot_acc(~p2_removed_ind, 100:300)';
p2_data_matrix = [diff_p2_pos_data_vec(:), diff_p2_vel_data_vec(:), diff_p2_acc_data_vec(:)];

p2_torque_data_matrix = diff_p2_plat_torqueimp(~p2_removed_ind, 100:300)';
p2_output_data_vec = p2_torque_data_matrix(:);

A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

p2_regressors = lsqlin(p2_data_matrix, p2_output_data_vec, A1, B1);

%%---P3---%%

%Outlier Rejection
if strcmp(outlier_removal_type, "POS")
    [p3_outliers_ind, p3_removed_ind] = rmoutliers(diff_p3_foot_pos(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
elseif strcmp(outlier_removal_type, "ACC_TOR")
    [p3_outliers_ind, p3_removed_ind] = rmoutliers(diff_p3_foot_acc(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
end

p3_removed_ind = p3_removed_ind | (abs(img3_pos') > POS_REJECTION_LIMIT);

%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
diff_p3_pos_data_vec = diff_p3_foot_pos(~p3_removed_ind, 100:300)';
diff_p3_vel_data_vec = diff_p3_foot_vel(~p3_removed_ind, 100:300)';
diff_p3_acc_data_vec = diff_p3_foot_acc(~p3_removed_ind, 100:300)';
p3_data_matrix = [diff_p3_pos_data_vec(:), diff_p3_vel_data_vec(:), diff_p3_acc_data_vec(:)];

p3_torque_data_matrix = diff_p3_plat_torqueimp(~p3_removed_ind, 100:300)';
p3_output_data_vec = p3_torque_data_matrix(:);

A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

p3_regressors = lsqlin(p3_data_matrix, p3_output_data_vec, A1, B1);

%%---P4---%%

%Outlier Rejection
if strcmp(outlier_removal_type, "POS")
    [p4_outliers_ind, p4_removed_ind] = rmoutliers(diff_p4_foot_pos(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
elseif strcmp(outlier_removal_type, "ACC_TOR")
    [p4_outliers_ind, p4_removed_ind] = rmoutliers(diff_p4_foot_acc(:, 100:300), 'ThresholdFactor', outlier_criterion_std);
end

p4_removed_ind = p4_removed_ind | (abs(img4_pos') > POS_REJECTION_LIMIT);

%1st Column: Diff Pos.
%2nd Column: Diff Vel.
%3rd Column: Diff Acc.
diff_p4_pos_data_vec = diff_p4_foot_pos(~p4_removed_ind, 100:300)';
diff_p4_vel_data_vec = diff_p4_foot_vel(~p4_removed_ind, 100:300)';
diff_p4_acc_data_vec = diff_p4_foot_acc(~p4_removed_ind, 100:300)';
p4_data_matrix = [diff_p4_pos_data_vec(:), diff_p4_vel_data_vec(:), diff_p4_acc_data_vec(:)];

p4_torque_data_matrix = diff_p4_plat_torqueimp(~p4_removed_ind, 100:300)';
p4_output_data_vec = p4_torque_data_matrix(:);

A1=[-1 0 0;0 -1 0;1 0 0;0 1 0;0 0 -1; 0 0 1];
B1=[0 ;0 ;1000;1000;-1*lim;u_lim];

p4_regressors = lsqlin(p4_data_matrix, p4_output_data_vec, A1, B1);


%% Storing final data



if(plot_figs==1)
    plot_regression_results(diff_p4_foot_pos(~p4_removed_ind, :), diff_p4_foot_vel(~p4_removed_ind, :), ...
                            diff_p4_foot_acc(~p4_removed_ind, :), diff_p4_plat_torqueimp(~p4_removed_ind, :), p4_regressors, '18')
    plot_regression_results(diff_p1_foot_pos(~p1_removed_ind, :), diff_p1_foot_vel(~p1_removed_ind, :), ...
                            diff_p1_foot_acc(~p1_removed_ind, :), diff_p1_plat_torqueimp(~p1_removed_ind, :), p1_regressors, '31')
    plot_regression_results(diff_p2_foot_pos(~p2_removed_ind, :), diff_p2_foot_vel(~p2_removed_ind, :), ...
                            diff_p2_foot_acc(~p2_removed_ind, :), diff_p2_plat_torqueimp(~p2_removed_ind, :), p2_regressors, '44')
    plot_regression_results(diff_p3_foot_pos(~p3_removed_ind, :), diff_p3_foot_vel(~p3_removed_ind, :), ...
                            diff_p3_foot_acc(~p3_removed_ind, :), diff_p3_plat_torqueimp(~p3_removed_ind, :), p3_regressors, '57')

end

%The order is from lower to greater stance phase percentage.
%The odd numbering is due to how the connections are made in the Simulink
%file.
imp_vals_br = [p4_regressors'; p1_regressors'; p2_regressors'; p3_regressors'];


[mean_p4, ci_p4, GoF_p4] = bootstrapping_regression(diff_p4_pos_data_vec', diff_p4_vel_data_vec', diff_p4_acc_data_vec', ...
                         p4_torque_data_matrix', 0.7, A1, B1);
[mean_p1, ci_p1, GoF_p1] = bootstrapping_regression(diff_p1_pos_data_vec', diff_p1_vel_data_vec', diff_p1_acc_data_vec', ...
                         p1_torque_data_matrix', 0.7, A1, B1);
          
[mean_p2, ci_p2, GoF_p2] = bootstrapping_regression(diff_p2_pos_data_vec', diff_p2_vel_data_vec', diff_p2_acc_data_vec', ...
                         p2_torque_data_matrix', 0.7, A1, B1);
[mean_p3, ci_p3, GoF_p3] = bootstrapping_regression(diff_p3_pos_data_vec', diff_p3_vel_data_vec', diff_p3_acc_data_vec', ...
                         p3_torque_data_matrix', 0.7, A1, B1);
                     
imp_vals_bs = [mean_p4; mean_p1; mean_p2; mean_p3];
imp_vals_bs_s = [ci_p4; ci_p1; ci_p2; ci_p3];
     
if(plot_figs==1)
    plot_regression_results(diff_p4_foot_pos(~p4_removed_ind, :), diff_p4_foot_vel(~p4_removed_ind, :), ...
                            diff_p4_foot_acc(~p4_removed_ind, :), diff_p4_plat_torqueimp(~p4_removed_ind, :), mean_p4, '18')
    plot_regression_results(diff_p1_foot_pos(~p1_removed_ind, :), diff_p1_foot_vel(~p1_removed_ind, :), ...
                            diff_p1_foot_acc(~p1_removed_ind, :), diff_p1_plat_torqueimp(~p1_removed_ind, :), mean_p1, '31')
    plot_regression_results(diff_p2_foot_pos(~p2_removed_ind, :), diff_p2_foot_vel(~p2_removed_ind, :), ...
                            diff_p2_foot_acc(~p2_removed_ind, :), diff_p2_plat_torqueimp(~p2_removed_ind, :), mean_p2, '44')
    plot_regression_results(diff_p3_foot_pos(~p3_removed_ind, :), diff_p3_foot_vel(~p3_removed_ind, :), ...
                            diff_p3_foot_acc(~p3_removed_ind, :), diff_p3_plat_torqueimp(~p3_removed_ind, :), mean_p3, '57')

end                     

GoF_matrix = [GoF_p4; GoF_p1; GoF_p2; GoF_p3]


[weight_m, weight_s] = process_weight(weightr15, weightr30, weightr45, weightr60);
[cop_m, cop_s] = process_cop(copr15, copr30, copr45, copr60);
                   