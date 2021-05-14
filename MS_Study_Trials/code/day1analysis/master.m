% The data analysis can be run by using master.m.
% Change the characters at the top of the code (perturb, subject) based on the trial and subject name.
% Data will be stored in variables dorsi_imp, plantar_imp, inver_imp and ever_imp.
% The colums are arranged as follows : Stiffness,Damping,Goodness,TA,SOL,PL,GCA,CoP

clc;
clear;

%%
% Select the perturbation : D- DOrsiflexion p-Plantarflexion I-Inversion E-Eversion
perturb='E';
% Use 1 to generate figures and 0 to not plot any figures
plotfig=1;
% Use the first letter of the first name for files
Subject='O';

%Select which platform to analyze:
platform_selection = 'R'; % 'L' or 'R'


DATA_FOLDER_REL_LOC = "./../../data/OmikDualTesting/Standing/Right/"

mvc_evaluation;

save('initial');
%% Analysis Begins
% IE_foot_gonio=33.3011;
if perturb=='D'
    dir=1;
    direction=1;
    signal=1;
    pfile=strcat(Subject,'DPP.DAT');
    
    clearvars -except dorsi_imp plantar_imp inver_imp ever_imp dir direction signal pfile eval_trial DATA_FOLDER_REL_LOC platform_selection
    load('initial');
    
    efile=strcat(Subject,'DP');
    etitle='DP standing';
    [DP_foot_gonio, DP_plat_gonio] = get_gonio_sf(DATA_FOLDER_REL_LOC, perturb); %Finding gains of goniometer
    dp_eval;
    dorsi_imp(1,:)=imp;
    
    
    
elseif perturb=='E'
    dir=1;
    direction=3;
    signal=1;
    pfile=strcat(Subject,'IEP');
    clearvars -except dorsi_imp plantar_imp inver_imp ever_imp dir direction signal pfile eval_trial DATA_FOLDER_REL_LOC platform_selection
    load('initial');
    efile=strcat(Subject,'IE');
    etitle='IE standing';
    
    [IE_foot_gonio, IE_plat_gonio] = get_gonio_sf(DATA_FOLDER_REL_LOC, perturb); %Finding gains of goniometer
    e_eval;
    ever_imp(1,:)=imp;
    
    
    
    
end

clearvars -except dorsi_imp plantar_imp inver_imp ever_imp
