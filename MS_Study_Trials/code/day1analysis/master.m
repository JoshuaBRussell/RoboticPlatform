% The data analysis can be run by using master.m.
% Change the characters at the top of the code (perturb, subject) based on the trial and subject name.
% Data will be stored in variables dorsi_imp, plantar_imp, inver_imp and ever_imp.
% The colums are arranged as follows : Stiffness,Damping,Goodness,TA,SOL,PL,GCA,CoP

mvc_evaluation;
%%
% Select the perturbation : D- DOrsiflexion p-Plantarflexion I-Inversion E-Eversion
perturb='E';
% Use 1 to generate figures and 0 to not plot any figures
plotfig=0;
% Use the first letter of the first name for files
Subject='P';

gonio_values; %Finding gains of goniometer
% DP_plat_gonio=-10.9640;
% IE_plat_gonio=-17.9921;
save('initial');
%% Analysis Begins
% IE_foot_gonio=33.3011;
if perturb=='D'
    dir=1;
    direction=1;
    signal=1;
    pfile=strcat(Subject,'DPP.DAT');
    
    clearvars -except dorsi_imp plantar_imp inver_imp ever_imp dir direction signal pfile eval_trial
    load('initial');
    
    efile=strcat(Subject,'DP');
    etitle='DP standing';
    dp_eval;
    dorsi_imp(1,:)=imp;
    
    
    
elseif perturb=='E'
    dir=1;
    direction=3;
    signal=1;
    pfile=strcat(Subject,'IEP');
    clearvars -except dorsi_imp plantar_imp inver_imp ever_imp dir direction signal pfile eval_trial
    load('initial');
    efile=strcat(Subject,'IE');
    etitle='IE standing';
    e_eval;
    ever_imp(1,:)=imp;
    
    
    
    
end

clearvars -except dorsi_imp plantar_imp inver_imp ever_imp
