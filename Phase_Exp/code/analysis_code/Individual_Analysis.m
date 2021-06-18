
clc;
clear;

%Exploratory script meant to peruse results of individual subject results.



%% Get Total Population Data Trend from Excel Sheet


POP_MEAN = readcell(strcat('./results/DIR_Summary_PrePert_TrqPos'),'Sheet',1,'Range','A66:P70');

pop_data_key = POP_MEAN(1, 2:end);

pop_data_value = {};        
for data_col = 2:length(pop_data_key)+1
    data_chunk = cell2mat(POP_MEAN(2:end, data_col));
    pop_data_value{data_col - 1} = data_chunk;
end

pop_data_map = containers.Map(pop_data_key, pop_data_value);





%% ---- Look at Summary Results ---- %%
SUBJ_NAME = 'Vu';
NUM_OF_DATA_POINTS = 16;
RESULTS_DIR = strcat('./results/', SUBJ_NAME, '/');


%Open Results Excel Sheet
M = readcell(strcat('./results/', SUBJ_NAME,'/', SUBJ_NAME, '_datasummary.xlsx'),'Sheet',1,'Range','A1:P5');
M = cellfun(@rmmissing, M, 'UniformOutput', false);

% Data_Key = {'Stiffness', 'Damping', 'Inertia', 'VAF', 'STIFF_CI', ...
%             'TA', 'PL', 'SOL', 'GCA',                             ...
%             'Weight', 'Weight_CI', 'CoP', 'CoP_CI',               ... 
%             'Pre_Count', 'VAFCount'};

data_key = M(1, 2:end);

data_value = {};        
for data_col = 2:length(data_key)+1
    data_chunk = cell2mat(M(2:end, data_col));
    data_value{data_col - 1} = data_chunk;
end

data_map = containers.Map(data_key, data_value);


%% Quick Summary Plots
% Stiffness vs. CoP
figure();
plot(data_map('CoP'), data_map('Stiffness')); hold on;
plot(pop_data_map('CoP'), pop_data_map('Stiffness'));
ylabel("Stiffness");
xlabel("CoP (m)");
title("Stiffness vs. CoP");
legend([SUBJ_NAME; "Population Avg."]);
saveas(gcf,strcat(RESULTS_DIR, 'Stiffness_vs_CoP.png'));


%Stiffness vs. EMG
figure();
plot(data_map('SOL'), data_map('Stiffness')); hold on;
plot(data_map('GCA'), data_map('Stiffness')); hold on;
plot(data_map('SOL') + data_map('GCA'), data_map('Stiffness')); hold on;
legend(["SOL"; "GCA"; "Triceps Surae"]);
ylabel("Stiffness");
xlabel("% MVC");
title("Stiffness vs. EMG");
saveas(gcf,strcat(RESULTS_DIR, 'Stiffness_vs_MVC.png'));

%Plot Stiffness vs EMG for POP.
figure();
plot(pop_data_map('SOL'),  pop_data_map('Stiffness')); hold on;
plot(pop_data_map('GCA'),  pop_data_map('Stiffness')); hold on;
plot(pop_data_map('SOL') + pop_data_map('GCA'), pop_data_map('Stiffness')); hold on;
legend(["SOL"; "GCA"; "Triceps Surae"]);
ylabel("Stiffness");
xlabel("% MVC");
title("Population: Stiffness vs. EMG");
saveas(gcf,strcat(RESULTS_DIR, 'Population_Stiffness_vs_MVC.png'));

%Stiffness vs. BW
figure();
plot(data_map('Weight'), data_map('Stiffness')); hold on;
plot(pop_data_map('Weight'), pop_data_map('Stiffness'));
ylabel("Stiffness");
xlabel("BW (N)");
title("Stiffness vs. BW");
legend([SUBJ_NAME; "Population Avg."])
saveas(gcf,strcat(RESULTS_DIR, 'Stiffness_vs_Weight.png'));

%% Plot Stiffness Bar Graph

CI_95 = data_map('Stiffness95%CI');
POP_CI_95 = pop_data_map('Stiffness 95% CI'); %Even though they look exactly the same in the Excel file, the the POP and Inidividual subject have a slightly different key string.

y = [data_map('Stiffness'), pop_data_map('Stiffness')];

figure();
bar_graph = bar(y); hold on;

[ngroups, nbars] = size(y);

%Get X coords of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
   x(i, :) = bar_graph(i).XEndPoints; 
end

errorbar(x', y, [CI_95, POP_CI_95], 'k','linestyle','none');

legend([SUBJ_NAME; "Population Avg."], 'Location', 'northwest');
xlabel("Stance Phase Point");
ylabel("Stiffness");
title("Stiffness Over Stance Phase (95% CI)");


saveas(gcf,strcat(RESULTS_DIR, 'Stiffness_BarGraph.png'));

%% Plot CoP Bar Graph
figure();
y = [data_map('CoP'), pop_data_map('CoP')];
bar(y);
legend([SUBJ_NAME; "Population Avg."], 'Location', 'northwest');

xlabel("Stance Phase Point");
ylabel("CoP");
title("CoP Over Stance Phase");

%saveas(gcf,strcat(RESULTS_DIR, 'CoP_BarGraph.png'));

%% Plot EMG Bar Graph
figure();
y = [data_map('SOL') + data_map('GCA'), pop_data_map('SOL') + pop_data_map('GCA')];
bar(y);
legend([SUBJ_NAME; "Population Avg."], 'Location', 'northwest');

xlabel("Stance Phase Point");
ylabel("% MVC (SOL + GCA)");
title("Triceps Suraue Over Stance Phase");

saveas(gcf,strcat(RESULTS_DIR, 'Triceps_Surae_BarGraph.png'));

%% ---- Look at PreSummary Results ---- %%
RESULTS_DIR = './results/';

load(strcat(RESULTS_DIR, "/", SUBJ_NAME, "/", SUBJ_NAME, "_bootstrap_vars.mat"));

%% Detailed Bootstrapp Result Plots

%Stiffness vs. CoP
figure();
scatter(sqrt(bio_factors_p1.CoP)', regress_coeffs_p1(:, 1)); hold on;
scatter(sqrt(bio_factors_p2.CoP)', regress_coeffs_p2(:, 1));
scatter(sqrt(bio_factors_p3.CoP)', regress_coeffs_p3(:, 1));
scatter(sqrt(bio_factors_p4.CoP)', regress_coeffs_p4(:, 1)); hold off;
ylabel("Stiffness");
xlabel("CoP (m)");
title("Stiffness vs. sqrt(CoP) ");

%Stiffness vs. EMG
figure();
p1_EMG = bio_factors_p1.EMG.SOL + bio_factors_p1.EMG.GCA;
p2_EMG = bio_factors_p2.EMG.SOL + bio_factors_p1.EMG.GCA;
p3_EMG = bio_factors_p3.EMG.SOL + bio_factors_p3.EMG.GCA;
p4_EMG = bio_factors_p4.EMG.SOL + bio_factors_p4.EMG.GCA;

scatter(p1_EMG, regress_coeffs_p1(:, 1)); hold on;
scatter(p2_EMG, regress_coeffs_p2(:, 1)); hold on;
scatter(p3_EMG, regress_coeffs_p3(:, 1)); hold on;
scatter(p4_EMG, regress_coeffs_p4(:, 1)); hold on;
legend(["Triceps Surae"]);
ylabel("Stiffness");
xlabel("% MVC");
title("Stiffness vs. EMG");

%Stiffness vs. BW
scatter(bio_factors_p1.Weight, regress_coeffs_p1(:, 1)); hold on;
scatter(bio_factors_p1.Weight, regress_coeffs_p2(:, 1)); hold on;
scatter(bio_factors_p1.Weight, regress_coeffs_p3(:, 1)); hold on;
scatter(bio_factors_p1.Weight, regress_coeffs_p4(:, 1)); hold on;
legend(["Weight"]);
ylabel("Stiffness");
xlabel("Body Weight (N)");
title("Stiffness vs. BW");


%% Perform Preliminary Stats Analysis on the Bootstrap Data
stiffness = [regress_coeffs_p4(:, 1), regress_coeffs_p1(:, 1), regress_coeffs_p2(:, 1), regress_coeffs_p3(:, 1)]
CoP = [bio_factors_p4.CoP, bio_factors_p1.CoP, bio_factors_p2.CoP, bio_factors_p3.CoP ]
BW  = [bio_factors_p4.Weight, bio_factors_p1.Weight, bio_factors_p2.Weight, bio_factors_p3.Weight];

p1_EMG = bio_factors_p1.EMG.SOL + bio_factors_p1.EMG.GCA;
p2_EMG = bio_factors_p2.EMG.SOL + bio_factors_p1.EMG.GCA;
p3_EMG = bio_factors_p3.EMG.SOL + bio_factors_p3.EMG.GCA;
p4_EMG = bio_factors_p4.EMG.SOL + bio_factors_p4.EMG.GCA;
EMG_Triceps_Surae = [p4_EMG, p1_EMG p2_EMG, p3_EMG];
corr_data_matrix = [stiffness(:), CoP(:), EMG_Triceps_Surae(:), BW(:)];
[R, P] = corrcoef(corr_data_matrix);

%Partial Correlation
[R_partial, P_partial] = partialcorr(corr_data_matrix);






%% Try Making a Regression Model
% mdl = fitlm(X, Y);
% anova(mdl);

%%Basic Model
X = [ones(length(stiffness(:)), 1), corr_data_matrix(:, 2:end)];
Y = stiffness(:);

part_coeff_vec = inv(X'*X)*X'*Y;

scatter(X(:, 2), X*part_coeff_vec - Y);
scatter(X(:, 3), X*part_coeff_vec - Y);
scatter(X(:, 4), X*part_coeff_vec - Y);

R_sqrd = (1 - var(X*part_coeff_vec - Y)/var(Y));

%% Sqrt of CoP
X = [ones(length(stiffness(:)), 1), corr_data_matrix(:, 2:end)];
X(:, 2) = sqrt(X(:, 2));
Y = stiffness(:);

part_coeff_vec = inv(X'*X)*X'*Y;

scatter(X(:, 2), X*part_coeff_vec - Y);
scatter(X(:, 3), X*part_coeff_vec - Y);
scatter(X(:, 4), X*part_coeff_vec - Y);

R_sqrd = (1 - var(X*part_coeff_vec - Y)/var(Y));




