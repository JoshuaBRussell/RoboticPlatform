

NUM_OF_PERT_POINTS = 4;
NUM_OF_SUBJS = 9;
NUM_OF_DATA_POINTS = 16;

SUBJ_BW = [87.2, 59.38, 64.03, 62.45, 52.57, 66.14, 76.32, 50.93, 71.54];


%% Open Excel File and Gather Data
SUBJ_NAME_COL    = 1;
STIFFNESS_COL    = 2;
DAMPING_COL      = 3;
INERTIA_COL      = 4;
VAF_COL          = 5;
STIFF_CI_COL     = 6;
TA_COL           = 7;
PL_COL           = 8;
SOL_COL          = 9;
GCA_COL          = 10;
WEIGHT_COL       = 11;
WEIGHT_CI_COL    = 12;
COP_COL          = 13;
COP_CI_COL       = 14;
PRE_COUNT_COL  = 15;
VAF_COUNT_COL    = 16;

M = readcell('./results/DIR_Summary.xlsx','Sheet',1,'Range','A2:P54');
M = cellfun(@rmmissing, M, 'UniformOutput', false);

%Locate Blank/Non-Numeric Cells so they can be removed
missing_loc = cellfun(@isempty, M);
str_loc = cellfun(@isstr, M);
loc_to_remove = ((missing_loc | str_loc));

%Isolate and organize the relevant data 

stiffness_data = M(~loc_to_remove(:, STIFFNESS_COL), STIFFNESS_COL);
stiffness_data = reshape(cell2mat(stiffness_data), NUM_OF_PERT_POINTS, NUM_OF_SUBJS);

Data_Key = {'Stiffness', 'Damping', 'Inertia', 'VAF', 'STIFF_CI', ...
            'TA', 'PL', 'SOL', 'GCA',                             ...
            'Weight', 'Weight_CI', 'CoP', 'CoP_CI',               ... 
            'Pre_Count', 'VAFCount'};

Data_Value = {};
for data_type_col = 2:(NUM_OF_DATA_POINTS)
   data_chunk = M(~loc_to_remove(:, data_type_col), data_type_col); 
   data_chunk = reshape(cell2mat(data_chunk), NUM_OF_PERT_POINTS, NUM_OF_SUBJS);
   
   Data_Value{data_type_col-1} = data_chunk;
end
        
Data_Map = containers.Map(Data_Key, Data_Value);

%% Get Correlation Results

%Stiffness Correlation with CoP, Triceps Surae, and Body-Weight
stiffness = Data_Map('Stiffness')
CoP = Data_Map('CoP')
EMG_SOL = Data_Map('SOL')
EMG_GCA = Data_Map('GCA');
EMG_Triceps_Surae = EMG_GCA + EMG_SOL;
BW = Data_Map('Weight');
corr_data_matrix = [stiffness(:), CoP(:), EMG_Triceps_Surae(:), BW(:)];
[R, P] = corrcoef(corr_data_matrix);

%Partial Correlation
[R_partial, P_partial] = partialcorr(corr_data_matrix);

%% Plots
TITLE_FONT_SIZE = 45;
LABEL_FONT_SIZE = 45;
MARKER_SIZE = 5000;
%---Overall Population---%
figure();
scatter(CoP(:), stiffness(:), MARKER_SIZE, 'r.')
title("Stiffness Vs CoP", 'FontSize', TITLE_FONT_SIZE);
ylabel("Stiffness (N*m/rad)", 'FontSize', LABEL_FONT_SIZE);
xlabel("Center of Pressure (cm)", 'FontSize', LABEL_FONT_SIZE);

figure();
scatter(EMG_Triceps_Surae(:), stiffness(:), MARKER_SIZE, 'r.')
title("Stiffness Vs EMG (Triceps Surae)", 'FontSize', TITLE_FONT_SIZE);
ylabel("Stiffness (N*m/rad)", 'FontSize', LABEL_FONT_SIZE);
xlabel("Muscle Activation (%MVC)", 'FontSize', LABEL_FONT_SIZE);

figure();
scatter(BW(:), stiffness(:), MARKER_SIZE, 'r.')
title("Stiffness Vs BodyWeight", 'FontSize', TITLE_FONT_SIZE);
ylabel("Stiffness (N*m/rad)", 'FontSize', LABEL_FONT_SIZE);
xlabel("Body Weight (N)", 'FontSize', LABEL_FONT_SIZE);

%---Individual---%
plot(CoP, stiffness)
plot(EMG_Triceps_Surae, stiffness)
plot(BW, stiffness)

for i = 1:NUM_OF_SUBJS
    scatter3(CoP(:, i), BW(:, i), stiffness(:, i)); hold on;
end

%% R^2 for each individual and CoP
R2_results = [];
for i = 1:NUM_OF_SUBJS
   corr_matrix = corrcoef(stiffness(:, i), CoP(:, i));
   r_sqr_i = (corr_matrix(1,2))^2;
   
   R2_results(i) = r_sqr_i;
end

%% ANOVA Analysis
anova1(stiffness') %Checks stiffness across stance phase percentages
%anova1(stiffness)


%% Create Linear Regression Models for Overall Population

%Absolute Values
A1 = [ones(length(stiffness(:)), 1), CoP(:), EMG_Triceps_Surae(:), BW(:)];
b = stiffness(:);

x1 = A1\b;

y1 = A1*x1;

A2 = [ones(length(stiffness(:)), 1), zscore(CoP(:)), zscore(EMG_Triceps_Surae(:)), zscore(BW(:))];
x2 = A2\b;
y2 = A2*x2;

%Body Weight Normalized
A_norm = [ones(length(stiffness(:)), 1), CoP(:), EMG_Triceps_Surae(:), BW(:)];

stiffness_bw_norm = stiffness./SUBJ_BW;
b_norm = stiffness_bw_norm(:);

x_norm = A_norm\b_norm;

y_norm = A_norm*x_norm;


%Differential 
diff_stiffness = stiffness - stiffness(1, :);
diff_CoP = CoP - CoP(1, :);
diff_EMG_TS = EMG_Triceps_Surae - EMG_Triceps_Surae(1, :);
diff_BW = BW - BW(1, :);

A_diff = [ones(length(diff_stiffness(:)), 1), diff_CoP(:), diff_EMG_TS(:), diff_BW(:)];
b_diff = diff_stiffness(:);
x1_diff = A_diff\b_diff;

y_diff = A_diff*x1_diff;

%Body Weight Normalized Differential
diff_stiffness = stiffness - stiffness(1, :);
diff_stiffness_norm = diff_stiffness./SUBJ_BW;
diff_CoP = CoP - CoP(1, :);
diff_EMG_TS = EMG_Triceps_Surae - EMG_Triceps_Surae(1, :);
diff_BW = BW - BW(1, :);

A_diff_norm = [ones(length(diff_stiffness(:)), 1), diff_CoP(:), diff_EMG_TS(:), diff_BW(:)];
b_diff_norm = diff_stiffness_norm(:);
x_diff_norm = A_diff_norm\b_diff_norm;

y_diff_norm = A_diff_norm*x_diff_norm;

%% Get Accuracy of Models

%VAF Absolute
VAF_Abs_1 = 100*(1 - var(stiffness(:) - y1)/var(stiffness(:)));
VAF_Abs_2 = 100*(1 - var(stiffness(:) - y2)/var(stiffness(:)));

%VAF BW Normalized
VAF_norm = 100*(1 - var(stiffness_bw_norm(:) - y_norm)/var(stiffness_bw_norm(:)));

%VAF Differential
VAF_Diff = 100*(1 - var(diff_stiffness(:) - y_diff)/var(diff_stiffness(:)));

%VAF Differential BW Normalized
VAF_diff_norm = 100*(1 - var(diff_stiffness_norm(:) - y_diff_norm)/var(diff_stiffness_norm(:)));



%% Get Linear Regression Models for Individuals
coeff_matrix = [];
model_out_matrix = [];
for i = 1:NUM_OF_SUBJS
   Ai = [ones(size(stiffness, 1), 1), CoP(:, i)];%, EMG_Triceps_Surae(:, i), BW(:, i)];
   bi = stiffness(:, i);

   xi = Ai\bi;

   yi = Ai*xi;
   coeff_matrix(:, i) = xi;
   model_out_matrix(:, i) = yi;
end

%% Get VAF for each individual subject
VAF_results = [];
for i = 1:NUM_OF_SUBJS
   VAF_i = 100*(1 - var(stiffness(:, i) - model_out_matrix(:, i))/var(stiffness(:, i)));
   VAF_results(i) = VAF_i;
end

%% Get Linear Regression Models for Individuals (BW and CoP)
coeff_matrix = [];
model_out_matrix = [];
for i = 1:NUM_OF_SUBJS
   Ai = [ones(size(stiffness, 1), 1), CoP(:, i), BW(:, i)];%, EMG_Triceps_Surae(:, i), BW(:, i)];
   bi = stiffness(:, i);

   xi = Ai\bi;

   yi = Ai*xi;
   coeff_matrix(:, i) = xi;
   model_out_matrix(:, i) = yi;
end

%% Get VAF for each individual subject
VAF_results = [];
for i = 1:NUM_OF_SUBJS
   VAF_i = 100*(1 - var(stiffness(:, i) - model_out_matrix(:, i))/var(stiffness(:, i)));
   VAF_results(i) = VAF_i;
end


%% Look at residuals

residual_abs = stiffness(:) - y1;
res_abs_chunk = reshape(residual_abs, 4,9);