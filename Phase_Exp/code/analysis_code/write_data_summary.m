%% Save Data Summary
base_filename = 'datasummary.xlsx';
filename = strcat(RESULTS_DIR, sub_name, '_', base_filename);
gait_phase_percentages = {num2str(PERT_POINT_1); num2str(PERT_POINT_2);num2str(PERT_POINT_3); num2str(PERT_POINT_4);};
varNames = {'Stiffness'; 'Damping'; 'Inertia'; 'Goodness'; 'Stiffness 95% CI'; ...
            'TA'; 'PL'; 'SOL'; 'GCA'; 'Weight'; 'Weight 95% CI'; 'CoP'; 'CoP 95% CI'; ...
            'Pre Count'};
section_1 = [imp_vals_bs, GoF_matrix, imp_vals_bs_s];
section_2 = [emg_final, weight_m, weight_s, cop_m, cop_s];
section_3 = [, pre_diff_selection_count];

data_table = table(imp_vals_bs(:, 1), imp_vals_bs(:, 2), imp_vals_bs(:, 3), ...
                   GoF_matrix, imp_vals_bs_s(:, 1), ...
                   emg_final(:,1), emg_final(:,2), emg_final(:,3), emg_final(:,4), ...
                   weight_m, weight_s, cop_m, cop_s,  ...
                   pre_diff_selection_count, ...
                   'RowNames', gait_phase_percentages, ...
                   'VariableNames',varNames);
writetable(data_table, filename, 'WriteRowNames',true);

%% Save Plots
hist_fig = figure()
hist([img0_pos, img1_pos, img2_pos, img3_pos, img4_pos]);
saveas(hist_fig, strcat(RESULTS_DIR, sub_name, "_foot_placement.png"));
