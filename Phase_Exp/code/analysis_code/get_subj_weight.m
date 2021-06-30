function [weight_list] = get_subj_weight(data_dir, weight_file, subj_list)
%Finds the subject weight using stored in weight_file in data_dir
%for the specified subjects list in subj_list

weight_list = [];
for subj_index = 1:length(subj_list) 
    
    weight_file_handle = fopen(strcat(data_dir, subj_list{subj_index}, '/', weight_file));
    
    weight_file_raw_data = fread(weight_file_handle);
    
    weight_file_data = SimulinkRealTime.utils.getFileScopeData(weight_file_raw_data);
    
    weight_list(subj_index, :) = nanmean(weight_file_data.data(:, 1));
    
end

end

