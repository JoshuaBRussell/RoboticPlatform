function [R_vec,x_vec, y_vec, TRANSFORMED_COORDS] = get_effective_shape(vertical_ankle_coord, ...
                                                    ankle_angle_profiles, ...
                                                    cop_profiles)

p_AF_x = zeros(size(cop_profiles));
p_AF_y = zeros(size(cop_profiles));

%Assumes each row is a seperate trial
for trial_index = 1:size(ankle_angle_profiles, 1)
   for time_index = 1:size(ankle_angle_profiles, 2)
       theta = ankle_angle_profiles(trial_index, time_index);
       R = [cos(theta), -sin(theta);
            sin(theta),  cos(theta)];
   
       p_F = [cop_profiles(trial_index, time_index); - vertical_ankle_coord]; 
       p_AF = R*p_F;
       
       p_AF_x(trial_index, time_index) = p_AF(1);
       p_AF_y(trial_index, time_index) = p_AF(2);
   
   end
end

TRANSFORMED_COORDS.p_AF_x = p_AF_x;
TRANSFORMED_COORDS.p_AF_y = p_AF_y;

%Find the radius of curvature for each trial
for trial_index = 1:size(ankle_angle_profiles, 1)
    [xc,yc,R,a] = circfit(p_AF_x(trial_index, 1:size(cop_profiles, 2)), p_AF_y(trial_index, 1:size(cop_profiles, 2)));
    R_vec(trial_index, :) = R;
    
    x_vec(trial_index, :) = xc;
    y_vec(trial_index, :) = yc;
end


end

