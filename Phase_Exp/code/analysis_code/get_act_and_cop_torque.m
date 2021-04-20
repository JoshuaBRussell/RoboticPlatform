function [act_torque,cop_torque] = get_act_and_cop_torque(foot_pos,force_plate_data)
%Inputs: 
%    - Foot Pos: vector containing foot displacement from foot axis.
%    - Force Plate Data: Struct containing all 4 force sensor data

%Outputs: 
%    - Act Torque 
%    - CoP Torque

for i=1:size(force_plate_data.F1, 1)
    act_torque(i,:)=force_plate_data.F1(i,:)*(0.315-(foot_pos(i)/100))+force_plate_data.F4(i,:)*(0.315-(foot_pos(i)/100))-force_plate_data.F2(i,:)*(0.105+(foot_pos(i)/100))-force_plate_data.F3(i,:)*(0.105+(foot_pos(i)/100))-force_plate_data.F5(i,:)*(0.095)-force_plate_data.F6(i,:)*(0.095);
    cop_torque(i,:)=force_plate_data.F1(i,:)*(0.315-(foot_pos(i)/100))+force_plate_data.F4(i,:)*(0.315-(foot_pos(i)/100))-force_plate_data.F2(i,:)*(0.105+(foot_pos(i)/100))-force_plate_data.F3(i,:)*(0.105+(foot_pos(i)/100))-force_plate_data.F5(i,:)*(0.025)-force_plate_data.F6(i,:)*(0.025);
end

end

