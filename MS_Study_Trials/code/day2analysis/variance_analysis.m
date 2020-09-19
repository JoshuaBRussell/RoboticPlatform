clear all
load('mshaptic.mat');



for i=1:size(rigid_foot)
foot_var(i,1)=var(rigid_foot(i,:));

foot_var(i,2)=var(comp_foot(i,:));



end

for i=1:size(rigid_plat)
plat_var(i,2)=var(comp_plat(i,:));
plat_var(i,1)=var(rigid_plat(i,:));

end



var_value=[mean(foot_var);mean(plat_var)]