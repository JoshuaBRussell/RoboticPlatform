for i=1:size(weight4)
    weight_15(i)=mean(weightr15(i,50:150));
weight_30(i)=mean(weightr30(i,50:150));
weight_45(i)=mean(weightr45(i,50:150));
weight_60(i)=mean(weightr60(i,50:150));
end
weight_15m=mean(weight_15);
weight_30m=mean(weight_30);
weight_45m=mean(weight_45);
weight_60m=mean(weight_60);
t=tinv([0.025 0.975],(size(weight4,1)))
weight_15s=t(2)*std(weight_15)/sqrt(size(weight4,1));
weight_30s=t(2)*std(weight_30)/sqrt(size(weight4,1));
weight_45s=t(2)*std(weight_45)/sqrt(size(weight4,1));
weight_60s=t(2)*std(weight_60)/sqrt(size(weight4,1));

for i=1:size(cop4)
cop_15(i)=mean(copr15(i,50:150))*100;
cop_30(i)=mean(copr30(i,50:150))*100;
cop_45(i)=mean(copr45(i,50:150))*100;
cop_60(i)=mean(copr60(i,50:150))*100;
end
cop_15m=mean(cop_15);
cop_30m=mean(cop_30);
cop_45m=mean(cop_45);
cop_60m=mean(cop_60);
t=tinv([0.025 0.975],(size(cop4,1)))
cop_15s=t(2)*std(cop_15)/sqrt(size(cop4,1));
cop_30s=t(2)*std(cop_30)/sqrt(size(cop4,1));
cop_45s=t(2)*std(cop_45)/sqrt(size(cop4,1));
cop_60s=t(2)*std(cop_60)/sqrt(size(cop4,1));
weight_rigid=[weight_15m, weight_15s,cop_15m,cop_15s;weight_30m, weight_30s,cop_30m,cop_30s;weight_45m, weight_45s,cop_45m,cop_45s;weight_60m, weight_60s,cop_60m,cop_60s]
% weight_haptic=[weighth_30m, weighth_30s,coph_30m,coph_30s;weighth_45m, weighth_45s,coph_45m,coph_45s;weighth_60m, weighth_60s,coph_60m,coph_60s]

