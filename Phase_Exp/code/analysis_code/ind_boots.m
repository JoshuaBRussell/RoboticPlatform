loops=100;
ideal=floor((analysis_value-1)*1);


for i=1:loops
    unique1=0;
    unique2=0;
    unique3=0;
    unique4=0;
    unique5=0;
    unique6=0;
    
    while(unique1<60)
        [bindp1,idx1]=datasample(p1imps(1:p1_fin2,:),p1_fin2);
        numb=length(unique(idx1));
        unique1=100*(numb/p1_fin2);
    end
    bootind1(i,:)=mean(bindp1);
    unique1=0;
    while(unique2<60)
        [bindp2,idx2]=datasample(p2imps(1:p2_fin2,:),p2_fin2);
        numb=length(unique(idx2));
        unique2=100*(numb/p2_fin2);
    end
    bootind2(i,:)=mean(bindp2);
    unique2=0;
    while(unique3<60)
        [bindp3,idx3]=datasample(p3imps(1:p3_fin2,:),p3_fin2);
        numb=length(unique(idx3));
        unique3=100*(numb/p3_fin2);
    end
    unique3=0;
    bootind3(i,:)=mean(bindp3);
    while(unique4<60)
        [bindp4,idx4]=datasample(p4imps(1:p4_fin2,:),p4_fin2);
        numb=length(unique(idx4));
        unique4=100*(numb/p4_fin2);
    end
    bootind4(i,:)=mean(bindp4);
    unique4=0;
    while(unique5<60)
        [bindp5,idx5]=datasample(p5imps(1:p5_fin2,:),p5_fin2);
        numb=length(unique(idx5));
        unique5=100*(numb/p5_fin2);
    end
    bootind5(i,:)=mean(bindp5);
    unique5=0;
    while(unique6<60)
        [bindp6,idx6]=datasample(p6imps(1:p6_fin2,:),p6_fin2);
        numb=length(unique(idx6));
        unique6=100*(numb/p6_fin2);
    end
    bootind6(i,:)=mean(bindp6);
    unique6=0;
    
end

boot_ind(1,:)=nanmean(bootind1);
boot_ind(2,:)=nanmean(bootind2);
boot_ind(3,:)=nanmean(bootind3);
boot_ind(4,:)=nanmean(bootind4);
boot_ind(5,:)=nanmean(bootind5);
boot_ind(6,:)=nanmean(bootind6);

save('bootstrap_ind.mat','boot_ind');
