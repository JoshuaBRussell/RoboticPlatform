function[actualpeaks]=find_all_peaks(emg_data,aa,ab,dir)

k=255;
ct=1;

%Finds the first jump in a series of jumps plus its location.
%If the dir variable is positive, a 1 is associated with it.
%If the dir variable is negative, a 2 is associated with it. 
for j=aa:ab
[a,b]=findpeaks(-1*dir*emg_data{1,j}.data(:,18));
sizep=size(a);

    for i=1:sizep(1,1)
 
        if(and(and(k~=a(i),a(i)~=0),dir>0))
          actualpeaks(ct,1)=j;
          actualpeaks(ct,2)=1;
          actualpeaks(ct,3)=b(i);
          k=a(i);
          ct=ct+1;
        end

        if(and(and(k~=a(i),a(i)~=0),dir<0))
          actualpeaks(ct,1)=j;
          actualpeaks(ct,2)=2;
          actualpeaks(ct,3)=b(i);
          k=a(i);
          ct=ct+1;
        end
    end
end
a=0;
b=0;
end

