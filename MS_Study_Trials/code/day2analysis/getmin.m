function[y]=getmin(peaks,img_st,Img)

k=round(peaks/20)-img_st;
if(k+100<length(Img))
[M,I]=min(Img(k:k+100,4));
y=Img(k+I-1,5);
else
    y=10;
end
    