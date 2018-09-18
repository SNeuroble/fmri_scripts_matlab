function imputeddata=imputeTRTdata(data)

imputeddata=data;

nsess=4;
nsub=12;

for i=1:nsub
allmeans(i)=nanmean(data((i-1)*nsess+1:i*nsess));
meanstbl((i-1)*4+1:i*4)=allmeans(i);
end

imputeddata(isnan(data))=meanstbl(isnan(data));