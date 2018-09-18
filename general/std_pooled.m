function pooledstd=std_pooled(data1,data2)
%two vectors of data

mu1=mean(data1);
std1=std(data1);
n1=length(data1);
mu2=mean(data2);
std2=std(data2);
n2=length(data2);

pooledstd=sqrt((n1*std1^2+n2*std2^2)/(n1+n2)+(n1*n2)/(n1+n2)^2*(mu1-mu2)^2);


% a kind of okay appx is just an n-weighted mean of the stds