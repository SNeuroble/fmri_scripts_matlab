function imgcorr_variablenruns(data,ftbl)
% this is a perm test

nsub=max(ftbl(:,1));
nperms=10;

length_tmp=sum(ftbl(:,1)==1);
halflength_tmp=floor(length_tmp/2);

for perm=[1:nperms]
for i=[1:nsub]
    for j=1:halflength_tmp % "window" size
        data_tmp=data(ftbl(:,1)==i);

        %choose half
        rand_ids=randperm(length(data_tmp));
        data_avghalf=mean(reshape(cell2mat(data_tmp(rand_ids(1:halflength_tmp))),[278,278,halflength_tmp]),3);

        data_other=mean(reshape(cell2mat(data_tmp(rand_ids(halflength_tmp+1:halflength_tmp+j))),[278,278,j]),3);
    
        allcorrs(i,j,perm)=corr(data_avghalf(:),data_other(:));
        
        
    end
%     disp(fprintf('subj %d done',i))
end
disp(fprintf('perm %d done',perm))
end
    

allcorrs_mean=mean(allcorrs,3);
allcorrs_mean=mean(allcorrs_mean,1);
plot(allcorrs_mean')

xlabel('no. runs (6 min)')
ylabel('Correlation to other half of the data (rM)')
    

end
