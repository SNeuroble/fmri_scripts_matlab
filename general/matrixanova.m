function [MS]=matrixanova(data2,ftbl)

nmatrices=length(data2);
matdim1=size(data2{1},1); matdim2=matdim1;% should be symmetric
for i=1:length(data2)
    data(:,i)=data2{i}(data2{i}==tril(data2{i},-1));
end
% data=reshape(cell2mat(data),[matdim1,matdim2,nmatrices]);


nfactors=size(ftbl,2);
if nfactors > 3
    error('Please specify 3 factors')
end


%% calculate means
X_mean=mean(data,2); % mean of all data
for fctr=1:nfactors
    % main effects
    for i=unique(ftbl(:,fctr))'
        X{1,fctr}(:,i)=mean(data(:,ftbl(:,fctr)==i),2);
        U{1,fctr}(:,i)= X{fctr}(i)-X_mean;
    end

end

% 2- and 3-way interactions
for fctr=1:nfactors
    
    if fctr<nfactors % 2-way
        for fctr2=(fctr+1):nfactors;
            sizef2=max(ftbl(:,fctr2));

            for i=unique(ftbl(:,fctr))'
                for j=unique(ftbl(:,fctr2))'
                    ids=(ftbl(:,fctr)==i) & (ftbl(:,fctr2)==j);
                    X{2,(fctr-2)+fctr2}(:,(i-1)*sizef2+j)=mean(data(:,ids),2);
                    U{2,(fctr-2)+fctr2}(:,(i-1)*sizef2+j)=X{2,(fctr-2)+fctr2}(:,(i-1)*sizef2+j)-(X_mean+U{1,fctr}(:,i)+U{1,fctr2}(:,j)); %% TODO: this subtraction operation doesn't need to be in a loop (use repmat)
                end
            end
        end
        
    else % 3-way
        fsize(1)=max(ftbl(:,1));
        fsize(2)=max(ftbl(:,2));
        fsize(3)=max(ftbl(:,3));
        
        for i=unique(ftbl(:,1))'
            for j=unique(ftbl(:,2))'
                for k=unique(ftbl(:,3))'
                    it=(i-1)*(fsize(2)+fsize(3))+(j-1)*fsize(2)+k;
                    X_expect=X_mean+U{1,1}(:,i)+U{1,2}(:,j)+U{1,3}(:,k)+U{2,1}(:,(i-1)*fsize(2)+j)+U{2,2}(:,(i-1)*fsize(3)+k)+U{2,3}(:,(j-1)*fsize(3)+k); % add interactions
                    U{3,1}(:,it)=data(:,it)-X_expect;
                end
            end
        end
    end
    
end

%% calculate variance - MSD
for i=1:size(U,1)
    for j=1:size(U,2)
        thesefactors=1:3;
        if i==1
            thesefactors(j)=[];
            thesefactors=prod(fsize(thesefactors));
        elseif i==2
            thesefactors=fsize(thesefactors(4-j));
        else
            thesefactors=2;
        end
        
        
        
%         var(i,j)=sum(U{i,j}(:).^2)/(size(U{i,j},2)-1)/size(U{i,j},2);
%         var(i,j)=sum(U{i,j}(:).^2)/(thesefactors*10);
%           MS(i,j)=sum(U{i,j}(:).^2)/(thesefactors-1)/10;
            MS(i,j)=sum(U{i,j}(:).^2)/size(U{i,j},2)/10;
    end
end

% var(1,1)=MS(1,1)-MS(2,1)-MS(2,2);
% var(1,2)=MS(1,2)-MS(2,1)-MS(2,3);
% var(1,3)=MS(1,3)-MS(2,2)-MS(2,3);
% var(2,1)=MS(2,1)-MS(3,1);


% to get R^2: https://en.wikipedia.org/wiki/Coefficient_of_determination#Definitions
