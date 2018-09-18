function [s_summary]=acrossruncorr(data,ftbl)
%just a little test
% add other HC sessions
% try across sessions
% which edges underlie within-run similarities? Do test-retest at the within-run level
% compare: decreasing similarity with time. Pool together observations 1 away, 2 away, &c.
% Then have 5 1 away, 4 2 away, 3 3 away, &c.

%% Setup
multisession=1;
docontrol=0;

% make it all ftbl=[subj sess run]
if size(ftbl,2)==5 % TRT
    ftbl=[ftbl(:,1) ftbl(:,5) ftbl(:,4)];
else % HC
    ftbl=[ftbl(:,1) ones(size(ftbl,1),1) ftbl(:,2)];
end

inter=0.85;
scaler=800;


%% Go

if multisession % Multisession
    
    for i=unique(ftbl(:,1))'
        for j=unique(ftbl(:,2))'
            id=find(ftbl(:,1)==i & ftbl(:,2)==j & ftbl(:,3)==1)-1;
            for k=unique(ftbl(:,3))'
                for k2=unique(ftbl(:,3))'
                    simil{i,j}(k,k2)=corr(data{id+k}(:),data{id+k2}(:));
                end
            end
        end
    end
    
    
    
    
    
    
    
    
    
else % Single Session
    
    for i=unique(ftbl(:,1))'
        id=find(ftbl(:,1)==i & ftbl(:,2)==1 & ftbl(:,3)==1)-1;
        for k=unique(ftbl(:,3))'
            for k2=unique(ftbl(:,3))';
                simil{i}(k,k2)=corr(data{id+k}(:),data{id+k2}(:));
            end
        end
    end
    
    
    
end


%% as a control
if docontrol
for i=unique(ftbl_new(:,1))'
    for j=1:4
        id(j)=find(ftbl_new(:,1)==i & ftbl_new(:,5)==j & ftbl_new(:,4)==j); 
    end
    
    for k1=1:4
        for k2=1:4
            simil{i}(k1,k2)=corr(data{id(k1)}(:),data{id(k2)}(:));
        end
    end
end
end

%% PLots


n=numel(simil);
nruns=size(simil{1,1},1);

if multisession
    s=reshape(simil,n,1);
    s=cell2mat(s);
    s=reshape(s',nruns,nruns,n);
    s_summary=mean(s,3);
else
    s=cell2mat(simil);
    s=reshape(s,nruns,nruns,n);
    s_summary=mean(s,3);
end

figure; image((s_summary-inter)*scaler)



