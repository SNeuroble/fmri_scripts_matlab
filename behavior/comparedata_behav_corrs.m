function [overall_corr,data_corr,behav_corr]=comparedata_behav_corrs(data,behavior,ftbl)
% compares matrix similarities (corr), behavior similarities (corr or diff),
% and then compares the two
% can be done across sub or within sub

dobysess=0;

nsub=max(ftbl(:,1));
nsess=max(ftbl(:,5));
thisnruns=max(ftbl(:,4));
matdim=size(data{1},1);
datamask=logical(tril(ones(matdim,matdim),-1));
% the sess mask ONLY includes within-subject similarities
sessmask=logical(tril(ones(nsub*nsess,nsub*nsess),-1));
submask=logical(tril(ones(nsub,nsub),-1));

dobehdiff=0;


if isa(behavior,'table')
    behavior=table2array(behavior);
end

if size(behavior,2) ~= nsub*nsess
    behavior=behavior';
end

if size(behavior,1) == 1
    dobehdiff=1;
    bmax=max(abs(behavior));
end

sessmask(:,:)=0;

if dobysess
    
    for sub=1:nsub
        it=(sub-1)*nsess;
        sessmask((sub-1)*nsess+1:sub*nsess,(sub-1)*nsess+1:sub*nsess)=1;
        for sess=1:nsess
            
            ids=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==sess);
            data_tmp1=data(ids);
            data_avg{it+sess}=mean(reshape(cell2mat(data_tmp1),[matdim,matdim,length(ids)]),3);
            data_tmp(:,it+sess)=data_avg{it+sess}(datamask);
            

        end
        % mean center behavior by sub
        behav_mean(:,sub)=mean(behavior(:,it+1:it+nsess),2);
        for thisb=1:size(behavior,1)
            behav_tmp(thisb,it+1:it+nsess)=behavior(thisb,it+1:it+nsess)-behav_mean(thisb,sub);
        end
        
    end
   
    
    for sub=1:nsub
        it=(sub-1)*nsess;
        for sess=1:nsess
            
            if dobehdiff % can't do corr, so just take difference
                for thisb1=1:length(behavior)
                    for thisb2=1:length(behavior)
                        behav_corr(thisb1,thisb2)=abs(behavior(thisb2)-behavior(thisb1));
                    end
                end
                
                bmax=max(max(behav_corr));
                behav_corr=bmax-behav_corr;
                
            else
                behav_corr=corr(behavior);
            end
            
        end
    end
        
mask=logical(tril(sessmask,-1));
    
else % do by sub
    
    for sub=1:nsub
        ids=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns);
        data_tmp1=data(ids);
        data_avg{sub}=mean(reshape(cell2mat(data_tmp1),[matdim,matdim,length(ids)]),3);
        
        data_tmp(:,sub)=data_avg{sub}(datamask);
        
        behav_tmp(:,sub)=mean(behavior(:,(sub-1)*nsess+1:sub*nsess));
    end
    
    
    if dobehdiff
        for thisb1=1:length(behav_tmp)
            for thisb2=1:length(behav_tmp)
                behav_corr(thisb1,thisb2)=abs(behav_tmp(thisb2)-behav_tmp(thisb1));
            end
        end
        bmax=max(max(behav_corr));
        behav_corr=bmax-behav_corr;
    else
        behav_corr=corr(behav_tmp);
    end
    
    mask=submask;
end

data_corr=corr(data_tmp);

figure; image((data_corr.*+mask-0.45)*30)
figure; image((behav_corr.*+mask+1))

behav_corr=behav_corr(mask);
data_corr=data_corr(mask);

% linear and rank corr
[overall_corr{2,1}(1), overall_corr{2,1}(2)]=corr(data_corr,behav_corr);
[overall_corr{2,2}(1), overall_corr{2,2}(2)]=corr(data_corr,behav_corr,'type','Spearman');
    

overall_corr{1,1}='Linear Corr';
overall_corr{1,2}='Spearman Corr';

% fprintf('Linear corr\n')
% fprintf('%1.4f ',overall_corr{1})
% fprintf('Spearman corr\n')
% fprintf('%1.4f ',overall_corr{2})
% fprintf('\n')

% controls
% compare with random matching
% compare with matching with [1 2 3 4] (perfect corr w session)
%   also just see how often [1 2 3 4] is the ranking
