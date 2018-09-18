function [overall_corr,data_corr,behav_corr]=comparedata_behav_corrs(data,behavior,ftbl)
% rank similarity of imgs, behavior

nsub=max(ftbl(:,1));
nsess=max(ftbl(:,5));
thisnruns=max(ftbl(:,4));
matdim=size(data{1},1);
datamask=logical(tril(ones(matdim,matdim),-1));
sessmask=logical(tril(ones(nsess,nsess),-1));
submask=logical(tril(ones(nsub,nsub),-1));

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

dobysess=1;


if dobysess
    
    for sub=1:nsub
        it=(sub-1)*nsess;
        for sess=1:nsess
            
            ids=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==sess);
            data_tmp1=data(ids);
            data_avg{(sub-1)*nsess+sess}=mean(reshape(cell2mat(data_tmp1),[matdim,matdim,length(ids)]),3);
            
            data_tmp(:,sess)=data_avg{it+sess}(datamask);
        end
        
        behav_tmp=behavior(:,it+1:it+4);
        
        if dobehdiff % can't do corr, so just take difference
            for thisb1=1:length(behav_tmp)
                for thisb2=1:length(behav_tmp)
                    behav_corr_tmp(thisb1,thisb2)=bmax-abs(behav_tmp(thisb2)-behav_tmp(thisb1));
                end
            end
        else
            behav_corr_tmp=corr(behav_tmp);
        end
        
        data_corr_tmp=corr(data_tmp);
        behav_corr{sub}=behav_corr_tmp(sessmask);
        data_corr{sub}=data_corr_tmp(sessmask);
        
        % linear and rank corr
        [overall_corr{2,1}(sub,1), overall_corr{2,1}(sub,2)]=corr(data_corr{sub},behav_corr{sub});
        [overall_corr{2,2}(sub,1), overall_corr{2,2}(sub,2)]=corr(data_corr{sub},behav_corr{sub},'type','Spearman');
        
    end

    
else
    
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
                behav_corr(thisb1,thisb2)=bmax-abs(behav_tmp(thisb2)-behav_tmp(thisb1));
            end
        end
    else
        behav_corr=corr(behav_tmp);
    end
    
    
    data_corr=corr(data_tmp);
    behav_corr=behav_corr(submask);
    data_corr=data_corr(submask);
    
    % linear and rank corr
    [overall_corr{2,1}(1), overall_corr{2,1}(2)]=corr(data_corr,behav_corr);
    [overall_corr{2,2}(1), overall_corr{2,2}(2)]=corr(data_corr,behav_corr,'type','Spearman');
    
    
end

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
