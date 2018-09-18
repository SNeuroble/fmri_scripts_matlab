function [r,q,sum_masked,data_avg]=simpleimgbehaviorcorr(data,behavior_data,ftbl,type,icc)
% type is either separatesessions or withinsession
% just correlate img with behavior

nsubs=max(ftbl(:,1));
nruns=max(ftbl(:,4));
nsess=max(ftbl(:,5));
ftbl_avg=ftbl(ftbl(:,4)==1,:);



dohalfmask=1;
maskcategory=3;


matdim1=size(data{1},1);
matdim2=size(data{1},2);

trimask=logical(tril(ones(matdim1,matdim2),-1));
% if matdim1==matdim2
%     alreadymasked=1;
% else
%     alreadymasked=1;
% end


if (~exist('maskcategory','var') || isempty(maskcategory)) || (~exist('icc','var') || isempty(icc))
    for i=1:length(data)
        data{i}=data{i}(trimask);
    end
    
    newmatdim=sum(sum(+trimask));
    
else % maskbyreliability
    
    % reliability masks
    mask{1}=icc{2,6}{51}>0.74;
    mask{2}=icc{2,6}{51}>0.6 & icc{2,6}{51}<0.74;
    mask{3}=icc{2,6}{51}>0.4 & icc{2,6}{51}<0.6;
    mask{4}=icc{2,6}{51}<0.4;
    if dohalfmask
        clearvars mask
        tmp=icc{2,6}{51};
        tmp=sort(tmp);
        midthresh=tmp(ceil(length(tmp)/2));
        mask{1}=tmp>=midthresh;
        mask{2}=tmp<=midthresh;
    end
    
    
    for i=1:length(data)
        data{i}=data{i}(mask{maskcategory} & trimask);
    end
    
    newmatdim=sum(sum(+(mask{maskcategory} & trimask)));
end

switch type
    case 'singlesession'
        interval=nsess;
        domeanwithinsub=0;
        
        
    case 'separatesessions'
        interval=1;
        domeanwithinsub=0;
        %         for sub=1:nsubs
        %             behavior_tmp=behavior_data(ftbl_avg(:,1)==sub);
        %             behavior_all(sub)=mean(behavior_tmp);
        %         end
        
    case 'meansession'
        interval=1;
        domeanwithinsub=1;
        
end


if domeanwithinsub
    behavior_all=mean(reshape(behavior_data,[nsess,nsubs]))';
else
    behavior_all=behavior_data(1:interval:end);
end


for thisnruns=1:6
    for sub=1:nsubs
        for sess=1:nsess
            it=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==sess);
            tmp=data(it);
            tmp=mean(reshape(cell2mat(tmp),[newmatdim,length(it)]),2);
            id=(sub-1)*nsess+sess;
            
            data_avg{thisnruns}(:,id)=tmp;
        end
    end
    
    if domeanwithinsub
        tmpdim1=size(data_avg{thisnruns},1);
        data_avg{thisnruns}=mean(reshape(data_avg{thisnruns},[tmpdim1,nsess,nsubs]),2);
        data_avg{thisnruns}=squeeze(data_avg{thisnruns})';
    else
        data_avg{thisnruns}=data_avg{thisnruns}(:,1:interval:end)';
    end
    
    
    
    qthresh=0.05;
    [r{thisnruns},p]=corr(data_avg{thisnruns},behavior_all);
    [p2,q{thisnruns}]=mafdr(p);
    q2_thresh=q{thisnruns}<qthresh;
    
%     if ~alreadymasked
%         for i=1:length(mask)
%             sum_masked(thisnruns,i)=sum(mask{i}.*q2_thresh);
%             perc_masked(thisnruns,i)=sum_masked(thisnruns,i)/sum(mask{i});
%         end
%     else
%         i=0;
%     end
%     
%     sum_masked(thisnruns,i+1)=sum(q2_thresh);
%     perc_masked(thisnruns,i+1)=sum_masked(thisnruns,i+1)/length(q2_thresh);
    
    sum_masked=0; perc_masked=0;

    disp(sprintf('%0.0f runs done', thisnruns))
end

% figure;
% hold on
% for i=1:size(perc_masked,2)
%     plot(1:6,perc_masked(:,i))
% end

