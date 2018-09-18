function [corrmat,moreinfo]=imgcorr_across_scans(data,ftbl,type,saveimg,dataiscm)
% e.g., [corrmat,avgd]=imgcorr_across_scans(data,ftbl,'session');
% does "fingerprinting" ID and perfect separability

if ~exist('type','var') | isempty(type)
  type = 'session';
end

if ~exist('saveimg','var') | isempty(saveimg)
  saveimg = 0;
end

if ~exist('dataiscm','var') | isempty(dataiscm)
  dataiscm = 0;
  % if this is true, data=corrmat (meaning the correlation of the cmats)
end

nsubs=max(ftbl(:,1));
nruns=max(ftbl(:,4));
nsess=max(ftbl(:,5));
nfactors=[nsubs nruns nsess];

if prod(nfactors) ~= size(ftbl,1)
    error('Amount of data does not match ftbl')
end
    

imgscaling=50;

%% do corr over sessions or runs
if strcmp(type,'session') % avg runs within sessions
   
    matdim=size(data{1},1);
    trilmask=logical(tril(ones(matdim,matdim),-1));

    ftbl_avg=ftbl(ftbl(:,4)==1,:); % single run
    
    % repeat for different n runs within session
    for thisnruns=1:max(ftbl(:,4))
        
        % get avg data over n runs within session
        for sub=1:nsubs
            for session=1:nsess
                it=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==session);
                data_tmp=data(it);
                data_avg{(sub-1)*4+session}=mean(reshape(cell2mat(data_tmp),[matdim,matdim,length(it)]),3);
            end
        end

        % do corr of each session to each other, across all subjects
        for i=1:nsubs*nsess
            for j=1:nsubs*nsess
                corrmat{thisnruns}(i,j)=corr(data_avg{i}(trilmask),data_avg{j}(trilmask));
            end
        end
        
        % do extra stats
        for i=1:nsubs
            for j=1:nsess
                it=(i-1)*nsess+j;
                tmp_corrs=corrmat{thisnruns}(it,:);
                tmp_corrs(it)=[];
                
                withinmask=ftbl_avg(:,1)==i;
                withinmask(it)=[];

                [winnercorr(thisnruns,it),maxid(thisnruns,it)]=max(tmp_corrs);
                correctid(thisnruns,it)=withinmask(maxid(thisnruns,it));
                minwithincorr(thisnruns,it)=min(tmp_corrs(withinmask));
                meanwithincorr(thisnruns,it)=mean(tmp_corrs(withinmask));
                maxbtwcorr(thisnruns,it)=max(tmp_corrs(~withinmask));
                meanbtwcorr(thisnruns,it)=mean(tmp_corrs(~withinmask));
                islarger(thisnruns,it)=minwithincorr(thisnruns,it)>maxbtwcorr(thisnruns,it);
                
            end
        end
        
        % do corr of each session to each other, within subjects
        for i=1:nsess
            for j=1:nsess
                t=corrmat{thisnruns}(i:nsess:end,j:nsess:end);
                corrmat_avg{thisnruns}(i,j)=mean(diag(t)); % could be faster
            end
        end

        figure; image(corrmat{thisnruns}*imgscaling)
        axis square
        if saveimg
            clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
            timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
            savefig(sprintf('corracrossscans_%s.fig',timeID));
            saveas(gcf,sprintf('corracrossscans_%s.png',timeID))
            close(gcf)
        end

    end
    
    moreinfo=[{'MaxID'},{'maxcorr'},{'iswithin'},{'meanwithin'},{'minwithin'},{'meanbetween'},{'maxbtw'};{maxid},{winnercorr},{correctid},{meanwithincorr},{minwithincorr},{meanbtwcorr},{maxbtwcorr}];
    
    
    
    
    
    
elseif strcmp(type,'run') % do across runs
    
    if ~dataiscm
        matdim=size(data{1},1);
        trilmask=logical(tril(ones(matdim,matdim),-1));
        
        for i=1:length(data)
            for j=1:length(data)
                corrmat(i,j)=corr(data{i}(trilmask),data{j}(trilmask));
            end
        end
    else
        corrmat=data;
    end
    
    % ** do between- v within-sub comparison **
    mask_w=zeros(size(corrmat)); % within-sub mask
    for sub=1:nsubs
        id=(sub-1)*nruns;
        mask_w((id+1):(id+nruns),(id+1):(id+nruns))=1;
    end
    mask_w=logical(tril(mask_w,-1));
    mask_b=ones(size(corrmat))-mask_w; % btw-sub mask
    mask_b=logical(tril(mask_b,-1));
    
    % get mean vals
    edges_w=corrmat(mask_w);
    edges_b=corrmat(mask_b);
    edges_tmp_b=edges_b(randperm(length(edges_b),length(edges_w)));
    

    moreinfo=[{'h'},{'p'},{'ci'},{'stats'},{'meaninfo'}];
    [moreinfo{2,1},moreinfo{2,2},moreinfo{2,3},moreinfo{2,4}]=ttest2(edges_w,edges_tmp_b);
    meaninfo=[{'mean within'},{'std within'},{'mean between'},{'std between'}];
    meaninfo{2,1}=mean(edges_w);
    meaninfo{2,2}=std(edges_w);
    meaninfo{2,3}=mean(edges_b);
    meaninfo{2,4}=std(edges_b);
    moreinfo{2,5}=meaninfo;
    
    % do figs
    clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
    timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
    
    % hist
    nbins=40;
    figure
    histogram(edges_w,nbins);
    hold on
    histogram(edges_tmp_b,nbins);
    hold off
    if saveimg
        savefig(sprintf('hist_within_between_%s.fig',timeID));
        saveas(gcf,sprintf('hist_within_between_%s.png',timeID))
        close(gcf)
    end
    
    % cmat
    figure; image(corrmat*50)
    axis square
    if saveimg
        savefig(sprintf('corracrossruns_%s.fig',timeID));
        saveas(gcf,sprintf('corracrossruns_%s.png',timeID))
        close(gcf)
    end
    
    
end
