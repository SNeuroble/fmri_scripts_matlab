function glm_results=fit_trav_glm_old(data);

% betas computed for one ROI (or edge) at a time
% data is nsessions x nROIs but can be input vice versa, or as cell

% TODO: pass in factor_tbl instead of hard-coding
nsubj=8;
nscanners=8;
ndays=2;

nsessions=nsubj*nscanners*ndays;

% make data usable if cell or if in wrong orientation 
if(iscell(data))
    data=cell2mat(data);
end

if(size(data,2)==nsessions)
    data=data';
end



%hard-coded design matrix
do_scanner_groups=0; % 0 means separate scanners

nsess_sc=nsubj*ndays; % num sessions per scanner

if(~do_scanner_groups)
    dmat=zeros(7,nsessions);
    
    for(i=1:1:nscanners)
        if(~(i==nscanners))
            dmat(i,((i-1)*nsess_sc+1):((i-1)*nsess_sc+nsess_sc))=1;
        else
            dmat(:,((i-1)*nsess_sc+1):((i-1)*nsess_sc+nsess_sc))=-1;
        end
    end
    
end

if(do_scanner_groups)
    
    dmat=zeros(1,nsessions);
    group='Siemens';
    
    if strcmp(group,'Siemens')
        group_factors=[1,2,3,5,8];
    elseif strcmp(group,'GE')
        group_factors=[4,6,7]; % FAKE; ONLY for when need to get other factor's beta
    end
        
    for(i=1:8)
        if(any(i==group_factors))
            dmat(1,((i-1)*nsess_sc+1):((i-1)*nsess_sc+nsess_sc))=1;
        else
            dmat(1,((i-1)*nsess_sc+1):((i-1)*nsess_sc+nsess_sc))=-1;
        end
    end
    
end

dmat=dmat';

clearvars -except dmat data;

% incomplete pre-allocation
beta=zeros(size(data,2));
dev_fit=zeros(size(data,2));
glm_stats{size(data,2)}=0;

for i=1:size(data,2)
    [beta(:,i),dev_fit(:,i),glm_stats{i}]=glmfit(dmat,data(:,i),'normal');
end

glm_results={'betas','dev fit','glm stats'; beta, dev_fit, glm_stats};


