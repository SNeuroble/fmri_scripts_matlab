function glm_results=fit_trav_glm(data,factor_tbl,this_group,ref_category)
% stats computed for one ROI (or edge) at a time
% data is nsessions x nROIs but can be input vice versa, or as cell
% 'this_group': 'all', 'GE' (only w ref_category==site), 'Siemens' (only w ref_category==site)
% 'ref_category': 'subj', 'site', 'day'

% make data usable if cell or wrong orientation 
if(iscell(data))
    data=cell2mat(data);
end

nsessions=prod(max(factor_tbl));
if nsessions ~=length(data)
    nsessions=length(data);
end

if(size(data,2)==nsessions);
    data=data';
end


%% Design Matrix - use specified factor category and group

if strcmp(ref_category,'subj')
    ref_category=1;
elseif strcmp(ref_category,'site')
    ref_category=2;
elseif strcmp(ref_category,'day')
    ref_category=3;
else
    error('Please enter a valid reference category.');
end

if strcmp(this_group,'all')
elseif (strcmp(this_group,'Siemens') || strcmp(this_group,'GE'))
    group_factors{1}=[1,2,3,5,8]; % Siemens
    group_factors{2}=[4,6,7]; % GE
    orig_factor_tbl=factor_tbl; % NEED to keep this
    
    for i=1:size(group_factors,2)
        factor_tbl(ismember(orig_factor_tbl(:,ref_category),group_factors{i}),ref_category)=i;
    end
else
    error('Please enter a valid group.');
end

codingtype=1;
if codingtype==1 % Coding: effects coding
    % create two design matrices; additional one to get last factor
    ref_factor(1)=max(factor_tbl(:,ref_category));
    ref_factor(2)=max(factor_tbl(:,ref_category))-1; % to switch reference
    nlevels=max(factor_tbl(:,ref_category));
        
     for n=1:length(ref_factor)
        dmat{n}=zeros(nsessions,max(factor_tbl(:,ref_category))-1);
        for i=1:nlevels
            dmat{n}((factor_tbl(:,ref_category)==i),i)=1;
        end 
        dmat{n}(factor_tbl(:,ref_category)==ref_factor(n),:)=-1;
        dmat{n}(:,ref_factor(n))=[];
     end

    
elseif codingtype==2 % Coding: factor X vs. mean of others
    ref_factor(1)=max(factor_tbl(:,ref_category));
    ref_factor(2)=max(factor_tbl(:,ref_category))-1; % to switch reference
    nlevels=max(factor_tbl(:,ref_category));
    
    for n=1:length(ref_factor)
        dmat{n}=ones(nsessions,max(factor_tbl(:,ref_category)));
        for i=1:nlevels
            dmat{n}((factor_tbl(:,ref_category)==i),i)=-(nlevels-1);
        end
        dmat{n}(:,ref_factor(n))=[];
    end
    
else
    error('Cannot proceed without coding specified.')
end

clearvars -except dmat data nsessions ref_factor;

%% Fit GLM

for n=1:length(dmat)
    % partial pre-allocation
    beta{n}=zeros(size(dmat{1},2)+1);
    dev_fit{n}=zeros(size(dmat{1},2)+1);
    
    for i=1:size(data,2)
        [beta{n}(:,i),dev_fit{n}(:,i),glm_stats{n}{i}]=glmfit(dmat{n},data(:,i),'normal');
        if i==1
            w = warning('query','last');
            if ~isempty(w)
                id = w.identifier;
                warning('off',id);
            end
        end
    end
end

glm_results={'betas','dev fit','glm stats'; beta, dev_fit, glm_stats};


