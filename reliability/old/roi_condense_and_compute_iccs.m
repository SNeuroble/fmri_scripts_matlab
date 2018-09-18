function [icc_results,icc_mean,sum_results] = roi_compute_iccs(roi_data,categories)
%% read inputs
% roi_data and mask assumed to be cell arrays;
% assuming rois order is: [scanner1_subj1_day1, scanner1_subj1_day2, scanner2_subj1_day1, &c]

%addpath('/Users/stephanie/Documents/MATLAB/ICC/'); % use if collapsing to 2


%% Actual (Incomplete) and Complete Data Categorization
% Actual - Incomplete
subj_actual=categories(:,1);
scanner_actual=categories(:,2);
day_actual=categories(:,3);

nsubj=max(subj_actual);
nscanners=max(scanner_actual);
ndays=max(day_actual);

sessions_actual=[subj_actual scanner_actual day_actual];


% Complete - note that these are the values used
subj_complete=repmat(1:nsubj,ndays,nscanners);
scanner_complete=repmat(1:nscanners,nsubj*ndays,1);
day_complete=repmat(1:ndays,1,nsubj*nscanners);

subj=subj_complete(:);
scanner=scanner_complete(:);
day=day_complete(:);

sessions_complete=[subj scanner day];



%% Fill in Missing ROI Data
temp = transpose(cell2mat(roi_data));
meandata=transpose(mean(temp));
clear temp;

for(i=1:(size(sessions_complete,1)))
    missing = ~ismember(sessions_complete(i,:),sessions_actual,'rows');
    
    if (missing)
        warning('Missing data in row:%s. Filling in now.', mat2str(sessions_complete(i,:)));
        roi_data=insertcell(roi_data,meandata,i);
    end
    
end

sessions=(1:size(sessions_complete,1));
edges = [1:size(roi_data{1}(:,1),1)];

clearvars -except roi_data subj scanner day edges sessions nsubj nscanners ndays;

%% setup

r0 = 0.5; % null hypothesis r (fraction variance explained by subj)
alpha = 0.05;
eps=1E-5;

icc_results = [];
sum_results = [];
more_results = 1;

startid = 1; % skip first empty row of roi

%% Determine Variance and ICC
for (i=edges)
    for(j=sessions)
        thisroi_icd(j) = roi_data{j}(i); % vector of length no_sessions
    end
    
   thisroi_icd=transpose(thisroi_icd);
    
%     % collapse days to create two factors
%      thisroi_icd=reshape(thisroi_icd,[ndays,nsubj*nscanners]);
%      thisroi_icd=mean(thisroi_icd);
%      thisroi_icd=reshape(thisroi_icd,[nsubj,nscanners]);
%      [r, LB, UB, F, df1, df2, p] = ICC(thisroi_icd,'1-k',alpha,r0);
%      icc_results = [icc_results; r, LB, UB, F, df1, df2, p];
    


    % get variance via anova
    [pvals, tbl, stats] = anovan(thisroi_icd, {subj scanner day}, ...
        'model', 3, 'random', [1 2 3], 'display', 'off',  'varnames', ...
        {'Subject' 'Scanner' 'Day'});
    
    
%     % get variance via fitting lme -- come back to this later
%      this_tbl=table(thisroi_icd, subj, scanner, day,'VariableNames', ...
%          {'ICD','Subject','Scanner','Day'});
%     
%      this_lme = fitlme(this_tbl, ...
%          'ICD~Subject+(1|Subject)+(1|Scanner)+(1|Day)+(Subject|Scanner)+(Subject|Day)+(Scanner|Day)', ...
%          'FitMethod', 'REML', 'DummyVarCoding', 'effects');
    
    
%     % pseudo-r based on mean squares
%     MSsubj=tbl{2,5};
%     MSscan=tbl{3,5};
%     MSday=tbl{4,5};
%     MSsu_sc=tbl{5,5};
%     MSsu_d=tbl{6,5};
%     MSsu_sc_d=tbl{8,5};
%     
%     r = MSsubj/(MSsubj+MSsu_sc/nscanners+MSsu_d/ndays+MSsu_sc_d/(nscanners*ndays));
    


    % r based on estimated variance
    if(tbl{2,13}>0) varsubj=tbl{2,13}; else varsubj=0; end
    if(tbl{5,13}>0) varsu_sc=tbl{5,13}; else varsu_sc=0; end
    if(tbl{6,13}>0) varsu_d=tbl{6,13}; else varsu_d=0; end
    if(tbl{8,13}>0) varsu_sc_d=tbl{8,13}; else varsu_sc_d=0; end
    
    
    r = varsubj/(varsubj+varsu_sc/nscanners+varsu_d/ndays+varsu_sc_d/(nscanners*ndays));
   
    %if (r < 0) r=0; end % for dealing with negatives

    icc_results = [icc_results; r];
    
    if (more_results == 1)
        sum_results{i} = [tbl(:,1) tbl(:,5:7) tbl(:,13)];
    end

end

bins=100;
hist(icc_results(:,1),bins) %single

icc_mean=mean(icc_results(:,1));

