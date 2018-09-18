function [roi_iccs,roi_stats] = calc_roi_iccs_old(roi_data,categories,reference,do_all_stats)
%% read inputs
% example: [iccs,stats] = calc_roi_iccs(data,categories,'subject',1)
% roi_data (single roi seed connectivity vector) assumed to be cell arrays
%   assuming roi_data order is: [scanner1_subj1_day1, scanner1_subj1_day2, scanner2_subj1_day1, &c]
%   TODO: if table, convert table -> matrix

%addpath('/Users/stephanie/Documents/MATLAB/ICC/'); % use if collapsing to 2


%% Define Categories
fill_in_missing_data=0; % "actual" data may be incomplete

% choose reference; kind of sloppy
switch reference
    case 'subject'
        subj_actual=categories(:,1);
        scanner_actual=categories(:,2);
        day_actual=categories(:,3);
    case 'scanner'
        subj_actual=categories(:,2);
        scanner_actual=categories(:,1);
        day_actual=categories(:,3);
    case 'day'
        subj_actual=categories(:,3);
        scanner_actual=categories(:,2);
        day_actual=categories(:,1);
end



nsubj=max(subj_actual);
nscanners=max(scanner_actual);
ndays=max(day_actual);

sessions_actual=[subj_actual scanner_actual day_actual];

if (size(sessions_actual,1) < (nsubj*nscanners*ndays))
    warning('Missing data: %s sessions.',mat2str((nsubj*nscanners*ndays)-(size(sessions_actual,1))));
end

% IMPORTANT: if true, MUST NOTE these values are hard-coded such that subj=col 1
if(fill_in_missing_data)
    % new definition of names
    subj_complete=repmat(1:nsubj,ndays,nscanners);
    scanner_complete=repmat(1:nscanners,nsubj*ndays,1);
    day_complete=repmat(1:ndays,1,nsubj*nscanners);

    subj=subj_complete(:);
    scanner=scanner_complete(:);
    day=day_complete(:);

    sessions_complete=[subj scanner day];


    %% Fill in Missing ROI Data

    % sometimes now input is mat, but code built for cell
    if(~iscell(roi_data(1,1)))
        roi_data = num2cell(roi_data,1);
    end

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
else
    subj=subj_actual;
    scanner=scanner_actual;
    day=day_actual;
    sessions_complete=sessions_actual;
    clearvars subj_actual scanner_actual day_actual sessions_actual
end

%sessions=(1:size(sessions_complete,1));
%edges = [1:size(roi_data{1}(:,1),1)];

clearvars -except do_all_stats roi_data subj scanner day edges sessions nsubj nscanners ndays categories sessions_complete;

%% setup for ICC calculation

r0 = 0.5; % null hypothesis r (fraction variance explained by subj)
alpha = 0.05;
eps=1E-5;

roi_iccs = {'r_G','r_G_mean','r_D','r_D_mean'};
roi_stats = [];
set_neg_var_zero=1;

startid = 1; % skip first empty row of roi

%% Determine Variance and ICC
for ( i = [1:size(roi_data{1},1)] )

    for ( j = [1:size(roi_data{1},2)] )

        for( k = [1:size(sessions_complete,1)] )
            thisroi_icd(k) = roi_data{k}(i,j); % vector of length no_sessions
        end

        thisroi_icd=thisroi_icd';

        % get variance via 3-way anova
        [pvals, tbl, stats] = anovan(thisroi_icd, {subj scanner day}, ...
            'model', 2, 'random', [1 2 3], 'display', 'off',  'varnames', ...
            {'Subject' 'Scanner' 'Day'});

    %     % get variance via fitting lme -- come back to this later
    %      this_tbl=table(thisroi_icd, subj, scanner, day,'VariableNames', ...
    %          {'ICD','Subject','Scanner','Day'});
    %      this_lme = fitlme(this_tbl, ...
    %          'ICD~Subject+(1|Subject)+(1|Scanner)+(1|Day)+(Subject|Scanner)+(Subject|Day)+(Scanner|Day)', ...
    %          'FitMethod', 'REML', 'DummyVarCoding', 'effects');



        % get variance estimates
        if(set_neg_var_zero) % set negative variances to zero
            if(tbl{2,13}>0) var_su=tbl{2,13}; else var_su=0; end
            if(tbl{5,13}>0) var_su_sc=tbl{5,13}; else var_su_sc=0; end
            if(tbl{6,13}>0) var_su_d=tbl{6,13}; else var_su_d=0; end
            if(tbl{8,13}>0) var_su_sc_d_e=tbl{8,13}; else var_su_sc_d_e=0; end

            if(tbl{3,13}>0) var_sc=tbl{3,13}; else var_sc=0; end
            if(tbl{4,13}>0) var_d=tbl{4,13}; else var_d=0; end
            if(tbl{7,13}>0) var_sc_d=tbl{7,13}; else var_sc_d=0; end
        else
            var_su=tbl{2,13};
            var_su_sc=tbl{5,13};
            var_su_d=tbl{6,13};
            var_su_sc_d_e=tbl{8,13};

            var_sc=tbl{3,13};
            var_d=tbl{4,13};
            var_sc_d=tbl{7,13};
        end

        % get r
        r_G = var_su/(var_su+var_su_sc/nscanners+var_su_d/ndays+var_su_sc_d_e/(nscanners*ndays));
        r_D = var_su/(var_su+var_sc/nscanners+var_d/ndays+var_su_sc/nscanners+var_su_d/ndays+var_sc_d/(nscanners*ndays)+var_su_sc_d_e/(nscanners*ndays));

    %     % pseudo-r based on mean squares
    %     MSsubj=tbl{2,5};
    %     MSscan=tbl{3,5};
    %     MSday=tbl{4,5};
    %     MSsu_sc=tbl{5,5};
    %     MSsu_d=tbl{6,5};
    %     MSsu_sc_d=tbl{8,5};
    %     r = MSsubj/(MSsubj+MSsu_sc/nscanners+MSsu_d/ndays+MSsu_sc_d/(nscanners*ndays));


        roi_iccs{2,1}(i,j) = [r_G];
        roi_iccs{2,3}(i,j) = [r_D];

        if (do_all_stats == 1)
            roi_stats{i,j} = [tbl(:,1) tbl(:,5:7) tbl(:,13)];
        end
    end

    if (mod(i*100/size(roi_data{1},1),1)==0)
        fprintf('%0.0f%% finished\n',i*100/size(roi_data{1},1));
    end
    
    if ~exist('amt_fin_old')
        amt_fin_old=0;
    end
    
    if (mod(round(i*100/size(roi_data{1},1)),5)==0)
        amt_fin_new=round(i*100/size(roi_data{1},1));
        if(~(amt_fin_new==amt_fin_old))
            amt_fin_old=amt_fin_new;
            fprintf('%0.0f%% finished\n',amt_fin_new);
        end
    end

end

icc_G_mean = mean(roi_iccs{2,1});
icc_D_mean = mean(roi_iccs{2,3});
roi_iccs{2,2}=icc_G_mean;
roi_iccs{2,4}=icc_D_mean;


%% Draw simple hist (see roi_compute_corrs for more complex)
% only draws for icc_G
draw_img=0;

if(draw_img)
    %bins = 100;

    figure();
    %hist(roi_iccs{2,1},bins) %single
    image(roi_iccs{2,1});

    report_subgroup=0;
    if(report_subgroup)
        % draw subgroup hists
        subgroup = [4,6,7];
        subgroup_IDs = [];

        for (i=1:size(subgroup,2))
            subgroup_IDs = [subgroup_IDs; find(categories(:,2)==subgroup(i))];
        end

        figure(2);
        image(roi_iccs{2,1}(1,subgroup_IDs))
        %hist(roi_iccs{2,1}(1,subgroup_IDs),bins)
        icc_mean_subgroup = mean(roi_iccs{2,1}(1,subgroup_IDs))
    end
end



