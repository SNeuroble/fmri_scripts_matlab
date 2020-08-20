function [icc_summary, var_comp_mean, selected_stats] = calc_roi_iccs(roi_data, factor_tbl, this_group)
% computes 1, 2- or 3-factor ICC (detects from nfactors in ftbl) 
% note: 2- and 3-factor ICC are G-Theory ICC (Webb and Shavelson, 2005)
% [icc_summary,var_comp_mean] = calc_roi_iccs(roi_data, factor_tbl, this_group)
% roi_data assumed to be cell array (x sessions) of 1-, 2-, or 3-D matrices

%   TODO, maybe: if type(roi_data)==table, convert table -> matrix


%% Define factors

num_factors=size(factor_tbl,2);
if(num_factors>4)
    error('Can only do < 5 factors');
else
    disp(sprintf('Doing %d-way ANOVA',num_factors))
end

factor_size=max(factor_tbl);

if (size(roi_data,2) < prod(factor_size))
    warning('Missing data: %s sessions.',mat2str(prod(factor_size)-size(roi_data,2)));
end



%% setup for ICC calculation

selected_stats = {zeros(size(roi_data{1},1),size(roi_data{1},2),size(roi_data{1},3))};


%% Do Subgroups if specified
if strcmp(this_group,'all')
    do_subgroups=0;
else
    do_subgroups=1;
    group{1}=[1,2,3,5,8];
    group{2}=[4,6,7];
end

if do_subgroups
    
    if strcmp(this_group,'Siemens')
        this_group=1;
    elseif strcmp(this_group,'GE')
        this_group=2;
    else
        warning('Choose a valid group: all, Siemens, or GE')
    end
    
    group_ids=find(ismember(factor_tbl(:,2),group{this_group}));
    factor_tbl=factor_tbl(group_ids,1:end);
    group_ids=group_ids';
else
    group_ids=[1:length(roi_data)];
    %     group_ids=group_ids';
end


%% Estimate Variance
for i = 1:size(roi_data{1},1) % roi1
    for j = 1:size(roi_data{1},2) % roi2
        for k = 1:size(roi_data{1},3)
            
            % get single vector across sessions
            it=1;
            for this_session = group_ids % sessions
                this_roidata(it) = roi_data{this_session}(i,j,k);
                it=it+1;
            end
            this_roidata=this_roidata';
             
            
            switch num_factors
                case 1
                    % 1-way ANOVA, residual
%                     CLASSICAL ICC
                    % note: if you want to calculate this, ftbl must be
                    % reduced to single column (e.g., only first = subj)
                    
                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1)}, ...
                        'model', 1, 'random', 1, 'display', 'off',  'varnames', ...
                        {'Reference_Factor'});
                    
                    %                 MS_comp=tbl(2:3,5); % variance components
                    %                 MS_comp=cell2mat(MS_comp);
                    %                 % composition:  1: var_su         2: var_e
                    %
                    %                 if(set_neg_var_zero) % set negative variances to zero
                    %                      MS_comp(MS_comp<0)=0;
                    %                 end
                    %
                    %                 icc_G=(MS_comp(1)-MS_comp(2))/(MS_comp(1)+(length(factor_tbl)/factor_size - 1)*MS_comp(2));
                    %                 % (MSb-MSw)/(MSb+(k-1)*MSw)
                    %                 % within subject = residual
                    %                 icc_D = icc_G;


                    
                case 2
                    % 2-way ANOVA, pairwise interaction+residual
                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1) factor_tbl(:,2)}, ...
                        'model', 1, 'random', 1:num_factors, 'display', 'off',  'varnames', ...
                        {'Ref_Factor' 'Factor2'});
                    
                    %                 var_comp=tbl(2:4,13); % variance components
                    %                 var_comp=cell2mat(var_comp);
                    %                 % composition:  1: var_su         2: var_sc       3: var_su_sc_e
                    %
                    %                 if(set_neg_var_zero) % set negative variances to zero
                    %                      var_comp(var_comp<0)=0;
                    %                 end
                    %
                    % %                 MS_comp=tbl(2:4,5);
                    % %                 MS_comp=cell2mat(MS_comp);
                    % %                 ICC=(MS_comp(1)-MS_comp(3))/(MS_comp(1)+(factor_size(2)-1)*MS_comp(3));
                    %
                    %                 icc_G=var_comp(1)/(var_comp(1)+var_comp(3)/factor_size(2));
                    %                 % r_G = var_su / ( var_su + var_su_sc_e/(n_sc) );
                    %                 icc_D=var_comp(1)/(var_comp(1)+var_comp(2)/factor_size(2)+var_comp(3)/factor_size(2));
                    %                 % r_D = var_su / ( var_su + var_sc/n_sc + var_su_sc_e/(n_sc));
                    
                    
                    
                case 3
                    % 3-way ANOVA, pairwise interactions, triple interaction+residual
                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1) factor_tbl(:,2) factor_tbl(:,3)}, ...
                        'model', 2, 'random', 1:num_factors, 'display', 'off',  'varnames', ...
                        {'Ref_Factor' 'Factor2' 'Factor3'});
                    
                    %                 var_comp=tbl(2:8,13); % variance components
                    %                 % composition:  1: var_su         2: var_sc       3: var_d
                    %                 %               4: var_su_sc      5: var_su_d     6: var_sc_d   7: var_su_sc_d_e
                    %
                    %                 var_comp=cell2mat(var_comp);
                    %                 if(set_neg_var_zero) % set negative variances to zero
                    %                      var_comp(var_comp<0)=0;
                    %                 end
                    %
                    %                 roi_iccs{2,1}(i,j)=var_comp(1)/(var_comp(1)+var_comp(4)/factor_size(2)+var_comp(5)/factor_size(3)+var_comp(7)/(factor_size(2)*factor_size(3)));
                    %                 % r_G = var_su / ( var_su + var_su_sc/n_sc + var_su_d/n_d + var_su_sc_d_e/(n_sc*n_d) );
                    %                 roi_iccs{2,3}(i,j)=var_comp(1)/(var_comp(1)+var_comp(2)/factor_size(2)+var_comp(3)/factor_size(3)+var_comp(4)/factor_size(2)+var_comp(5)/factor_size(3)+var_comp(6)/(factor_size(2)*factor_size(3))+var_comp(7)/(factor_size(2)*factor_size(3)));
                    %                 % r_D = var_su / ( var_su + var_sc/n_sc + var_d/n_d + var_su_sc/n_sc + var_su_d/n_d+var_sc_d/(n_sc*n_d) + var_su_sc_d_e/(n_sc*n_d));

                    
                    
                case 4
                    % 4-way ANOVA, 2- and 3-way interactions except run x session, 4-way interaction+residual
                    % can't do session instead of day

                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1) factor_tbl(:,2) factor_tbl(:,3) factor_tbl(:,4)}, ...
                        'model', 3, 'random', 1:num_factors, 'display', 'off',  'varnames', ...
                        {'Ref_Factor' 'Factor2' 'Factor3' 'Factor4'});
                    
                    
                    % p r sc se
                    %model  =[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 0 1; 0 0 1 1 ; 1 1 1 0; 1 0 1 1];
                    %[pvals, tbl, stats, terms] = anovan(this_roidata,{factor_tbl(:,1) factor_tbl(:,2) factor_tbl(:,4) factor_tbl(:,3)},...
                    %     'model',model,'random',1:num_factors,'display','off','varnames',{'Ref_Factor','Factor2','Factor3','Factor4'}); % identical; p r sess sc
                    
                    
                    %                 var_comp=tbl(2:15,13); % variance components
                    %                 % composition:  1: var_su        2:var_r       3: var_sc      4. var_d
                    %                 %               5: var_su_r      6:var_su_sc     7: var_su_d     8: var_r_sc    9: var_sc_d       [NESTED: var_r_d]
                    %                 %               10: var_su_r_sc     11:var_su_sc_d                [NESTED: var_su_r_d]
                    %                 %               12: var_su_r_sc_d
                    
            end
            
            selected_stats{i,j,k} = [tbl(:,1) tbl(:,5:7) tbl(:,13)];
        end
    end
    
    
    
    
    clear tbl pvals
    
    % report percent completion
    if ~exist('amt_fin_old')
        amt_fin_old=0;
    end
    if (mod(round(i*100/size(roi_data{1},1)),5)==0)
        amt_fin=round(i*100/size(roi_data{1},1));
        if(~(amt_fin==amt_fin_old))
            amt_fin_old=amt_fin;
            fprintf('%0.0f%% finished\n',amt_fin);
        end
    end
    
    
end


%% Estimate ICC

Dstudy_range=5;
Dstepsz=0.125;

[icc_summary,var_comp_mean]=stats_to_icc(selected_stats,factor_tbl,Dstudy_range,Dstepsz);


%% Draw simple figure (see roi_compute_corrs for more complex)

% only for icc_G
do_graph=0;

if(do_graph)
    
    figure
    image(roi_iccs{2,1})
    
    % draw subgroup figure
    do_subgroup=0;
    if(do_subgroup)
        subgroup = [4,6,7];
        subgroup_IDs = [];
        for (i=1:size(subgroup,2))
            subgroup_IDs = [subgroup_IDs; find(factor_tbl(:,2)==subgroup(i))];
        end
        figure(2);
        image(roi_iccs{2,1}(1,subgroup_IDs))
        icc_mean_subgroup = mean(roi_iccs{2,1}(1,subgroup_IDs))
    end
    
end


