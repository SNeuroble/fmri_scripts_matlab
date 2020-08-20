function [icc_summary, var_comp_mean, selected_stats] = calc_roi_iccs(roi_data, factor_tbl, this_group, modeltype)
% computes 1, 2- or 3-factor ICC (detects from nfactors in ftbl) 
% modeltype: 'full' (all interactions) or 'linear'
% note: 2- and 3-factor ICC are G-Theory ICC (Webb and Shavelson, 2005)
% [icc_summary,var_comp_mean] = calc_roi_iccs(roi_data, factor_tbl, this_group,'full')
% roi_data assumed to be cell array (x sessions) of 1-, 2-, or 3-D matrices

%   TODO, maybe: if type(roi_data)==table, convert table -> matrix


%% Define factors
if ~exist('modeltype','var') | isempty(modeltype); modeltype='full'; end

num_factors=size(factor_tbl,2);

if(num_factors>4)
    error('Can only do < 5 factors');
else
    disp(sprintf('Doing %d-way ANOVA, %s model \n(full=all interactions, linear=no interactions)',num_factors,modeltype))
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
    %% DONT DO THE STUFF IN THIS LOOP. Just do:
    % ftbl2=ftbl;
    % ftbl2(ismember(ftbl(:,2),[1,2,3,5,8]),2)=1;
    % ftbl2(ismember(ftbl(:,2),[4,6,7]),2)=2;
    % [icc,var,stats,sigmask]=run_reliability('none',data,ftbl2);
    
    
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
        for k = 1:size(roi_data{1},3) % if it's 3D
            
            % get single vector across sessions
            it=1;
            for this_session = group_ids % sessions
                this_roidata(it) = roi_data{this_session}(i,j,k);
                it=it+1;
            end
            this_roidata=this_roidata';
             
            
            switch num_factors
                case 1 % 1-way ANOVA - CLASSICAL ICC
                    % note: if you want to use this, ftbl must be reduced to single column (e.g., only first = subj)
                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1)}, ...
                        'model', modeltype, 'random', 1, 'display', 'off',  'varnames', ...
                        {'Reference_Factor'});
                    
                    %                 % composition:  1: var_su         2: var_e
                    %                 % (MSb-MSw)/(MSb+(k-1)*MSw)
                    %                 % within subject = residual


                    
                case 2 % 2-way ANOVA
                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1) factor_tbl(:,2)}, ...
                        'model', modeltype, 'random', 1:num_factors, 'display', 'off',  'varnames', ...
                        {'Ref_Factor' 'Factor2'});
                    
                    %                 % composition:  1: var_su         2: var_sc       3: var_su_sc_e
                    %                 % r_G = var_su / ( var_su + var_su_sc_e/(n_sc) );
                    %                 % r_D = var_su / ( var_su + var_sc/n_sc + var_su_sc_e/(n_sc));
                    
                    
                    
                case 3 % 3-way ANOVA - THIS IS THE DEFAULT FOR US
                     [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1) factor_tbl(:,2) factor_tbl(:,3)}, ...
                         'model', modeltype, 'random', 1:num_factors, 'display', 'off',  'varnames', ...
                         {'Ref_Factor' 'Factor2' 'Factor3'});
                    

                    %                 % composition, if suxscxd:  1: var_su         2: var_sc       3: var_d
                    %                 % r_G = var_su / ( var_su + var_su_sc/n_sc + var_su_d/n_d + var_su_sc_d_e/(n_sc*n_d) );
                    %                 % r_D = var_su / ( var_su + var_sc/n_sc + var_d/n_d + var_su_sc/n_sc + var_su_d/n_d+var_sc_d/(n_sc*n_d) + var_su_sc_d_e/(n_sc*n_d));

                    
                    % NEW - needs to be cleaned
                    % repeated measures model
                    userm=0;
                    if userm
                    % if using, move all this above
                    Sub=factor_tbl(factor_tbl(:,3)==1,1);
                    Sess=factor_tbl(factor_tbl(:,3)==1,2);
                    Run=unique(factor_tbl(:,3))-1; % MUST start at 0
                    nrun=length(Run);
                    nscan=size(factor_tbl,1)/nrun;

                    Sub=num2str(Sub);
                    Sess=num2str(Sess);
                    % end the stuff we move above
                    
                    datarm=reshape(this_roidata,nrun,nscan)';
                    
                    t = table(Sub,Sess,datarm(:,1),datarm(:,2),datarm(:,3),datarm(:,4),datarm(:,5),datarm(:,6),...
                    'VariableNames',{'Sub','Sess','data1','data2','data3','data4','data5','data6',});
                
                    
                    rm = fitrm(t,'data1-data6 ~ Sub + Sess','WithinDesign',Run);
                    
                    ranovatbl=ranova(rm);
                    end

                
                
                    
                    
                    
                case 4 % 4-way ANOVA
                    [pvals, tbl, anova_stats] = anovan(this_roidata, {factor_tbl(:,1) factor_tbl(:,2) factor_tbl(:,3) factor_tbl(:,4)}, ...
                        'model', modeltype, 'random', 1:num_factors, 'display', 'off',  'varnames', ...
                        {'Ref_Factor' 'Factor2' 'Factor3' 'Factor4'});
                    
                    % % Alternate model:
                    % % p r sc se
                    %model  =[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 0 1; 0 0 1 1 ; 1 1 1 0; 1 0 1 1];
                    %[pvals, tbl, stats, terms] = anovan(this_roidata,{factor_tbl(:,1) factor_tbl(:,2) factor_tbl(:,4) factor_tbl(:,3)},...
                    %     'model',model,'random',1:num_factors,'display','off','varnames',{'Ref_Factor','Factor2','Factor3','Factor4'}); % identical; p sc r d
                    
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

% Dstudy_range=5;
Dstudy_range=6;
Dstepsz=0.125;

[icc_summary,var_comp_mean]=stats_to_icc(selected_stats,factor_tbl,Dstudy_range,Dstepsz);
% icc_summary=0; var_comp_mean=0;

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


