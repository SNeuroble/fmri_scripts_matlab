function [data,residual_img_pos,residual_img_neg]=count_TPs__separate_pos_neg(task_cope,group_size_full,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before starting, mount HCP instance: sshfs -o IdentityFile=~/.ssh/MRCInstance1.pem
%    ec2-user@52.87.169.145: mnt/
% Task can be: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION
% This script summarizes TPR v effect size data.
% Summarization: bins voxels by effect size, gets mean TPR within bins,
% and fits spline to effect size vs. mean TPR
% Plot: binned effect size, binned TPR,  d v. TPR spline, d v. TPR residual map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variables

group_size_subset__default='20';
doRandomise__default=1;
doTFCE__default=0;
testing__default=0;
cloud_dir__default='/Users/stephanie/Documents/data/mnt/data/hcpTask/';
out_dir_prefix__default=('/Users/stephanie/Google Drive 2/Academic/Lab/Steph - Lab/cluster failure TPR/images/emp figs/esz_summary');
data_default=[];
make_figs_default=1;
do_log_default=1;

p = inputParser;
addRequired(p,'task_cope',@ischar);
addRequired(p,'group_size_full',@ischar);
addOptional(p,'data',data_default);
addOptional(p,'group_size_subset',group_size_subset__default);
addOptional(p,'doRandomise',doRandomise__default);
addOptional(p,'doTFCE',doTFCE__default);
addOptional(p,'testing',testing__default);
addOptional(p,'cloud_dir',cloud_dir__default);
addOptional(p,'out_dir_prefix',out_dir_prefix__default);
addOptional(p,'make_figs',make_figs_default);
addOptional(p,'do_log',do_log_default);
parse(p,task_cope,group_size_full,varargin{:});

task_cope=p.Results.task_cope;
group_size_full=p.Results.group_size_full;
data=p.Results.data;
group_size_subset=p.Results.group_size_subset;
doRandomise=p.Results.doRandomise;
doTFCE=p.Results.doTFCE;
testing=p.Results.testing;
cloud_dir=p.Results.cloud_dir;
out_dir_prefix=p.Results.out_dir_prefix;
make_figs=p.Results.make_figs;
do_log=p.Results.do_log;

%% SETUP
% Set up dirs
if doRandomise; suffix_CC='__randomise'; else; suffix_CC='__FLAME'; end
if doTFCE; suffix_TFCE='TFCE'; else; suffix_TFCE=''; end
if testing; suffix_temp='__temp'; else; suffix_temp=''; end
out_dir=strcat(out_dir_prefix,suffix_CC,suffix_TFCE,suffix_temp,'/');
nperms=500;

% plot param
adorn_plot=1;
adaptive_binning=0;
interval=0.1;
nbins=75; % for histograms; redefined below for EMOTION and RELATIONAL to avoid empty (NAN) bins
ax_xmin=-2.5; ax_xmax=2.5; ax_ymin=0; ax_ymax_esz=0.15; ax_ymax_tp=100;
fontsz=25;
spline_smoothing=0.995;

% thresholds
thresh_high=0.8;
thresh_low=0.5;

if ~ exist(out_dir, 'dir'); mkdir(out_dir); end

root_dir=strcat(cloud_dir,task_cope);
out_prefix=strcat(out_dir,task_cope);
root_dir_full=strcat(root_dir,'/GroupSize',group_size_full);
root_dir_sub=strcat(root_dir,'/GroupSize',group_size_subset,suffix_CC,suffix_TFCE,'/Summary');
logfile=strcat(out_prefix,'_thresh_esz.txt');

%% GET DATA
data_file=strcat(out_prefix,'_data.mat');
[d,~,~] = fileparts(out_dir_prefix);
sample_esz_file=strcat(d,'/sample_esz_file.mat');
if isempty(data)
    if exist(data_file, 'file') == 2
        load(data_file,'data');
    else
        esz_full=load_nii(strcat(root_dir_full,'/dcoeff.nii.gz'));
        pos_tp_counts=load_nii(strcat(root_dir_sub,'/all_clusters_Pos_sum.nii.gz'));
        neg_tp_counts=load_nii(strcat(root_dir_sub,'/all_clusters_Neg_sum.nii.gz'));
        mask=esz_full.img~=0;
        esz_full=esz_full.img(mask);
        pos_tp_counts=pos_tp_counts.img(mask);
        neg_tp_counts=neg_tp_counts.img(mask);
        data=[esz_full(:), pos_tp_counts(:)*100/nperms, neg_tp_counts(:)*100/nperms];
        data(:,4)=data(:,2)+data(:,3);
        save(data_file,'data','mask');
    end
end
nvox=size(data,1);

%% SUMMARIZE

% create bins
if adaptive_binning % set limits based on min and max esz
    interval=(max(data(:,1))-min(data(:,1)))/nbins;
    bin_edges=[min(data(:,1)):interval:max(data(:,1))];
    bin_edges=[bin_edges,max(data(:,1))+interval]; % to include last interval
else
    bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
end

% bin TPR
for i=1:length(bin_edges)-1
    ids=data(:,1)>=bin_edges(i) & data(:,1)<bin_edges(i+1);
    counts(i,1)=nanmean(data(ids,2));
    sd(i,1)=std(data(ids,2));
    counts(i,2)=nanmean(data(ids,3));
    sd(i,2)=std(data(ids,3));
end
counts(isnan(counts))=0;
sd(isnan(sd))=0;

% get voxels within thresholds (vox < thresh_high; vox < thresh_high & > thresh_low)
ids1=data(:,1) <= thresh_high & data(:,1) >= -thresh_high;
ids2=abs(data(:,1)) <= thresh_high & abs(data(:,1)) >= thresh_low;

% calc percent voxels within thresholds
perc_vox_lt_thr_high=sum(+ids1) * 100 / nvox;
perc_vox_btw_thr_high_and_low=sum(+ids2) * 100 / nvox;

% mean TPR within thresholds
tpr_lt_thr_high = ( sum(data(ids1, 4)) / sum(+ids1) );
tpr_btw_thr_high_and_low = ( sum(data(ids2, 4)) / sum(+ids2) ) ;

% mean TPR "at" (around) thresholds
half_bin_width=0.05; % ad hoc bin size
t=[thresh_high, thresh_low];
for i=1:length(t)
    ids_pos=data(:,1) <= (t(i)+half_bin_width) & data(:,1) >= (t(i)-half_bin_width);
    ids_neg=data(:,1) <= (-t(i)+half_bin_width) & data(:,1) >= (-t(i)-half_bin_width);
    tpr_at_thr(3*i-2)=( sum(data(ids_pos, 4)) / sum(+ids_pos) );
    tpr_at_thr(3*i-1)=( sum(data(ids_neg, 4)) / sum(+ids_neg) );
    tpr_at_thr(3*i)=mean(tpr_at_thr((3*i-2):(3*i-1)));
end


%% FIT EFFECT SIZE vs. TPR

% fit spline to d v TPR -> calculate residuals
% if ~ exist(strcat(out_prefix,'_pos_resid.nii.gz'), 'file')
load(data_file,'mask');
[y_predict_pos,res_pos,SSE_pos]=fit_spline(data(:,1),data(:,2),spline_smoothing,strcat(out_prefix,'_esz_v_TPR_pos'));
[y_predict_neg,res_neg,SSE_neg]=fit_spline(data(:,1),data(:,3),spline_smoothing,strcat(out_prefix,'_esz_v_TPR_neg'));


residual_img_pos=structure_data(res_pos,mask);
residual_img_neg=structure_data(res_neg,mask);

%% PLOT

if make_figs
    
    % 1. Plot esz hist
    if adaptive_binning
        h=histogram(data(:,1),nbins,'Normalization','probability');
    else
        h=histogram(data(:,1),bin_edges,'Normalization','probability');
    end
    hold on;
    plot(h.BinEdges(1:end-1) + h.BinWidth/2, h.BinCounts/nvox)
    hold off;
    
    % h=histfit(data(:,1),nbins,'kernel');
    % delete(h(1))
    % ytl=yticks*100/nvox;
    % for i=1:length(ytl); ytl2{i}=num2str(ytl(i),'%2.1f'); end
    % yticklabels(ytl2)
    
    % add stuff to esz hist
    if adorn_plot
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
        set(gca,'fontsize',fontsz)
        % highlight
        hold on
        rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_esz],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
        hold off
    end
    
    % save esz hist
    saveas(gcf,strcat(out_prefix,'_esz_hist'),'png')
    
    
    % 2. Plot TPs by esz (binned)
    plot_type=1;
    if plot_type==1
        % one way, scatterplot w fitted line (best)
        figure
        hold on
        [~,ind]=sort(data(:,1));
        yyaxis left
        scatter(data(:,1),data(:,2),1,'b.')
        scatter(data(:,1),data(:,3),1,'b.')
        plot(data(ind,1),y_predict_pos(ind),'k-','LineWidth',2)
        plot(data(ind,1),y_predict_neg(ind),'k-','LineWidth',2)
        hold off
    elseif plot_type==2
        % another way - bar plots (second best)
        figure
        histogram('BinEdges',bin_edges,'BinCounts',counts(:,1))
        hold on
        histogram('BinEdges',bin_edges,'BinCounts',counts(:,2))
        hold off
    elseif plot_type==3
        % a third way - error bar plots
        figure
        errorbar(bin_edges(1:end-1),counts(:,1),sd(:,1));
        hold on
        errorbar(bin_edges(1:end-1),counts(:,2),sd(:,2));
        hold off
    end
    
    % add stuff to TPs by esz (binned)
    if adorn_plot
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
        set(gca,'fontsize',fontsz)
        % add trace of previous hist
        hold on
        yyaxis right
        axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
        plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts/nvox,'--','LineWidth',2)
%         scaling=ax_ymax_tp/ax_ymax_esz;
%         plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts*scaling/nvox,'--','LineWidth',2)
        % highlight
        rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_tp],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
        hold off
    end
    
    % save TPs by esz (binned)
    saveas(gcf,strcat(out_prefix,'_perc_esz_detected'),'png')
    
    
    % 3. Plot esz v. TPR
    [~,ind]=sort(data(:,1));
    
    figure;
    hold on;
    scatter(data(:,1),data(:,2),1,'b.') % pos
    scatter(data(:,1),data(:,3),1,'b.') % neg
    plot(data(ind,1),y_predict_pos(ind),'k-','LineWidth',2)
    plot(data(ind,1),y_predict_neg(ind),'k-','LineWidth',2)
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
    hold off;

    % save
    saveas(gcf,strcat(out_prefix,'_esz_v_TPR'),'png')
    
    % plot residuals for diagnostics
    figure;
    hold on;
    scatter(data(:,1),res_pos,1,'b.')
    scatter(data(:,1),res_neg,1,'b.')
    % highlight points less or greater than 2 std
    n_std=2;
    % pos residuals
    std_thresh_pos=n_std*std(res_pos);
    idx=res_pos>std_thresh_pos | res_pos<(-std_thresh_pos);
    scatter(data(idx,1),res_pos(idx),1,'.')
    % neg residuals
    std_thresh_neg=n_std*std(res_neg);
    idx=res_neg>std_thresh_neg | res_neg<(-std_thresh_neg);
    scatter(data(idx,1),res_neg(idx),1,'.')
    plot(data(ind,1),zeros(size(data(:,1))),'k-','LineWidth',2) % plot zero residual line
    hold off;
    
%     subplot(3,1,3) % NEG PLOT
%     hold on;
%     scatter(data(:,1),res_neg,1,'.')
%     plot(data(ind,1),zeros(size(data(:,1))),'k-','LineWidth',2)
%     % highlight points less or greater than 2 std
%     n_std=2;
%     std_thresh_neg=n_std*std(res_neg);
%     idx=res_neg>std_thresh_neg | res_neg<(-std_thresh_neg);
%     scatter(data(idx,1),res_neg(idx),1,'.')
%     hold off;

    %[r,p]=corr(residuals,x); % residuals should now be uncorrelated w x

    % save plot
    saveas(gcf,strcat(out_prefix,'_esz_v_TPR__residuals'),'png')

    
end


%% Log Summaries: perc esz and TP at thresholds
if do_log
    fid=fopen(logfile,'w');
    fprintf(fid,'Percent less than d=%1.1f: %f; ALSO greater than %1.1f: %f',thresh_high,perc_vox_lt_thr_high,thresh_low,perc_vox_btw_thr_high_and_low);
    fprintf(fid,'\nAvg percent detected between d=+/-%1.1f: %f; ALSO greater than %1.1f: %f',thresh_high,tpr_lt_thr_high,thresh_low,tpr_btw_thr_high_and_low);
    fprintf(fid,'\nPercent detected at d=+/-%1.1f: %f (+), %f (-), %f (mean) ',thresh_high,tpr_at_thr(1:3));
    fprintf(fid,'\nPercent detected at d=+/-%1.1f: %f (+), %f (-), %f (mean) ',thresh_low,tpr_at_thr(4:6));
    fprintf(fid,'\n(%1.0f total permutations)',nperms);
    fclose(fid);
end

%% FINISH

fprintf('Results saved to %s...\n',out_prefix)
% close all

end
