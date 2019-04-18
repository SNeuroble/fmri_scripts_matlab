function [data,residual_img_tpr]=count_TPs(task_cope,group_size_full,varargin)

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
cloud_dir__default='/Users/steph/Documents/data/mnt/data/hcpTask/';
out_dir_prefix__default='/Users/steph/Google Drive File Stream/My Drive/Academic/Lab/Steph - Lab/cluster failure TPR/images/emp figs/esz_summary';

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
template_nii=strcat(d,'/template_esz_file.mat');
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
        data(:,4)=data(:,2)+data(:,3); % this is replaced in line 105 (Better combined TPR)
        save(data_file,'data','mask');
    end
end
nvox=size(data,1);

%% BETTER COMBINED TPR
data(data(:,1)<0,2)=0;
data(data(:,1)>0,3)=0;
data(:,4)=data(:,2)+data(:,3); % data(:,1)=esz, data(:,2)=TPR

%% SUMMARIZE

% create bins
if adaptive_binning % set limits based on min and max esz
    interval=(max(data(:,1))-min(data(:,1)))/nbins;
    bin_edges=[min(data(:,1)):interval:max(data(:,1))];
    bin_edges=[bin_edges,max(data(:,1))+interval]; % to include last interval
else
    bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
end

% bin TPR - might not need
for i=1:length(bin_edges)-1
    ids=data(:,1)>=bin_edges(i) & data(:,1)<bin_edges(i+1);
    counts(i)=nanmean(data(ids,4));
    sd(i)=std(data(ids,4));
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
[y_predict_tpr,res_tpr,~]=fit_spline(data(:,1),data(:,4),spline_smoothing,strcat(out_prefix,'_esz_v_TPR'));
load(data_file,'mask');
residual_img_tpr=structure_data(res_tpr,mask);

% save residual img
if ~ exist(strcat(out_prefix,'_resid.nii.gz'), 'file')
    load(template_nii)
    nii=full_file;
    nii.img=residual_img_tpr;
    save_nii(nii,strcat(out_prefix,'_resid.nii.gz'));
end
    

%% PLOT

if make_figs
    
    % 1. Plot esz hist
    figure;
    if adaptive_binning
        h=histogram(data(:,1),nbins,'Normalization','probability');
    else
        h=histogram(data(:,1),bin_edges,'Normalization','probability');
    end
    
    % add stuff: curve & rectangular highlight
    hold on;
    plot(h.BinEdges(1:end-1) + h.BinWidth/2, h.BinCounts/nvox)
    rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_esz],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
    hold off

    % prettify plot
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
    set(gca,'fontsize',fontsz)
    
    % save
    saveas(gcf,strcat(out_prefix,'_esz_hist'),'png')
    
    
    % 2. Plot esz v. TPR
    
    plot_type=1;
    fig=figure;
    new_palette=[[45,170,180]/256; [255,166,0]/256];
    set(fig,'defaultAxesColorOrder',new_palette);
    if plot_type==1
        % one way, scatterplot w fitted line (best)
        hold on
        [~,ind]=sort(data(:,1));
        yyaxis left
        scatter(data(:,1),data(:,4),1,'.') %or lighter: [57,213,185]
        plot(data(ind,1),y_predict_tpr(ind),'-','Color',[0, 96, 113]/256,'LineWidth',2) % [0, 80, 112]
        hold off
    elseif plot_type==2
        % another way - bar plots (second best)
        histogram('BinEdges',bin_edges,'BinCounts',counts)
    elseif plot_type==3
        % a third way - error bar plots
        errorbar(bin_edges(1:end-1),counts,sd);
    end
    
    % add stuff: curve of previous hist & rectangular highlight
    hold on
    yyaxis right
    plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts/nvox,'--','LineWidth',2)
    rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_tp],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
    hold off
        
    % prettify plot
    yyaxis left
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
    yyaxis right
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
    set(gca,'fontsize',fontsz)

    % save
    saveas(gcf,strcat(out_prefix,'_esz_v_TPR'),'png')
    
    
    % 3. Plot residuals

    figure;
    scatter(data(:,1),res_tpr,1,new_palette(1,:),'.')
    
    % highlight points less or greater than 2 std
    hold on;
    n_std=2;
    std_thresh_pos=n_std*std(res_tpr);
    idx=res_tpr>std_thresh_pos | res_tpr<(-std_thresh_pos);
    scatter(data(idx,1),res_tpr(idx),1,'.')
    plot(data(:,1),zeros(size(data(:,1))),'k-','LineWidth',2) % plot zero residual line
    hold off;
    
    % prettify
    axis([ax_xmin,ax_xmax,-ax_ymax_tp,ax_ymax_tp])
    set(gca,'fontsize',fontsz)
    
    
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
