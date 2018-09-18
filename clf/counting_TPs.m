function [data]=counting_TPs(task,data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% before starting, run: sshfs -o IdentityFile=~/.ssh/MRCInstance1.pem
%    ec2-user@52.87.169.145:data/hcpTask/ mnt/
% task can be: SOCIAL; WM; GAMBLING; RELATIONAL; EMOTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SETUP

% User-defined
group_size_subset='20'; 
suffix_CC='__randomise'; % or FLAME
% suffix_temp='__temp';
suffix_temp='';
cloud_dir=strcat('/Users/stephanie/Documents/data/mnt/');
out_dir=strcat('/Users/stephanie/Google Drive 2/Academic/Lab/Steph - Lab/cluster failure TPR/images/hcp/esz_summary',suffix_CC,suffix_temp,'/');
nperms=500;

% plot param
adorn_plot=1;
adaptive_binning=0;
interval=0.1;
nbins=75; % for histograms; redefined below for EMOTION and RELATIONAL to avoid empty (NAN) bins
ax_xmin=-2.5; ax_xmax=2.5; ax_ymin=0; ax_ymax_esz=0.15; ax_ymax_tp=100;
fontsz=25;

% thresholds
thresh_high=0.8;
thresh_low=0.5;

mkdir(out_dir)

% cope - task pairs
%  FLAME nperms: % soc: 78, WM: ?, gamb: 100, rel: 73,  emot: 100
switch task
    case 'SOCIAL'
        cope='SOCIAL_cope6';
        group_size_full='484';
    case 'WM'
        cope='WM_cope20';
        group_size_full='493';
    case 'GAMBLING'
        cope='GAMBLING_cope6';
        group_size_full='491'; 
    case 'RELATIONAL'
        cope='RELATIONAL_cope4'; 
        group_size_full='480';
    case 'EMOTION'
        cope='EMOTION_cope3';
        group_size_full='482';
    otherwise
        fprintf('No task specified.');
        return
end

root_dir=strcat(cloud_dir,cope);
out_prefix=strcat(out_dir,cope);
root_dir_full=strcat(root_dir,'/GroupSize',group_size_full);
root_dir_sub=strcat(root_dir,'/GroupSize',group_size_subset,suffix_CC,'/Summary');
logfile=strcat(out_prefix,'_thresh_esz.txt');

%% GET DATA
if ~ exist('data','var')
    esz_full=load_nii(strcat(root_dir_full,'/dcoeff.nii.gz'));
    pos_tp_counts=load_nii(strcat(root_dir_sub,'/all_clusters_Pos_sum.nii.gz'));
    neg_tp_counts=load_nii(strcat(root_dir_sub,'/all_clusters_Neg_sum.nii.gz'));
    mask=esz_full.img~=0;
    esz_full=esz_full.img(mask);
    pos_tp_counts=pos_tp_counts.img(mask);
    neg_tp_counts=neg_tp_counts.img(mask);
    data=[esz_full(:), pos_tp_counts(:)*100/nperms, neg_tp_counts(:)*100/nperms];
    data(:,4)=data(:,2)+data(:,3);
end
nvox=size(data,1);
    
%% SUMMARIZE

% bin detections
if adaptive_binning
    interval=(max(data(:,1))-min(data(:,1)))/nbins;
    bin_edges=[min(data(:,1)):interval:max(data(:,1))];
    bin_edges=[bin_edges,max(data(:,1))+interval]; % to include last interval
else
    bin_edges=linspace(ax_xmin,ax_xmax,nbins+1);
end

for i=1:length(bin_edges)-1
    ids=data(:,1)>=bin_edges(i) & data(:,1)<bin_edges(i+1);
    counts(i,1)=nanmean(data(ids,2));
    sd(i,1)=std(data(ids,2));
    counts(i,2)=nanmean(data(ids,3));
    sd(i,2)=std(data(ids,3));
end
counts(isnan(counts))=0;
sd(isnan(sd))=0;

% threshold vox by esz (# vox < T2; between T1 & T2)
ids1=data(:,1) <= thresh_high & data(:,1) >= -thresh_high;
ids2=abs(data(:,1)) <= thresh_high & abs(data(:,1)) >= thresh_low;
perc_esz_lt=sum(+ids1) * 100 / nvox;
perc_esz_btw=sum(+ids2) * 100 / nvox;

% threshold TPs by esz (avg # TPs < T2; btw T1 & T2)
% first_thresh=find(bin_edges>thresh_low,1);
perc_tp_lt = ( sum(data(ids1, 4)) * 100 / sum(+ids1) ) / nperms;
perc_tp_btw = ( sum(data(ids2, 4)) * 100 / sum(+ids2) ) / nperms;

%% PLOT & LOG

% equalize axes
% draw yellow square highlight
% thin bins - replace w probability distr?

% Plot esz
if adaptive_binning
    h=histogram(data(:,1),nbins,'Normalization','probability');
else
    h=histogram(data(:,1),bin_edges,'Normalization','probability');
end
hold on;
plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts/nvox)
hold off;

% h=histfit(data(:,1),nbins,'kernel');
% delete(h(1))
% ytl=yticks*100/nvox;
% for i=1:length(ytl)
%     ytl2{i}=num2str(ytl(i),'%2.1f');
% end
% yticklabels(ytl2)

% add stuff
if adorn_plot
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_esz])
    set(gca,'fontsize',fontsz)
    % highlight
    hold on
    rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_esz],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
    hold off
end

saveas(gcf,strcat(out_prefix,'_esz_hist'),'png')

% Plot TPs by esz
plot_type=1;
if plot_type==1
    % one way (best)
    figure
    histogram('BinEdges',bin_edges,'BinCounts',counts(:,1))
    hold on
    histogram('BinEdges',bin_edges,'BinCounts',counts(:,2))
    hold off
elseif plot_type==2
    % another way
    figure
    scatter(data(:,1),data(:,2),'.')
    hold on
    scatter(data(:,1),data(:,3),'.')
    hold off
elseif plot_type==3
    % a third way
    figure
    errorbar(bin_edges(1:end-1),counts(:,1),sd(:,1));
    hold on
    errorbar(bin_edges(1:end-1),counts(:,2),sd(:,2));
    hold off
end

% add stuff
if adorn_plot
    axis([ax_xmin,ax_xmax,ax_ymin,ax_ymax_tp])
    set(gca,'fontsize',fontsz)
    % add trace of previous hist
    hold on
    scaling=ax_ymax_tp/ax_ymax_esz;
    plot(h.BinEdges(1:end-1)+ h.BinWidth/2,h.BinCounts*scaling/nvox,'--','LineWidth',2)
    % highlight
    rectangle('Position',[-thresh_high,ax_ymin,2*thresh_high,ax_ymax_tp],'FaceColor',[1 1 0 0.2],'EdgeColor','none')
    hold off
end

saveas(gcf,strcat(out_prefix,'_perc_esz_detected'),'png')

% log perc esz and TP within threshold
fid=fopen(logfile,'w');
fprintf(fid,'Percent less than %1.1f: %f; ALSO greater than %1.1f: %f',thresh_high,perc_esz_lt,thresh_low,perc_esz_btw);
fprintf(fid,'\nAvg percent detected under %1.1f: %f; ALSO greater than %1.1f: %f',thresh_high,perc_tp_lt,thresh_low,perc_tp_btw);
fprintf(fid,'\n(%1.0f total permutations)',nperms);
fclose(fid);

%% FINISH

fprintf('Results saved to %s...\n',out_prefix)
close all

end