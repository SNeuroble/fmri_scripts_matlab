% create clf imgs for all tasks
% before starting, run: sshfs -o IdentityFile=~/.ssh/MRCInstance1.pem******
%    ec2-user@52.87.169.145:data/hcpTask/ mnt/

out_dir='/Users/steph/Google Drive File Stream/My Drive/Academic/Lab/Steph - Lab/cluster failure TPR/images/emp figs/';
template_nii=strcat(out_dir,'/template_esz_file.mat');
atlas_img='/Users/steph/Google Drive File Stream/My Drive/Academic/Lab/Steph - Lab/Misc/software/data/bioimagesuite/images/shenetal_neuroimage2013/shen_1mm_268_parcellation.nii.gz';

task={'SOCIAL'; 'WM'; 'GAMBLING'; 'RELATIONAL'; 'EMOTION'};
summarize_all_results=1;
do_tfce=0;

for i=1:length(task)
    % cope - task pairs
    %  FLAME nperms: % soc: 78, WM: ?, gamb: 100, rel: 73,  emot: 100
    
    switch task{i}
        case 'SOCIAL'
            task_cope='SOCIAL_cope6';
            group_size_full='484';
        case 'WM'
            task_cope='WM_cope20';
            group_size_full='493';
        case 'GAMBLING'
            task_cope='GAMBLING_cope6';
            group_size_full='491';
        case 'RELATIONAL'
            task_cope='RELATIONAL_cope4';
            group_size_full='480';
        case 'EMOTION'
            task_cope='EMOTION_cope3';
            group_size_full='482';
        otherwise
            fprintf('No task specified.');
            return
    end
    
    % [tp_summary2{i}]=count_TPs(task{i},'doTFCE',1); % ,'data'
    [processed_data{i},res_tpr(:,:,:,i)]=count_TPs(task_cope,group_size_full,'doTFCE',do_tfce,'do_log',0,'make_figs',0); % ,'data'
end

%% Save mean residual maps
if summarize_all_results
    
    if do_tfce
        tfce_suffix='TFCE';
    else
        tfce_suffix='';
    end
    
    % calculate mean residuals    % TODO: T1_2mm_stripped has wrong dim - check stripped for other BIS
    load(template_nii)
    nii=full_file;
    res_avg=mean(res_tpr,4);

    % save mean residual img
    nii.img=res_avg;
    save_nii(nii,strcat(out_dir,'esz_summary__randomise',tfce_suffix,'/MEAN_ALL_resid.nii.gz'));
    
    % average in atlas
    atlas=load_nii(atlas_img);
    atlas=atlas.img;
    atlas=atlas(1:2:end,1:2:end,1:2:end); % subsample to match residual image size
    
    % save atlas residual img
    nii.img=average_within_3d_atlas_simple(res_avg,atlas);
    save_nii(nii,strcat(out_dir,'esz_summary__randomise',tfce_suffix,'/MEAN_ALL_ATLAS_resid.nii.gz'));
    
end



%% plot residual data by networks - testing
% no need to run/load anything before these lines
atlas_nets=load_nii('/Volumes/GoogleDrive/My Drive/Academic/Lab/Steph - Lab/Misc/Software/data/bioimagesuite/images/shenetal_neuroimage2013/shen_1mm_268_parcellation__in_subnetworks.nii.gz');
atlas_nets=atlas_nets(1:2:end,1:2:end,1:2:end); % subsample to match residual image size
n_nets=max(atlas_nets.img(:));

addpath('/Volumes/GoogleDrive/My Drive/Academic/Lab/Steph - Lab/Misc/Software/scripts/Matlab/general/linspecer/')

task={'SOCIAL'; 'WM'; 'GAMBLING'; 'RELATIONAL'; 'EMOTION'};
for i=1:length(task)
    switch task{i}
        case 'SOCIAL'
            task_cope='SOCIAL_cope6';
            group_size_full='484';
        case 'WM'
            task_cope='WM_cope20';
            group_size_full='493';
        case 'GAMBLING'
            task_cope='GAMBLING_cope6';
            group_size_full='491';
        case 'RELATIONAL'
            task_cope='RELATIONAL_cope4';
            group_size_full='480';
        case 'EMOTION'
            task_cope='EMOTION_cope3';
            group_size_full='482';
        otherwise
            fprintf('No task specified.');
            return
    end

    load(['/Volumes/GoogleDrive/My Drive/Academic/Lab/Steph - Lab/cluster failure TPR/images/emp figs/esz_summary__randomise/',task_cope,'_data.mat'])
    emotion_net_ids=atlas_nets.img(mask);

    emotion_resid_nii=load_nii(['/Volumes/GoogleDrive/My Drive/Academic/Lab/Steph - Lab/cluster failure TPR/images/emp figs/esz_summary__randomise/',task_cope,'_resid.nii.gz']);
    emotion_resid=emotion_resid_nii.img(mask);

    % Plot
    ax_xmin=-2.5; ax_xmax=2.5; ax_ymin=0; ax_ymax_tp=30;
    fontsz=25;
    C = linspecer(n_nets+1); 

    figure;
    hold on;
    for this_net=0:n_nets
        % plot
        idx=emotion_net_ids==this_net;
        scatter(data(idx,1),emotion_resid(idx),1,'.','MarkerEdgeColor',C(this_net+1,:))
    end
    plot(data(:,1),zeros(size(data(:,1))),'k-','LineWidth',2) % plot zero residual line
    hold off;

    % prettify
    axis([ax_xmin,ax_xmax,-ax_ymax_tp,ax_ymax_tp])
    set(gca,'fontsize',fontsz)

    % save plot
    %saveas(gcf,strcat(out_prefix,'_esz_v_TPR__residuals'),'png')

end

% %  summary
% % std
% det=[7.025648, 9.442790, 8.419462, 10.355855, 9.626895];
% mean_det=mean(det);
% 
% % tfce
% det=[19.325865, 28.438047, 11.194561, 27.054356, 19.272800];
% mean_det=mean(det);
