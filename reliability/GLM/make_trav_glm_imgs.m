function make_trav_glm_imgs(glm_data,mask)
% start w glm_data, e.g.,
%load('rmc_glm_site_GE_11H40M58S.mat','glm_data')
% load('rmc_mask.mat','gmmask_composite')

ps = get_multiseed_ps(glm_data,'p');
[fdr_ps,qs]=fdr_correct_multiseed_pvals(ps);
thresh_ps=threshold_multiseed_pvals(qs);

for i=2:size(thresh_ps,2)      
    struct_data=structure_data(thresh_ps{3,i},gmmask_composite); 
    save_mat_as_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz',struct_data,sprintf('glm_data_tmp%d.nii.gz',i-1)) % visualize
end


% to count them up:
% for i=2:size(thresh_ps_site,2)  
% sum(thresh_ps_site{3,i})
% end