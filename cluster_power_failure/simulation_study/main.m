% SAVE YOURSELF THE TROUBLE - I already made this in Python (see ~/Documents/scripts/python/img_corrs.py )
% before starting, mount data2_smn33/cluster_failure/RandomGroupAnalyses/Results/Oulu/6mm/boxcar30_REML/ to "mnt"
% or anything right before "Software"

error('Deprecated. Use ~/Documents/scripts/python/img_corrs.py .')
% Define paths
reference_base_path='~/Documents/data/mnt2/';
smoothness_base_path='~/Documents/data/mnt2/';
target_base_path='~/Documents/data/mnt2/';
mask_base_path='~/Documents/data/sens_spec/';
Software='FSL_Perm';
CDT='2.3';
radius_sq='9';
esz='0.8';
%target_path_suffix1={'esz0.2_noblur','esz0.2';'esz0.5_noblur','esz0.5';'esz1_noblur','esz1';'esz1.5_noblur','esz1.5'};

target_paths=strcat(target_base_path,Software,'/GroupSize20/Twosamplettest/ClusterThreshold',CDT,'/TPR_wFPR/rsq',radius_sq,'/EffectSize',esz,'/NoBlur/Summary/all_clusters_activationadded_sum.nii.gz');

mask_path=strcat(mask_base_path,'TT_N27_3mm_mask.nii'); % copy locally clusterfailure/talairach

reference_paths=strcat(reference_base_path,Software,'/GroupSize20/Twosamplettest/ClusterThreshold',CDT,'/FPR/Summary/all_clusters_sum.nii.gz');

smoothness_path=strcat(smoothness_base_path,'smoothness_allsubs_sph21+tlrc.HEAD');

[corrs,meanFPR]=calc_spatial_similarity(reference_paths,smoothness_path,target_paths,mask_path);


