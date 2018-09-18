function make_icc_visualization(composite_mask,icc_summary,sig_mask,name)
% icc_summary structure is 3x4 cell
% icc_summary{3,1}{1} should not already be masked by sig_mask i.e., both
% matrices must share same dimensions
% e.g., name='pcc'


img=structure_data(composite_mask,icc_summary{3,1}{1}.*(+sig_mask));
img2=img;
img2(img2>0.74)=4;
img2(img2>0.6 & img2<0.74)=3;
img2(img2>0.4 & img2<0.60)=2;
img2(img2>0 & img2<0.40)=1;
img2(img2==0)=0;

save_mat_as_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz',img2,sprintf('/Users/stephanie/Documents/data/traveling_subs/visualization_maps/%s_reliability.nii.gz',name));