function save_mat_as_nii(reference_filename,new_mat,output_filename)
% input has same dimensions as loaded reference file
% example: save_mat_as_nii('/Users/stephanie/bioimagesuite35/images/shenetal_neuroimage2013/fconn_atlas_150_2mm.nii',mat,'temp.nii.gz')
% example 3D: save_mat_as_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz',mat,'temp.nii.gz')
% reference should have same properties as desired output file

reference_img=load_untouch_nii(reference_filename);
output_img=reference_img;
output_img.img=new_mat;
output_img.hdr.dime.datatype=64;
output_img.hdr.dime.bitpix=64;
save_untouch_nii(output_img,output_filename);

end
