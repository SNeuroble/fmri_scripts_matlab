%-----------------------------------------------------------------------
% Job saved on 10-Apr-2017 16:22:04 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------



matlabbatch{1}.spm.tools.MRI.T1T2_calc.t1w = {'/Users/stephanie/Documents/data/mnt/TRT001_1_TA_stack3d_MPR_crop.nii,1'};
matlabbatch{1}.spm.tools.MRI.T1T2_calc.t2w = {'/Users/stephanie/Documents/data/mnt/TRT001_1_TA_stack3d_SPC.nii,1'};
matlabbatch{1}.spm.tools.MRI.T1T2_calc.res_dir = {'/Users/stephanie/Documents/MATLAB/workspaces/testretest/myelination/volume/TRT001_1'};
matlabbatch{1}.spm.tools.MRI.T1T2_calc.t1_template = {'/Users/stephanie/Documents/MATLAB/fmri/spm12/toolbox/MRTool/template/mni_icbm152_t1_tal_nlin_sym_09a.nii,1'};
matlabbatch{1}.spm.tools.MRI.T1T2_calc.t2_template = {'/Users/stephanie/Documents/MATLAB/fmri/spm12/toolbox/MRTool/template/mni_icbm152_t2_tal_nlin_sym_09a.nii,1'};
matlabbatch{1}.spm.tools.MRI.T1T2_calc.tpm_file = {'/Users/stephanie/Documents/MATLAB/fmri/spm12/toolbox/MRTool/template/TPM.nii,1'};
matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasreg_t1 = 0.0001;
matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasfwhm_t1 = 60;
matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasreg_t2 = 0.0001;
matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasfwhm_t2 = 60;
matlabbatch{1}.spm.tools.MRI.T1T2_calc.calibration_flag = 1;
matlabbatch{1}.spm.tools.MRI.T1T2_calc.DelIntermediate = 1;


% can't use environment variables bc this script is passed into another

% matlabbatch{1}.spm.tools.MRI.T1T2_calc.t1w = T1img;
% matlabbatch{1}.spm.tools.MRI.T1T2_calc.t2w = T2img;
% matlabbatch{1}.spm.tools.MRI.T1T2_calc.res_dir = thisoutdir;
% matlabbatch{1}.spm.tools.MRI.T1T2_calc.t1_template = T1template;
% matlabbatch{1}.spm.tools.MRI.T1T2_calc.t2_template = T2template;
% matlabbatch{1}.spm.tools.MRI.T1T2_calc.tpm_file = tpmimage;

