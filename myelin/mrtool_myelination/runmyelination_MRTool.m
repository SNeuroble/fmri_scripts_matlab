
% define templates
T1template={'/Users/stephanie/Documents/MATLAB/fmri/spm12/toolbox/MRTool/template/mni_icbm152_t1_tal_nlin_sym_09a.nii,1'};
T2template={'/Users/stephanie/Documents/MATLAB/fmri/spm12/toolbox/MRTool/template/mni_icbm152_t2_tal_nlin_sym_09a.nii,1'};
tpmimage={'/Users/stephanie/Documents/MATLAB/fmri/spm12/toolbox/MRTool/template/TPM.nii,1'};
    
% define data directories
% sshfs smn33@172.23.202.134:smn33_data/test_retest/ Documents/data/mnt
datadir='/Users/stephanie/Documents/data/mnt';
outdir=sprintf('%s/myelin_MRTool',datadir);

% define scans
% removed 10 and 11 for missing data
scans={'TRT001_1_TA';
'TRT001_2_TB';
'TRT001_3_TB';
'TRT001_4_TA';
'TRT002_1_TB';
'TRT002_2_TA';
'TRT002_3_TB';
'TRT002_4_TA';
'TRT003_1_TA';
'TRT003_2_TA';
'TRT003_3_TB';
'TRT003_4_TB';
'TRT004_1_TA';
'TRT004_2_TB';
'TRT004_3_TA';
'TRT004_4_TB';
'TRT005_1_TA';
'TRT005_2_TB';
'TRT005_3_TB';
'TRT005_4_TA';
'TRT006_1_TB';
'TRT006_2_TB';
'TRT006_3_TA';
'TRT006_4_TA';
'TRT007_1_TB';
'TRT007_2_TB';
'TRT007_3_TA';
'TRT007_4_TA';
'TRT008_1_TA';
'TRT008_2_TA';
'TRT008_3_TB';
'TRT008_4_TB';
'TRT009_1_TA';
'TRT009_2_TA';
'TRT009_3_TB';
'TRT009_4_TB';
'TRT012_1_TA';
'TRT012_2_TA';
'TRT012_3_TB';
'TRT012_4_TB'};

% setupfile={'/Users/stephanie/Documents/MATLAB/myscripts/mrtool_myelination/myelination_setupfile_MRTool.m'};


% quick file check
T2test=sprintf('%s/%s/%s_stack3d_SPC_crop.nii',datadir,scans{1},scans{1});
if ~exist(T2test,'file'); error(sprintf('File %s missing. Troubleshooting: (1) Check niftis are gunzipped, (2) Check connection, (3) Check paths.\n',T2test)); end

for i=4:length(scans) 
    thisscan=scans{i};
    
    T1img={sprintf('%s/%s/%s_stack3d_MPR_crop.nii,1',datadir,thisscan,thisscan)};
    T2img={sprintf('%s/%s/%s_stack3d_SPC_crop.nii,1',datadir,thisscan,thisscan)};
    thisoutdir={sprintf('%s/%s/',outdir,thisscan)};
    mkdir(thisoutdir{1})
    thisoutdir={sprintf('%s',thisoutdir{1})};
    
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.t1w = T1img;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.t2w = T2img;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.res_dir = thisoutdir;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.t1_template = T1template;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.t2_template = T2template;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.tpm_file = tpmimage;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasreg_t1 = 0.0001;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasfwhm_t1 = 60;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasreg_t2 = 0.0001;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.biasfwhm_t2 = 60;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.calibration_flag = 1;
    matlabbatch{1}.spm.tools.MRI.T1T2_calc.DelIntermediate = 1;
    
    spm('defaults', 'FMRI');
    fprintf('Done initializing inputs. Starting job. \n')
    
    spm_jobman('run', matlabbatch);
%     spm_jobman('run', thisscript);
end



% original generated script:
%jobfile = {'/Users/stephanie/Documents/MATLAB/myscripts/mrtool_myelination/myelination_setupfile_MRTool.m'};
%jobs = repmat(jobfile, 1,1);
%inputs = cell(0, 1);
%spm('defaults', 'FMRI');
%spm_jobman('run', jobs, inputs{:});
