% Plot by group

data_folder='/Users/stephanie/Google Drive 2/Academic/Lab/Steph - Lab/cluster failure TPR/images/emp figs/esz_summary__randomise';
task_cope='EMOTION_cope3';
pos=true;
if pos; pos_string='pos'; else; pos_string='neg'; end
data_prefix=strcat(data_folder,'/',task_cope);
data_file=strcat(data_prefix,'_data.mat');
img_file=strcat(data_prefix,'_',pos_string,'_resid.nii.gz');
assignments_file=strcat(data_prefix,'_',pos_string,'_group_assignments.txt');

load(data_file,'data','mask');
full_file=load_nii(img_file);
assignments = csvread(assignments_file);

% get mapping
full_file.img(mask)=assignments;
save_nii(full_file, strcat(data_prefix,'_',pos_string,'_group_assignments.nii.gz'));


