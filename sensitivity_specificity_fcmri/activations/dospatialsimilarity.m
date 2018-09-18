%dospatialsimilarity
% extract in TT mask
% mount data2_smn33/cluster_failure/RandomGroupAnalyses/Results/Oulu/6mm/boxcar30_REML/ to "mnt2"

referencepath='~/Documents/data/mnt2/finalimgs/GroupAnalyses_';

a_paths={'esz0.5_noblur','esz0.2','esz0.5_noblur','esz0.5','esz1_noblur','esz1','esz1.5_noblur','esz1.5'};
a_pathending='.nii.gz';
nblur=2;
nesz=length(a_paths)/nblur;

mask_path='~/Documents/data/sens_spec/TT_N27_3mm_mask.nii'; % copy locally clusterfailure/talairach

reference1_path=sprintf('%s%s',referencepath,'activationonly.nii.gz');
reference2_path=sprintf('%s%s',referencepath,'FPR_norm10000.nii.gz');

mask=load_nii(mask_path);
ref{1}=load_nii(reference1_path);
ref{2}=load_nii(reference2_path);

mask=mask.img>0;
ref{1}=double(ref{1}.img(mask));
ref{2}=double(ref{2}.img(mask));


% load smoothness
% Vm_v is 21 long (7x3) - 3 directions (xyz), 7 subbriks (0: full fstat, 1: model_coefficient, 2: model_tstat, 3: model_fstat, 4: activity_coeff, 5: activity_tstat, 6: activity_fstat)
addpath(genpath('../fmri/afni_matlab'))
reference3_path='../../data/mnt2/finalimgs/smoothness_allsubs_sph20+tlrc.HEAD';
[~, Vm_v, Info, ~] = BrikLoad (reference3_path, 'matrix');
Vm_v=Vm_v(:,:,:,4:6); % choose only model coefficients - see above
Vm_v=mean(Vm_v,4); % avg over 3 
ref{3}=Vm_v(mask);

ncomparisons=length(ref);
corrs=zeros(length(a_paths),ncomparisons); % 3 from refs

% docorrs
for it=1:length(a_paths)
    
    thisdatapath=sprintf('%s%s%s',referencepath,a_paths{it},a_pathending);
    
    a=load_nii(thisdatapath);
    a=double(a.img(mask));

    for i=1:ncomparisons
        corrs(it,i)=corr(a,ref{i});
    end
    
    meanFPR(it)=mean(a);

end


% plot
xmin=0; xmax=nesz+1;
ymin=0; ymax=0.7;

figure;
for i=1:ncomparisons
    subplot(1,ncomparisons,i)
    bar(reshape(corrs(:,i),nblur,nesz)')
    axis([xmin xmax ymin ymax])
end
