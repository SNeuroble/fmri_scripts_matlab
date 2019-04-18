function [corrs,meanFPR]=calc_spatial_similarity(reference_paths,smoothness_path,target_paths,mask_path)

% Get mask
mask=load_nii(mask_path);
mask=mask.img>0;

% Get + mask reference data
for i=1:length(reference_paths)
	ref{i}=load_nii(reference_paths{i});
	ref{i}=double(ref{i}.img(mask));
end

% Get + mask smoothness
% Note: Vm_v has 21 elements (7x3) - 3 directions (xyz), 7 subbriks (0: full fstat, 1: model_coefficient, 2: model_tstat, 3: model_fstat, 4: activity_coeff, 5: activity_tstat, 6: activity_fstat)
addpath(genpath('~/Documents/MATLAB/fmri/afni_matlab'))
[~, Vm_v, Info, ~] = BrikLoad (smoothness_path, 'matrix');
Vm_v=Vm_v(:,:,:,4:6); % choose only model coefficients - see above - TODO: is this what we want???
Vm_v=mean(Vm_v,4); % avg over 3 
ref{i+1}=Vm_v(mask);

nref_imgs=length(ref);
ntarget_imgs=length(target_paths(:));
corrs=zeros(ntarget_imgs,nref_imgs); % 3 from refs

% Calculate correlations
for it=1:ntarget_imgs
    a=load_nii(target_paths{i});
    a=double(a.img(mask));
    for i=1:nref_imgs
        corrs(it,i)=corr(a,ref{i});
    end
    meanFPR(it)=mean(a);
end


% Plot
%xmin=0; xmax=target_dim2+1;
%ymin=0; ymax=0.7;

%figure;
%corrs=corrs'; % TODO: I think this is necessary for this script
%for i=1:nref_imgs
%    subplot(1,nref_imgs,i)
%    bar(reshape(corrs(:,i),target_dim1,target_dim2)')
%    axis([xmin xmax ymin ymax])
%end
