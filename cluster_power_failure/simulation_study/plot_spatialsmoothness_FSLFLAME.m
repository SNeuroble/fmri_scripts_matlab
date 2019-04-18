function [mean_FWHM,residuals]=plot_spatialsmoothness_FSLFLAME(residuals)

addpath('/Users/stephanie/Documents/MATLAB/myscripts/cluster_power_failure/simulation_study/')
data_path='~/Documents/data/mnt/data/cluster_failure/RandomGroupAnalyses/Results/Beijing/6mm/randomEvent_REML/SubjectAnalyses';
% mount FSL2

% Load the anatomical template
% mask = load_nii(['/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz']);
% mask = load_nii(['~/bioimagesuite35/images/FSL_standard/MNI152_T1_2mm_brain_mask.nii.gz']);
% mask = double(mask.img);

% [sy sx sz] = size(mask); % AFNI edit
sx=64; sy=54; sz=50;
st = 20;

voxel_size = 3; % 3 mm voxels

mean_FWHM = zeros(sy,sx,sz);
mean_FWHMx = zeros(sy,sx,sz);
mean_FWHMy = zeros(sy,sx,sz);
mean_FWHMz = zeros(sy,sx,sz);

% files = 1000;

% SMN
subjects=dir([data_path,'/*results']);
subjects={subjects.name};
subjects = cellfun(@(x) x(1:end-8), subjects, 'UniformOutput', 0);
st = length(subjects);
files = 1;

%figure(1)

% Loop over group residuals
% for file = 1:files
%     
%     file
if ~exist('residuals')
    for i = 1:st

        % Load a new group difference map
    %     residuals = load_nii(['/home/andek/Research_projects/RandomGroupAnalyses/Results/FSL/AllGroupResidualsFLAME/group_residuals_Beijing_6mm_Event2_comparison' num2str(file) '.nii.gz']);
    %     residuals = load_nii(file); % SMN

        file=[data_path,'/',subjects{i},'.results/','stats.',subjects{i},'.nii.gz'];

        r = load_nii(file);
        residuals(:,:,:,i) = double(r.img)/10000;

        normalized_residuals(:,:,:,i) = zeros(size(residuals(:,:,:,i)));
    end
end

%% Start
% Normalize each "timeseries"
mask = zeros(sy,sx,sz);
for x = 1:sx
    for y = 1:sy
        for z = 1:sz
%             if (residuals(y,x,z,1) ~= 0) % only check first image
%             if (sum(+(residuals(y,x,z,:) ~= 0)) > 0.25*st) % percent nonzero
%             if (all(residuals(y,x,z,:) ~= 0)) % all nonzero 
                mask(y,x,z) = 1;
                timeseries = squeeze(residuals(y,x,z,:));
                timeseries = timeseries(:);
                normalized_residuals(y,x,z,:) = timeseries / sqrt(timeseries' * timeseries);
%             end
        end
    end
end

spatial_derivatives_xt = zeros(sy,sx,sz,st);
spatial_derivatives_yt = zeros(sy,sx,sz,st);
spatial_derivatives_zt = zeros(sy,sx,sz,st);

% Calculate spatial derivatives
brainvoxels = 0;
allvoxelsmask = zeros(size(mask));
for x = 2:sx-1
    for y = 2:sy-1
        for z = 2:sz-1
            allvoxelsinside = mask(y,x,z) + mask(y+1,x,z) + mask(y-1,x,z) + mask(y,x+1,z) + mask(y,x-1,z) + mask(y,x,z+1) + mask(y,x,z-1);
            if (allvoxelsinside == 7)
                allvoxelsmask(y,x,z) = 1;
                brainvoxels = brainvoxels + 1;
                for subjects = 1:st
                    spatial_derivatives_xt(y,x,z,subjects) = (normalized_residuals(y,x+1,z,subjects) -  normalized_residuals(y,x-1,z,subjects)) / 2;
                    spatial_derivatives_yt(y,x,z,subjects) = (normalized_residuals(y+1,x,z,subjects) -  normalized_residuals(y-1,x,z,subjects)) / 2;
                    spatial_derivatives_zt(y,x,z,subjects) = (normalized_residuals(y,x,z+1,subjects) -  normalized_residuals(y,x,z-1,subjects)) / 2;
                end
            end
        end
    end
end

spatial_derivatives_xt(isnan(spatial_derivatives_xt)) = 0;
spatial_derivatives_yt(isnan(spatial_derivatives_yt)) = 0;
spatial_derivatives_zt(isnan(spatial_derivatives_zt)) = 0;

spatial_derivatives_x = sum(spatial_derivatives_xt.^2,4);
spatial_derivatives_y = sum(spatial_derivatives_yt.^2,4);
spatial_derivatives_z = sum(spatial_derivatives_zt.^2,4);

spatial_derivatives_x_ = sum(spatial_derivatives_xt(:).^2) / brainvoxels;
spatial_derivatives_y_ = sum(spatial_derivatives_yt(:).^2) / brainvoxels;
spatial_derivatives_z_ = sum(spatial_derivatives_zt(:).^2) / brainvoxels;

lambda = [spatial_derivatives_x_ 0 0;
    0 spatial_derivatives_y_ 0;
    0 0 spatial_derivatives_z_];

W = inv(2*lambda + eps*eye(3));
FWHMx_ = sqrt(8*log(2)*W(1,1))*voxel_size
FWHMy_ = sqrt(8*log(2)*W(2,2))*voxel_size
FWHMz_ = sqrt(8*log(2)*W(3,3))*voxel_size
FWHM_ = (FWHMx_*FWHMy_*FWHMz_)^(1/3)

FWHM = zeros(sy,sx,sz);
FWHMx = zeros(sy,sx,sz);
FWHMy = zeros(sy,sx,sz);
FWHMz = zeros(sy,sx,sz);

% Calculate FWHM
for x = 2:sx-1
    for y = 2:sy-1
        for z = 2:sz-1
            if (allvoxelsmask(y,x,z) == 1)

                lambda = [spatial_derivatives_x(y,x,z) 0 0;
                    0 spatial_derivatives_y(y,x,z) 0;
                    0 0 spatial_derivatives_z(y,x,z)];

                W = inv(2*lambda + eps*eye(3));
                FWHMx(y,x,z) = sqrt(8*log(2)*W(1,1));
                FWHMy(y,x,z) = sqrt(8*log(2)*W(2,2));
                FWHMz(y,x,z) = sqrt(8*log(2)*W(3,3));
                FWHM(y,x,z) = (FWHMx(y,x,z)*FWHMy(y,x,z)*FWHMz(y,x,z))^(1/3);

            end
        end
    end
end


mean_FWHM = mean_FWHM + FWHM;
mean_FWHMx = mean_FWHMx + FWHMx;
mean_FWHMy = mean_FWHMy + FWHMy;
mean_FWHMz = mean_FWHMz + FWHMz;

%figure(1)
%imagesc(mean_FWHM(:,:,51) * voxel_size / file); colormap gray; colorbar; title('Smoothness (mm FWHM)'); drawnow



% Summarize all results

mean_FWHM = mean_FWHM/files;
mean_FWHMx = mean_FWHMx/files;
mean_FWHMy = mean_FWHMy/files;
mean_FWHMz = mean_FWHMz/files;

% Show some results
slice = 32;

figure
imagesc(mean_FWHM(:,:,slice) * voxel_size); colormap gray; colorbar; title('Smoothness AFNI OLS (mm FWHM)'); axis image; axis off

% figure
% imagesc(mean_FWHMx(:,:,slice) * voxel_size); colormap gray; colorbar; title('Smoothness in x direction')
% 
% figure
% imagesc(mean_FWHMy(:,:,slice) * voxel_size); colormap gray; colorbar; title('Smoothness in y direction')
% 
% figure
% imagesc(mean_FWHMz(:,:,slice) * voxel_size); colormap gray; colorbar; title('Smoothness in z direction')

% % To save - but WARNING - still must reorient to LPS for comparison w other results (TPR, FPR)
% tmp=load_nii(file);
% tmp.img=mean_FWHM*voxel_size;
% save_nii(tmp,'Mean_FWHM.nii')

end


