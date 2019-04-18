function [struct_data,mask]=structure_data(unstruct_data,mask,varargin)
% restructure data that has been pulled out of a mask
% varargin is for structure type. Use defaults if voxelwise or matrix. Otherwise, if 1-D, specify 'none'.
%
% Ex: from glm data
% load glm_data; load gmmask_composite
% ps = get_multiseed_vals(glm_data,'p'); % another example: unstruct_data=pcc_icc_summary{3,1}{1};
% [fdr_ps,qs]=fdr_correct_multiseed_pvals(ps);
% thresh_ps=threshold_multiseed_pvals(qs);
% for i=2:size(thresh_ps,2)
%    struct_data{i-1}=structure_data(thresh_ps{3,i},gmmask_composite);
% end
% THE FOLLOWING ONLY SAVES FOR SITE 3:
% save_mat_as_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz',struct_data{3},'temp.nii.gz') % visualize

% tips:
% see run_glm line 30 for obtaining gmmask_composite
%   use make_glm_imgs after to convert to structured niis
% for quick view: nii=load_nii('temp.nii.gz'); view_nii(nii)


if ~exist('mask','var') | isempty(mask)
    %assuming want matrix
    ismatrix=1;
    r=int16(roots([1 -1 -2*length(unstruct_data)]));
    mask=ones(r(1),r(1));
    mask=tril(mask,-1);
else
    if size(mask,2)==1 % mask is array; may have to convert to matrix
        if ~isempty(varargin)
            if strcmp(varargin{1},'matrix');
                ismatrix=1;
            elseif strcmp(varargin{1},'none');
                ismatrix=0;
            else error('Incorrect structure type. Must choose ''matrix'' or ''none''.')
            end
        else % if not specified, assume matrix
            ismatrix=1;
        end
            if ismatrix
                warning('Mask is array. Converting to matrix.')

                % compute original matrix size and structure
                % y=x(x-1)/2; y=tril size, x=atlas size
                r=int16(roots([1 -1 -2*length(mask)]));
                mask=ones(r(1),r(1));
                mask=tril(mask,-1);
            end
        
    elseif sum(sum(sum(mask))) ~= length(unstruct_data)
        error('Mask - data dimension mismatch')
    else
        ismatrix=0;
    end
end


t=find(mask);
struct_data=+mask;


for i=1:size(t,1)
    struct_data(t(i))=unstruct_data(i);
end

if ismatrix
    % reflect entries across diag bc matrix
    struct_data=tril(struct_data,-1)+struct_data';
end

% to subsequently mask by threshold: d=data*(+data>0.74);

end