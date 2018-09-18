function run_GLM(data,ftbl)

% IMPORTANT NOTE: previous version used fdr_ps, but should have used qs.
% Results nearly identical but technically q's should be used.

% look for saved output

% USER-DEFINED
% either manually replace or use pwd (after cd to datadir)
procedure='tmpdata';
% datadir=sprintf('/Users/stephanie/Documents/data/traveling_subs/results_%s/',procedure);
% [data,ftbl]=load_trav_files(datadir,'*nii*');
% % [data,ftbl]=load_trav_files('/Users/stephanie/Documents/data/traveling_subs/results_matrix','*.nii.gz');
% Note: if want to use any other data besides image or matrix, edit below and must make into a vector 


%% main

sz=size(data{1},1);
nd=ndims(data{1});

if(nd==3)
    doing_voxelwise=1; % note: if doing matrix, don't load mask (voxelwise=0)
elseif (nd==2)
    doing_voxelwise=0; % note: if doing matrix, don't load mask (voxelwise=0)
else
    error('weird dimensions')
end

if(doing_voxelwise) % mask bc voxels
    gmmask = load_untouch_nii('/Users/stephanie/Documents/data/traveling_subs/5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz');
    gmmask=gmmask.img;
    if isa(gmmask,'uint8') % if mask is a uint8, convert to single
        im2single(gmmask)
    end
    gmmask(gmmask<3)=0;
    gmmask(gmmask==3)=1;
    [masked_data,gmmask_composite]=get_masked_data_special_trav(data,gmmask);
    % not doing sigedgemask here
    
else % no gmmask bc matrix; apply trimask instead
    if sz>1
        for i=1:length(data)
    %         masked_data{i}=data{i}(:);
            masked_data{i}=data{i}(logical(tril(ones(sz,sz),-1)));
        end
    else
        masked_data=data;
    end
end

% if just want composite mask, uncomment these:
% save('degree_unnorm_gmmask_composite.mat', 'gmmask_composite')
% break;


data_rearr=rearrange_reliability_cellmat(masked_data);
[factor_tbl{1:length(data_rearr)}]=deal(ftbl);

factors={'subj','site','day'};
refcategories={'all'};
factors={'site'};
refcategories={'GE'};

for thisfactor=factors
    for thisref=refcategories
       
        % don't do if 'GE' unless also 'site'
        if strcmp(thisref{1},'GE') && (strcmp(thisfactor{1},'subj') || strcmp(thisfactor{1},'day'))
        else
            glm_data = multiseed_GLM(@fit_trav_glm,data_rearr,factor_tbl,thisref{1},thisfactor{1});
            ps = get_multiseed_glm_vals(glm_data,'p');
            
            if sz > 1
                [fdr_ps,qs]=fdr_correct_multiseed_pvals(ps);
                thresh_ps=threshold_multiseed_pvals(qs);
            else
                thresh_ps=threshold_multiseed_pvals(cell2mat(ps));
            end
            
%             bs = get_multiseed_betas(glm_data); mm=do_for_sig(bs,fdr_ps,'var');

            % count # q's <0.05
            counts=0;
            thisthresh=3; % 2: p<0.1; 3: p<0.05
            for i=2:size(thresh_ps,2)
                counts(i-1)=sum(thresh_ps{thisthresh,i});
            end
            
            clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
            save(sprintf('%s_glm_%s_%sH%sM%sS.mat',procedure,thisfactor{1},clockinfo{1},clockinfo{2},clockinfo{3}),'glm_data','counts')
        end
        
    end

end
