function x = get_multiseed_glm_vals(msglm,thisfield)
% get pvals, betas, &c from msglm struct
% options for thisfield: beta, dfe, sfit, s, estdisp, covb: ,se, coeffcorr, t,  p,  resid, residp, residd, resida,  wts
% msglm is a cell vector output from multiseed_analysis ({1:nrois} {2,3} {1:nref_factors} {1}.(thisfield))
% checks for effect coding weirdness

x{length(msglm{1}{2,3})}{length(msglm)}=0;

for i=1:length(msglm)
    for j=1:length(msglm{1}{2,3})
        x{j}{i}=msglm{i}{2,3}{j}{1}.(thisfield);
    end
end

for i=1:length(x)
    x{i}=rearrange_reliability_cellmat(x{i});
end


% if effect coding, reassign multiseed ps to include original and switched reference factor
if length(x)==2
    x{1}{length(x{1})+1}=x{2}{length(x{1})};
end

x=x{1};

end