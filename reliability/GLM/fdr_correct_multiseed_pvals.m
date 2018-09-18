function [fdr_ps,qs] = fdr_correct_multiseed_pvals(ps)
% ps are cell vector

for j=2:size(ps,2) % skip first val, intercept
    [fdr_ps(:,j-1),qs(:,j-1)]=mafdr(ps{j}(:));
end
