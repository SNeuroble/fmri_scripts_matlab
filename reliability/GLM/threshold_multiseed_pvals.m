function [sig_ps] = threshold_multiseed_pvals(ps)
% output is 6xnfactors; 6 p-thresholds

if ismatrix(ps)
    for i=1:size(ps,2)
        ps2{i}=ps(:,i);
    end
    clear ps
    ps=ps2;
end

threshold=[0.1,0.05,0.01,0.005,0.001];

for i=1:length(threshold);
    sig_ps{i+1,1}=sprintf('p<%.3f',threshold(i));
    for j=1:size(ps,2)
        p_mask=(ps{j}<threshold(i));
        sig_ps{1,j+1}=sprintf('level %d',(j));
        sig_ps{i+1,j+1}=p_mask;
    end    
end

end