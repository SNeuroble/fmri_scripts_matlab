function [stats,means]=imgcorr_within_v_between(cmat)
% cmat is correlation between every run and every other run
% right now, this just returns correlations within and between

nruns_per_sub=4*6;
nsub=size(cmat,1)/nruns_per_sub;
dims=size(cmat);

mask_w=zeros(size(cmat)); % within-sub mask
for sub=1:nsub
    id=(sub-1)*nruns_per_sub;
    mask_w((id+1):(id+nruns_per_sub),(id+1):(id+nruns_per_sub))=1;
end

mask_b=ones(size(cmat))-mask_w; % btw-sub mask

% do only tril
mask_w=logical(tril(mask_w,-1));
mask_b=logical(tril(mask_b,-1));

%% get mean vals
edges_w=cmat(mask_w);
edges_b=cmat(mask_b);
edges_tmp_b=edges_b(randperm(length(edges_b),length(edges_w)));

figure
histogram(edges_w,50);
hold on
histogram(edges_tmp_b,50);
hold off

stats=[{'h'},{'p'},{'ci'},{'stats'}];
[stats{2,1},stats{2,2},stats{2,3},stats{2,4}]=ttest2(edges_w,edges_tmp_b);

means=[{'mean within'},{'std within'},{'mean between'},{'std between'}];
means{2,1}=mean(edges_w);
means{2,2}=std(edges_w);
means{2,3}=mean(edges_b);
means{2,4}=std(edges_b);


% % split half - in progress, may never be used
% 
% for sub=1:nsub
%     t=randperm(nruns_per_sub);
%     ids=sub2ind([nruns_per_sub,nruns_per_sub],t(1:nruns_per_sub/2),t(nruns_per_sub/2+1:nruns_per_sub));
%     
%     id=(sub-1)*nruns_per_sub;
%     
%     tmp=cmat((id+1):(id+nruns_per_sub),(id+1):(id+nruns_per_sub));
%     
%     rel=2*r_half/(1+r_half);
% end
% 
% n_replications=rel*(1-r_half)/(r_half*(1-rel));