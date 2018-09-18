function [r,p,pmask,mdl]=cpm_train(x,y,pthresh)
% train model

% feature selection from all train subs
[r,p]=corr(x',y);

pmask=(+(r>0))-(+(r<0)); % BINARIZED
pmask=pmask.*(+(p<pthresh));

% get summary features for each train sub
for i=1:size(x,2)
    summary_feature(i)=mean(x(pmask>0,i))-mean(x(pmask<0,i));
end

% fit behavior to summary features
mdl=polyfit(summary_feature,y',1);
    
    