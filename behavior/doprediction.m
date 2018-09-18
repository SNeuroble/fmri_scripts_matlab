function [r,p,posnegmask,pmdl,nmdl,fullmdl,ppredict,npredict,fullpredict]=doprediction(x_train,y_train,x_test,pthresh,nsubs,nsess,doloo,mdl1,mdl2)
% call from predict_behavior
% can provide previously defined models (mdl1,mdl2) as input
% can convert binarized to weighted--see 3 chunks below

usemultilevelmodel=0; % USER DEFINED - don't need here if mean centering (only wanted to model sep intercepts)
if nsess==1; usemultilevelmodel=0; end

mdl1exists=~isempty(mdl1);
mdl2exists=~isempty(mdl2);


% build model
if ~(mdl1exists && mdl2exists) % if either model doesn't exist, need to rebuild both
    
    % remove missing data
    trainids=~isnan(y_train);
    y_train=y_train(trainids);
    x_train=x_train(:,trainids);
    
    % feature selection - get significantly predictive edges    
    x_train=double(x_train);
    if ~usemultilevelmodel  % build normal corr model
        [r,p]=corr(x_train',y_train);
        % data.connectivity=x_train'; % [r,p]=corr(data.connectivity,data.behav); %tmp
    else fitmymlm
    end
    
    tmp=p<pthresh;
    posnegmask=(+(r>0))-(+(r<0)); % BINARIZED
%     posnegmask=r*100; % WEIGHTED
%     minval=min(abs(posnegmask(tmp))); % WEIGHTED
%     posnegmask(posnegmask>0)=posnegmask(posnegmask>0)-minval; % WEIGHTED
%     posnegmask(posnegmask<0)=posnegmask(posnegmask<0)+minval; % WEIGHTED
    posnegmask=posnegmask.*(+tmp);
    
    % strength of thresholded edges (also separated into pos and neg corrs)
    for i=1:size(x_train,2)
        
        
        % BINARIZED
        pstrength(i)=sum(x_train(posnegmask>0,i));
        nstrength(i)=sum(x_train(posnegmask<0,i));
        fullstrength(i)=pstrength(i)+-nstrength(i); % probably change to GLM
        
%         % WEIGHTED
%         pstrength(i)=sum(x_train(:,i).*(posnegmask.*+(posnegmask>0)));
%         nstrength(i)=sum(x_train(:,i).*(posnegmask.*+(posnegmask<0)));
%         fullstrength(i)=sum(x_train(:,i).*(posnegmask)); % !!!
        
    end
    
    % fit linear model - behavior ~ strength
    pmdl=polyfit(pstrength,y_train',1);
    nmdl=polyfit(nstrength,y_train',1);
    fullmdl=polyfit(fullstrength,y_train',1);
    
else
    
    sub=length(mdl2.pos_mdl);
    sess=length(mdl2.pos_mdl{sub});
    
    pmdl=mdl2.pos_mdl{sub}{sess};
    nmdl=mdl2.neg_mdl{sub}{sess};
    fullmdl=mdl2.full_mdl{sub}{sess};
    posnegmask=mdl1.posnegmask{sub}{sess};
    r=0; p=0;
    
end


 % predict
if doloo
    
    % remove missing data
%     testids=~isnan(y_test);
%     y_test=y_test(testids);
%     x_test=x_test(testids);
    
    % BINARIZED
    pstrength_test=sum(x_test(posnegmask>0));
    nstrength_test=sum(x_test(posnegmask<0));
    fullstrength_test=pstrength_test+-nstrength_test;

%     % WEIGHTED
%     pstrength_test=sum(x_test.*(+posnegmask.*+(posnegmask>0)));
%     nstrength_test=sum(x_test.*(+posnegmask.*+(posnegmask<0)));
%     fullstrength_test=sum(x_test.*(+posnegmask)); % !!!

    ppredict=polyval(pmdl,pstrength_test);
    npredict=polyval(nmdl,nstrength_test);
    fullpredict=polyval(fullmdl,fullstrength_test);
else
    pstrength_test=nan; nstrength_test=nan; fullstrength_test=nan;
    ppredict=nan; npredict=nan; fullpredict=nan;
end

