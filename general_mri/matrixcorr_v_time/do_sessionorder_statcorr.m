function corrmat=do_sessionorder_statcorr(data,ftbl)
% shows whether session order significantly correlated w correlation val

study='trav';

if size(data{1},1)>1
    doconn=1; 
else
    doconn=0;
end

atlas_len=size(data{1},1);

nsubs=max(ftbl(:,1));


%% Trav or TRT
if strcmp(study,'trav')
    load('trav_sessions.mat','session_order')
    
    nsessions=max(ftbl(:,2)); % this is site, which == nsessions
    
    trimask=logical(tril(ones(size(data{1})),-1));
    
    for sub=1:nsubs
        
        % get mean over all runs, then reorder
        for session=1:nsessions
            it=find(ftbl(:,1)==sub & ftbl(:,2)==session_order(sub,session) & ftbl(:,3)<3);
            order_ids(sub,2*session-1:2*session)=it;
            
            data_tmp=data(it);
            data_avg{sub,session}=mean(reshape(cell2mat(data_tmp),[atlas_len,atlas_len,length(it)]),3);
        end
        
        if doconn
             for s1=1:nsessions
                 for s2=1:nsessions
                     corrmat{sub}(s1,s2)=corr(data_avg{(sub),s1}(trimask),data_avg{(sub),s2}(trimask));
                 end
             end            
        else
            tmp=cell2mat(data_avg(sub,1:nsessions));
            mdl=fitlm(1:8,tmp);
            corrmat(1,sub)=mdl.Coefficients.pValue(2);
            corrmat(2,sub)=mdl.Rsquared.Ordinary;
            corrmat(3,sub)=mdl.Coefficients.Estimate(2);
        end
        
        
    end
    
    

    
    
elseif strcmp(study,'trt')
     warndlg('Check that same number of runs for each session. This script does not check this.')
     
     nruns=max(ftbl(:,4));
     nsessions=max(ftbl(:,5));
     
     ftbl_avg=ftbl(ftbl(:,4)==1,:);
     
     trimask=logical(tril(ones(size(data{1})),-1));
     
     for sub=1:nsubs
         for session=1:nsessions
             it=find(ftbl(:,1)==sub & ftbl(:,4)<=nruns & ftbl(:,5)==session);
             data_tmp=data(it);
             data_avg{(sub-1)*4+session}=mean(reshape(cell2mat(data_tmp),[atlas_len,atlas_len,length(it)]),3);
         end
         
         for s1=1:nsessions
             for s2=1:nsessions
                 corrmat{sub}(s1,s2)=corr(data_avg{(sub-1)*4+s1}(trimask),data_avg{(sub-1)*4+s2}(trimask));
             end
         end
         
    end
    
end
