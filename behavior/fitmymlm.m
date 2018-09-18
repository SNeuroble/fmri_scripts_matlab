% only necessary if it gets called from doprediction
% relies heavily on variable definitions there

% format data into table with variable names - necessary for MLM
if doloo; nsubs_leftin=nsubs-1;
else nsubs_leftin=nsubs; end
nsess_leftin=nsess;
t=repmat(1:nsess_leftin,1,nsubs_leftin)';
t2=repmat(1:nsubs_leftin,nsess_leftin,1); t2=reshape(t2,nsess_leftin*nsubs_leftin,1);
data=table(y_train,t2,t,x_train(1,:)');
data.Properties.VariableNames={'behav','subj','sess','connectivity'};

nedges=size(x_train,1);
amt_fin_old=0;

for i=1:nedges % have to fit a separate MLM for each edge
    data.connectivity=x_train(i,:)'; % TODO: speed this up by pre-defining
    mdl=fitlme(data,'behav ~ 1 + sess * connectivity + (1|subj)');
    p(i)=mdl.Coefficients(3,6).pValue;
    r(i)=mdl.Coefficients(3,2).Estimate; % THESE ARE BETAS, NOT R - just named r for consistency with non-multilevel model
    % report percent complete
    if (mod(round(i*100/nedges),10)==0); amt_fin=round(i*100/nedges);
        if(~(amt_fin==amt_fin_old)); amt_fin_old=amt_fin;
            fprintf('edge%0.0f %0.0f%% finished\n',i,amt_fin);
        end
    end
end