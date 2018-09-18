function run_GLM_multirun(data,ftbl)
% uses run_GLM, but after averaging over different numbers of runs

% warning('Check that same number of runs for each session. This script does not check this.')

nsubs=max(ftbl(:,1));
nruns=max(ftbl(:,4));
nsess=max(ftbl(:,5));

nfactors=[nsubs nruns nsess];

ftbl_avg=ftbl(ftbl(:,4)==1,:); % single run
ftbl_avgnew=[ftbl_avg(:,1) ftbl_avg(:,5) ftbl_avg(:,2)];

data_dim=size(data{1},1);

for thisnruns=1:max(ftbl(:,4))
    
    for sub=1:nsubs
        for session=1:nsess
            it=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==session);
            data_tmp=data(it);
            data_avg{(sub-1)*4+session}=mean(reshape(cell2mat(data_tmp),[data_dim,data_dim,length(it)]),3);
        end
    end
    
    
    run_GLM(data_avg,ftbl_avgnew);
    
end




