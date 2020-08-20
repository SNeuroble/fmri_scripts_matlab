function [icc,var,stats,data_avg]=run_rel_variableruns(correctiontype,data,ftbl)
% type is 'none', 'Bonf', 'FDR'
% if you want f-tests, set modeltype to 'linear'
modeltype='linear';

warndlg('Check that same number of runs for each session. This script does not check this.')

nsubs=max(ftbl(:,1));

if size(ftbl,2)==5
    nruns=max(ftbl(:,4));
    nsessions=max(ftbl(:,5));
end

atlasdim=size(data{1});
ftbl_avg=ftbl(ftbl(:,4)==1,:);
data_avg{nsubs*nsessions}=0;

for thisnruns=1:nruns
% for runs=5:6
    
    for sub=1:nsubs
        for thissess=1:nsessions
            it=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==thissess);
            data_tmp=data(it);
            data_avg{(sub-1)*nsessions+thissess}=mean(reshape(cell2mat(data_tmp),[atlasdim,length(it)]),3);
        end
    end

  
    [icc,var,stats,~]=run_reliability(correctiontype,data_avg,ftbl_avg(:,1:3),modeltype);
    
end