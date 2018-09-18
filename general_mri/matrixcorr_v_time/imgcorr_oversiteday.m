 function [cmat_mean]=imgcorr_oversiteday(data,ftbl,saveimg)
%do_sessionorder_statcorr is probably similar, but a little more complex

if ~exist('saveimg','var') | isempty(saveimg)
  saveimg = 0;
end

ftbl_avg=ftbl(ftbl(:,4)==1,:);

nsubs=max(ftbl(:,1));
nsites=max(ftbl(:,2));
ndays=max(ftbl(:,3));
nruns=max(ftbl(:,4)); thisnruns=nruns;
nsess=max(ftbl(:,5));

matdim=size(data{1},1);
trimask=logical(tril(ones(matdim,matdim),-1));


% do mean over all runs
for sub=1:nsubs
    for session=1:nsess
        it=find(ftbl(:,1)==sub & ftbl(:,4)<=thisnruns & ftbl(:,5)==session);
        data_tmp=data(it);
        data_avg{(sub-1)*4+session}=mean(reshape(cell2mat(data_tmp),[matdim,matdim,length(it)]),3);
    end
end
    

% re-sort
for s=1:nsubs
    for i=1:nsites
        for j=1:ndays
            it=find(ftbl_avg(:,1)==s & ftbl_avg(:,2)==i & ftbl_avg(:,3)==j);
            data_resort((s-1)*nsites*ndays+((i-1)*nsites+j))=data_avg(it);
            % sorts like this: subj(1 1 1 1) site(1 1 2 2) day(1 2 1 2)
        end
    end
end

% within sub, across site and day correlation
% mean cmat
for s=1:nsubs
    for i=1:nsites*ndays
        for j=1:nsites*ndays
            s1=(s-1)*nsess+i;
            s2=(s-1)*nsess+j;
            cmat{s}(i,j)=corr(data_resort{s1}(trimask),data_resort{s2}(trimask));
        end
    end
end

cmatdim=size(cmat{1},1);
n=length(cmat);
t=reshape(cell2mat(cmat),cmatdim,cmatdim,n);

cmat_mean{1}=mean(t,3); % mean
for i=1:nsites*ndays
    for j=1:nsites*ndays
        cmat_mean{2}(i,j)=std(t(i,j,:)); % stdev
    end
end

figure; image(cmat_mean{1}*50);
axis square
if saveimg
    clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
    timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
    savefig(sprintf('corrsitevday_%s.fig',timeID));
    saveas(gcf,sprintf('corrsitevday_%s.png',timeID))
    close(gcf)
end

