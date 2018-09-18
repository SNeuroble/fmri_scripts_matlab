function [dataavg] = averagecellmatrices(data,ftbl,allx,ally,allz,dosubnetworks)
% x, y, and z specify the number of levels for each factor (x,y,z) to average/iterate over
% ftbl must be reordered to match x,y,z
% IMPORTANT DISTINCTION:
%       if want to average over a list (e.g., avg over runs 1:4), provide allz=[1:4]';
%       if want to iterate over a list (e.g., matrix for each of runs 1:4), provide allz=1:4;
% Remember: returns lower triangle mask
% e.g., averagecellmatrices(data,[ftbl(:,1) ftbl(:,4), ftbl(:,5)],1:nsubs,[1:nruns]',1:nsess)

% options:  iterate over each sub/sess/day
%           use fixed nsub/sess/day as limit (if single is provided)

if ~exist('dosubnetworks','var') | isempty(dosubnetworks)
    dosubnetworks = 0;
end

if ~iscell(data); % shouldn't really provide non-cell non-matrices, though...
    for i=1:length(data); tmp(i)={data(i)}; end
    data=tmp;
end

matdim=size(data{1});

if (matdim(1) == matdim(2)) && (matdim(1)>1)
    trimask=logical(tril(ones(size(data{1})),-1));
end

datait=1;

for x=allx
    for y=ally
        for z=allz
            
%             if length(allx)==1; x=1:allx; end
%             if length(ally)==1; y=1:ally; end
%             if length(allz)==1; z=1:allz; end
            
            it=ismember(ftbl(:,1),x) & ismember(ftbl(:,2),y) & ismember(ftbl(:,3),z);
            tmp=data(it);
            
            if iscell(tmp)
                tmp=reshape(cell2mat(tmp),[matdim(1),matdim(2),sum(it)]);
            end
            
            if ~(length(x)==1 && length(y)==1 && length(z)==1)
                tmp=mean(tmp,ndims(tmp));
            end
            
            if dosubnetworks
                tmp=domatrixsummary_avg(tmp,'subnetwork',0,1,0);
            end
            
            if (matdim(1) == matdim(2)) && (matdim(1)>1)
                tmp=tmp(trimask);
            else
                tmp=tmp(:);
            end
            
            
            dataavg(:,datait)=tmp;
            datait=datait+1;
            
        end
    end
end


