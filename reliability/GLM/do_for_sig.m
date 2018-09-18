function masked_mean=do_for_sig(data,fdr_ps,this_fun)
% data is an 1xn cell
% get mean, magnitude, or variance of sig values (e.g., bs) for each group
% this_fun: 'mean' - meanval,'magn' - mean absval, or 'vari' - variance

data_vec=zeros(size(data{1},1)*size(data{1},2),length(data)-1);
for i=2:length(data)
    data_vec(:,i-1)=data{i}(:);
end

pmask_vec=fdr_ps<0.05;

masked_mean=zeros(1,size(data_vec,2));
switch this_fun
    case 'mean'
        for i=1:size(data_vec,2)
            masked_mean(i)=mean(data_vec(pmask_vec(:,i)));
        end
    case 'magn'
        for i=1:size(data_vec,2)
            masked_mean(i)=mean(abs(data_vec(pmask_vec(:,i))));
        end
    case 'vari'
        for i=1:size(data_vec,2)
            masked_mean(i)=var(data_vec(pmask_vec(:,i)));
        end
    otherwise
        error('Please enter a valid function.')
        
end