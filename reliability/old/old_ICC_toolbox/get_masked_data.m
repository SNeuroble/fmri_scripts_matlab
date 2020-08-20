function [masked_data]=get_masked_data(data,mask)
% data is cell array 1xnsessions; mask is just one matrix

mask=logical(mask); % NEW

data_new=[];
for(i=1:length(data))
    data_new=data{i}(mask);
    masked_data{i}=data_new;
end