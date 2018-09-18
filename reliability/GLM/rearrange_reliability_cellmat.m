function rearranged_data = rearrange_reliability_cellmat(data)
% seed or edge vector: nsessions vector of 1xnvars -> nvars vector of nsessionsx1

nsess=length(data);
nvars=length(data{1}); % also, nvars=size(d2,1);
d2=cell2mat(data);
d2=d2';
rearranged_data=mat2cell(d2,nsess,ones(1,nvars));

% full matrix: nsessions vector of nROIxnROI -> nROI vector of nsessionsxnROI matrices
% IN PROGRESS: if there's a full matrix
% (usually shouldn't be bc symmetric, so should have done trimask)
% nroi1=size(data{1},1);
% nroi2=size(data{1},2);
% nsess=length(data);
% d2=cell2mat(data);
% d3=mat2cell(d2,ones(1,nroi1));
% % do some stuff to reshape matrix within each cell
% d3=d3';