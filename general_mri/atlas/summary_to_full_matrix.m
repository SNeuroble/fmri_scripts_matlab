function fullmat=summary_to_full_matrix(summarymat)
%% TODO: update with map
% if have separate categorized matrices (1,2,3,4) and want to convert to single matrix
% to make categorized:
% for i=1:4
%     summat_separatecat{i}=zeros(size(summarymat));
%     summat_separatecat{i}(summarymat==(i+1))=1; % poor, then fair, ... excellent
% end
% then pass in summat_separatecat as summarymat

lobe_mapping = [0,30,41,44,64,84,100,113,129,136,139,167,179,181,202,217,234,249,266,274,278];
lobe_mapping=lobe_mapping+1; % C++ -> Matlab indexing

for i=1:(length(lobe_mapping)-1)
    nodes(i)=lobe_mapping(i)+ceil((lobe_mapping(i+1)-lobe_mapping(i))/2);
end

if ~iscell(summarymat)
    summarymat={summarymat};
end

for k=1:length(summarymat)
    fullmat{k}=zeros(lobe_mapping(end)-1,lobe_mapping(end)-1);
    for i=1:length(nodes)
        for j=1:length(nodes)
            fullmat{k}(nodes(i),nodes(j))=summarymat{k}(i,j);
            
            % plot autocorrelation
            if nodes(i)==nodes(j)
% TODO: buggy snippet; no results
%                diagnode1=lobe_mapping(i)+ceil((1/3)*(lobe_mapping(i+1)-lobe_mapping(i)));
%                diagnode2=lobe_mapping(i)+ceil((2/3)*(lobe_mapping(i+1)-lobe_mapping(i)));

               diagnode2=lobe_mapping(i)+ceil((1/3)*(lobe_mapping(i+1)-lobe_mapping(i)));
               diagnode1=lobe_mapping(i)+ceil((2/3)*(lobe_mapping(i+1)-lobe_mapping(i)));
               
% TODO: buggy snippet here too, BSM connects to PFC
%                diagnode1=nodes(i)+1;
%                diagnode2=nodes(i)-1;

               
               fullmat{k}(diagnode1,diagnode2)=summarymat{k}(i,j);
               
            end
            
        end
    end
    % reflect over diag (summary is usually lower triang)
    fullmat{k}=tril(fullmat{k})+tril(fullmat{k},-1)';
end




end
