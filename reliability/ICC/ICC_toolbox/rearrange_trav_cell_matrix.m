function rearranged_data = rearrange_trav_cell_matrix(data)
% session vector of ROIxROI matrices -> ROI vector of ROIxsession matrices

[rearranged_data{1:size(data{1},1)}] = deal(zeros(size(data,2),size(data{1},2)));

for seed=1:size(data,2) %roi seeds
    for edge=1:size(data{1},2) %roi edges for seed
        for site=1:size(data{1},1) %sites
            rearranged_data{site}(seed,edge)=data{seed}(site,edge);
        end
    end
end

