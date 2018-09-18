function [newdata,stats] = remove_dropout(ref_data,olddata)
% edit 12/22/16: ???
% data = num_regions x num_sessions OR num_sessions of region-containing cells

if(iscell(ref_data(1,1)))
    ref_data = cell2mat(ref_data);
end

if(iscell(olddata(1,1)))
    olddata = cell2mat(olddata);
end


% make mean map
mean_data=transpose(ref_data);
mean_data=mean(mean_data);

avg_roi=mean(mean_data);

% remove rois smaller than 2 stdevs from mean as dropout regions
fractn_to_remove=0.05;
threshold=2.1;
num_rois_excluded=0;
while(num_rois_excluded < (fractn_to_remove*size(mean_data,2)))
    threshold=threshold-0.1;
    excluded_rois=(mean_data<(avg_roi-threshold*std(mean_data)));
    num_rois_excluded=sum(excluded_rois);
end

stats=[num_rois_excluded, threshold, mean(mean_data(excluded_rois)), avg_roi, {excluded_rois}, {mean_data(excluded_rois)}, {find(excluded_rois)}];
newdata=olddata(~excluded_rois,:);