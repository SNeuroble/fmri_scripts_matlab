function [iccs_in_categories]=categorize_iccs(data)
% input = row or column

% criteria via Cicchetti and Sparrow (1981)
poor=[0,0.4];
fair=[0.4,0.59];
good=[0.6,0.74];
excellent=[0.74,1];

poor_iccs=find(data>=poor(1) & data<poor(2));
fair_iccs=find(data>=fair(1) & data<fair(2));
good_iccs=find(data>=good(1) & data<good(2));
excellent_iccs=find(data>=excellent(1) & data<excellent(2));

iccs_in_categories=[{poor_iccs},{fair_iccs},{good_iccs},{excellent_iccs}];

% % can also section by stdev
% top=data>(mean(data)+threshold*std(data)); % binary
% bottom=data<(mean(data)+threshold*std(data));
% top=find(top);
% bottom=find(bottom);