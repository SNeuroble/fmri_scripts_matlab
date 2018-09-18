function compare_model_fits_Ftest(stats1,stats2)

% Example:
% compare whether node size and location relate to reliability
% (load nodesize v. icc fig) or load nodesize and icc
% getdatafromfig
% (reorder x and y data so that first 208 are ctx, second are ctx, e.g., x_reord=reorder_matrix_atlas(xdata','subnetwork');)
% location=[ones(208,1); zeros(60,1)]; % ones correspond with ctx
% (IMPORTANT! zscore both predictors and outcome to get standardized beta coefficients)
% stats1=fitlm(zscore(location), zscore(icc));
% stats2=fitlm([zscore(nodesize), zscore(location)], zscore(icc)); % analogous to icc ~ nodesize + location
% compare fits

% if higher DF in stats 2, must switch 1 and 2 (1 should have higher DF and higher SSE)
if (stats1.DFE < stats2.DFE)
    stats=stats1;
    stats1=stats2;
    stats2=stats;
end

% get F-stat
F=((stats1.SSE-stats2.SSE)/(stats1.DFE-stats2.DFE))/...
    (stats2.SSE/stats2.DFE);

% get p
p=1-fcdf(F,stats1.DFE,stats2.DFE);

fprintf('F=%2.4f (p=%2.4E)\n',F,p)
