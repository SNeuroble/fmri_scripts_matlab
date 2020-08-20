function [roi_iccs,icc_mean]=calc_roi_iccs_excl_dropout(ref_data,processed_data,categories)
% data must be cell arrays; ref data is used for dropout calculation,
% processed_data for analysis
% hint: use mean functional data for ref_data

[new_processed_data,stats] = remove_dropout(ref_data,processed_data);
[roi_iccs,icc_mean,roi_stats] = calc_roi_iccs(new_processed_data,categories);