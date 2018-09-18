function [roi_iccs,roi_var,roi_stats] = run_trav_loo(data,ftbl)
% data from load_trav_data
% may want to turn saving off in 'run_reliability'

this_factor=2; % site
n_levels=max(ftbl(:,this_factor));

    for this_level=1:n_levels
        these_ids=find(ftbl(:,this_factor)~=this_level);
        ftbl2=ftbl(these_ids,:);
        data2=data(these_ids);

    %     [roi_iccs{this_level},roi_var{this_level},roi_stats{this_level}]=calc_roi_iccs(data2,ftbl2,'all');
        [roi_iccs{this_level},roi_var{this_level},roi_stats{this_level}]=run_reliability('none',data2,ftbl2);

        
    end


end

