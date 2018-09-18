function [percentchange_vsall,outliers]=calc_loo_percchangeinvar(loo_var,overall_var)
% e.g., lt_percchange=calc_loo_percchangeinvar(lt_var_loo,lt_var,2);
% loo_var should be n var comps x n levels   % matrix=[n var comp x n stuff left out]

if size(overall_var,1)~=size(loo_var,1)
    if size(overall_var,2)==size(loo_var,1) % overall_var is row vector
        overall_var=overall_var';
    else
        error('Dimension mismatch or transpose loo_var mat')
    end
end



meanloovar=mean(loo_var')';

for i=1:size(loo_var,1) % want size==7 bc want to go through ncomp
    percentchange_vsall(i,:)=100*((loo_var(i,:)-overall_var(i))./overall_var(i));
    percentchange_vsloomean(i,:)=100*((loo_var(i,:)-meanloovar(i))./meanloovar(i));
end


%% calculate outliers
outlier_dev1=1*std(percentchange_vsall')'; %1sd
outlier_dev2=2*std(percentchange_vsall')'; %2sd
component_mean=mean(percentchange_vsall')';
ul1=component_mean+outlier_dev1;
ll1=component_mean-outlier_dev1;
ul2=component_mean+outlier_dev2;
ll2=component_mean-outlier_dev2;
for i=1:length(component_mean)
    outliers(i,:)= +(percentchange_vsall(i,:)>ul1(i) | percentchange_vsall(i,:)<ll1(i)) + +(percentchange_vsall(i,:)>ul2(i) | percentchange_vsall(i,:)<ll2(i));
end


%% Plot [-site, subj, day]
% may want to standardize
bar([-percentchange_vsall(1,:)',percentchange_vsall(2,:)',percentchange_vsall(3,:)'])
figure(2)
bar([-percentchange_vsloomean(1,:)',percentchange_vsloomean(2,:)',percentchange_vsloomean(3,:)'])

end
