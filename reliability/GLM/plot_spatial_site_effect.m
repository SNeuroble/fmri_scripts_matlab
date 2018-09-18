function plot_spatial_site_effect(effect_summary)
% effect_mat is         {2x6}               {1x9}   {1x14}
% corresponding with    {pvals(label+data)} {sites} {anatomical_groups}

n_pthresh=size(effect_summary,2)-1; % remove titles
n_sites=length(effect_summary{2,1})-1; % remove intercept
n_anat_groups=length(effect_summary{2,1}{1});

pthresh_to_plot=[1,3,5,6]; % 

it=1;
figure()
% for this_pthresh=1:n_pthresh
for this_pthresh=pthresh_to_plot
    effect_mat{this_pthresh}=cell2mat(effect_summary{2,this_pthresh});
    effect_mat{this_pthresh}=reshape(effect_mat{this_pthresh},n_anat_groups,n_sites+1); %+1 to include intercept
    effect_mat{this_pthresh}=effect_mat{this_pthresh}(1:end,2:end);
    subplot(length(pthresh_to_plot),1,it)
    bar(effect_mat{this_pthresh})
    title(sprintf('p<%s',effect_summary{1,this_pthresh}))
    it=it+1;
end

% xlabel('Sites (S1, S2, S3, G4, S5, G6, G7, S8)')
[ax1,h1]=suplabel('% edges in anatomical group','y');
[ax2,h2]=suplabel('1-R frontal, 2-R insula, 3-R parietal, 4-R temporal, 5-R occipital, 6-R limbic, 7-R subcort.+cbl, 8-L frontal, 9-L insula, 10-L parietal, 11-L temporal, 12-L occipital, 13-L limbic, 14-L subcort.+cbl)');

set(h1,'FontSize',14)



