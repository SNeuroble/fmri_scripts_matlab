% procedures={'pcc_stats','pcc_none_stats','rmc_stats','rmc_none_stats','lt_stats','lt_none_stats','icd_stats','icd_none_stats','matrix_stats','matrix_none_stats'};
% procedures={'tmp_icc_1run','tmp_icc_2runs','tmp_icc_3runs','tmp_icc_4runs','tmp_icc_5runs','tmp_icc_6runs',...
%     'tmp_icc_14H16M','tmp_icc_14H22M','tmp_icc_14H30M','tmp_icc_14H38M','tmp_icc_14H47M','tmp_icc_14H56M'};
%procedures={'tmp_icc_23H5M','tmp_icc_23H6M','tmp_icc_23H7M','tmp_icc_23H8M','tmp_icc_23H9M','tmp_icc_23H10M'};
procedures={'tmp_icc_09H31M48S','tmp_icc_09H40M33S','tmp_icc_09H49M14S','tmp_icc_09H57M42S','tmp_icc_10H06M10S','tmp_icc_10H14M06S'};
% unthresh then thresh

% load('ftbl_new.mat')
ftbl=ftbl_new;
n_prime=1;
n_inc=1;

it=1;
for procedure=procedures
	disp(sprintf('Doing %s',procedure{1}))
	load(sprintf('%s.mat',procedure{1}),'stats');
    for i=1:3
        [fps, qs]=mafdr(getvalsfromcellarray(stats,i+1,4,1));
        qsummary(it,i)=sum(qs<0.05);
    end
% 	[icc,var]=stats_to_icc(stats,ftbl,n_prime,n_inc);
it=it+1;
end

co = get(gca,'ColorOrder');
set(gca, 'ColorOrder', [1 0 0; 0 1 0; 0 0 1], 'NextPlot', 'replacechildren');
plot(qsummary/length(stats)*100)

%{
% load data and exclude subs
for procedure=procedures
	disp(sprintf('Doing %s, excl sub 4',procedure{1}))
	[data,ftbl] = load_reliability_data(sprintf('/Users/stephanie/Documents/data/traveling_subs/results_%s/',procedure{1}),'.*\.nii\.gz');
	data_excl=data(ftbl(:,1)~=4);
	ftbl_excl=ftbl(ftbl(:,1)~=4,:);
	[icc,var,stats,sigmask]=run_reliability('none',data_excl,ftbl_excl);

	stats_bonf=stats(sigmask);
	[icc,var]=stats_to_icc(stats_bonf,ftbl_excl,5,5,.25);
end
%}
