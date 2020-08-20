function mask = create_sigedge_mask(data,p_thresh,correctiontype)
% data should be raw data via load_trav_files or gm-masked data
% suggestion: p_thresh=0.05

if strcmp(correctiontype,'none')
	mask = ones(size(data{1}));	
	return;
end


data=rearrange_trav_cell_matrix(data);

m=length(data); % num units of brain (bc rearranged)

% Bonferroni
if strcmp(correctiontype,'Bonf')
for(i=1:m)
H(i)=ttest(data{1,i},0,'alpha',p_thresh/m);
% returns separate H (0 or 1) for each col in data{1,i}, e.g., 278 for mat

end

% FDR
% TODO: just assign H(i) instead of H{i}, same w p, same w for bonferroni

elseif strcmp(correctiontype,'FDR')
for(i=1:m)
[H(i),ps(i)]=ttest(data{1,i},0,'alpha',p_thresh);
% returns separate H (0 or 1) for each col in data{1,i}, e.g., 278 for mat
end
fdr_ps=mafdr(ps);
H2=fdr_ps<0.05;
H=H2;
end

% H=cell2mat(H);
mask=reshape(H,m,1);
mask(isnan(mask))=0;
mask=logical(mask);

