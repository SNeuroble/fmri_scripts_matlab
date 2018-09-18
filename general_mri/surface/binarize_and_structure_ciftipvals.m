% input is pvals derived from cifti, output is nifti
% mask and make data for subject, scanner, and day effects
% p_masked_corrected, mask_all, gleft, and gright (any right and left func.gii file) should be loaded

l=length(mask_all{1});

% takes 1 and 2 are different masking procedures
%% TAKE 1 BEGIN
% 
% p_bin_struct{1}=zeros(l,3);  p_bin_struct{2}=zeros(l,3);
% 
% for m=1:2
%     p_masked_bin{m}(p_masked_corrected{m}>0.05)=1;
%     p_masked_bin{m}(p_masked_corrected{m}<0.05)=2;
%     for i=1:3
%         p_bin_struct{m}(:,i)=structure_data(p_masked_bin{m}(:,i),mask_all{m},'none');
%     end
% end



%% TAKE 2 BEGIN
% this also allows correction

% corrected
% for i=1:3
% [~,p_corrected(:,i)]=mafdr(p(:,i));
% end
% 
% p_sep{1}=p_corrected(1:l,:);  p_sep{2}=p_corrected(l+1:2*l,:);

% non-corrected
p_sep{1}=p(1:l,:);  p_sep{2}=p(l+1:2*l,:);
% corr/noncorr ends

p_bin_struct{1}=zeros(l,3);  p_bin_struct{2}=zeros(l,3);

for m=1:2
    p_bin_struct{m}(p_sep{m}>0.05)=1;
    p_bin_struct{m}(p_sep{m}<0.05)=2;
    for i=1:3
        p_bin_struct{m}(:,i)=p_bin_struct{m}(:,i).*+mask_all{m};
    end
end



%% Saving for Take 1 and 2
gleft.cdata=p_bin_struct{1}(:,1);
save(gleft,'pmask_L_scanner.func.gii')
gright.cdata=p_bin_struct{2}(:,1);
save(gright,'pmask_R_scanner.func.gii')

gleft.cdata=p_bin_struct{1}(:,2);
save(gleft,'pmask_L_scanner.func.gii')
gright.cdata=p_bin_struct{2}(:,2);
save(gright,'pmask_R_scanner.func.gii')

gleft.cdata=p_bin_struct{1}(:,3);
save(gleft,'pmask_L_day.func.gii')
gright.cdata=p_bin_struct{2}(:,3);
save(gright,'pmask_R_day.func.gii')