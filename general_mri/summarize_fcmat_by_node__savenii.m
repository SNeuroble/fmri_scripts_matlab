function summarize_fcmat_by_node__savenii(matrix)
% input is nxn connectivity matrix

if size(matrix,1)~=size(matrix,2)
    matrix=structure_data(matrix);
end

nnodes=size(matrix,1);

reference_img=load_untouch_nii('~/bioimagesuite35/images/shenetal_neuroimage2013/shen_1mm_268_parcellation.nii.gz');
output_mat=zeros(size(reference_img.img));

for i=1:nnodes
    output_mat(reference_img.img==i)=mean(matrix(:,i));
end

reference_img.img=output_mat;
% reference_img.hdr.dime.datatype=64;
% reference_img.hdr.dime.bitpix=64;
save_untouch_nii(reference_img,'tmp.nii.gz');
