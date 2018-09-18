%procedure
 
% normal
load('sandbox_cmat.mat')
[r,p]=imgcorr_v_interscaninterval(cmat,1,0);


%tril data
trilmask=logical(tril(ones(268),-1));
%{
for i=1:length(data_new)
    data2{i}=data_new{i}(trilmask);
end
 
%}
%mask defined as in isi script
mask=zeros(48);
for i=1:12
    mask((i-1)*4+1:i*4,(i-1)*4+1:i*4)=1;
end
mask=logical(tril(mask,-1));
ids=find(mask);
[ids_sub(:,1),ids_sub(:,2)]=ind2sub(size(mask),ids);
 %{
for i=1:size(ids_sub,1)
    mini{i}=minicorrs(data2{ids_sub(i,1)},data2{ids_sub(i,2)});
end
%}

load('~/Documents/data/mnt/minicorrs_isi.mat')
 
nedges=size(mini{1},1);
thiscmat=double(mask);
minicorr_w_isi=zeros(nedges,1);
minicorr_w_isi__sig=zeros(nedges,1);
% get vals
for i=1:nedges
    vals=getvalsfromcellarray(mini,i,1,1);
    for j=1:size(ids_sub,1)
        thiscmat(ids_sub(j,1),ids_sub(j,2))=vals(j);
    end
    [minicorr_w_isi(i),minicorr_w_isi__sig(i)]=imgcorr_v_interscaninterval(thiscmat,0,0);
end
domatrixsummary_avg(minicorr_w_isi*10,'reorderimg',1);
drawmatrix_atlas(minicorr_w_isi,'reorderimg',1)


% same, except start from mean of edges
%{
for i=1:length(mini)
    miniavg(:,:,i)=domatrixsummary_avg(mini{i},'suppressimg',1);
end
 
newids=find(~isnan(miniavg(:,:,1)));
[cmatids_sub(:,1),cmatids_sub(:,2)]=ind2sub(size(miniavg(:,:,1)),newids);

miniavgcorr_w_isi=zeros(size(miniavg(:,:,1)));
miniavgcorr_w_isi__sig=zeros(size(miniavg(:,:,1)));
thiscmat=double(mask);
for i=1:length(newids)
    vals=squeeze(miniavg(cmatids_sub(i,1),cmatids_sub(i,2),:));
    for j=1:size(ids_sub,1)
        thiscmat(ids_sub(j,1),ids_sub(j,2))=vals(j);
    end
    [miniavgcorr_w_isi(cmatids_sub(i,1),cmatids_sub(i,2)),miniavgcorr_w_isi__sig(cmatids_sub(i,1),cmatids_sub(i,2))]=imgcorr_v_interscaninterval(thiscmat,0,0);
end

% prob crap
t=reshape(miniavg(2,1,:),6,12);
t=mean(t,2);
order=[1,4,6,2,5,3]; % 1~=4~=6<2~=5<3
scatter(1:6,t(order));
corr([1:6]',t(order))
%}