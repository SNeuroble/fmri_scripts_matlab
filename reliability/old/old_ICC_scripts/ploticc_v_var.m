function ploticc_v_var(stats,icc)
% makeplot(stats,icc{2,4}{1})
% icc_ctx=getctx(icc{2,4}{1},'cortex','subnetwork');
% quick script

% icc - ctx, subctx
icc_ctx=getctx(icc{2,4}{1},'cortex','subnetwork');
icc_subctx=getctx(icc{2,4}{1},'subcortex','subnetwork');

% var - ctx, subctx
nvarcomps=size(stats{1},1)-3;
for i=1:nvarcomps
    var{i}=getvalsfromcellarray(stats,i+1,5,1);
end

var=cell2mat(var);
var(:,nvarcomps+1)=sum(var,2);


nvar2plot=nvarcomps+1;

for i=1:nvar2plot
    
    v_ctx=getctx(var(:,i),'cortex','subnetwork');
    v_subctx=getctx(var(:,i),'subcortex','subnetwork');
    
    subplot(2,ceil(nvar2plot/2),i)
    scatter(v_ctx(:),icc_ctx(:),'.')
    xlim([-0.01 0.2])
    hold on
    scatter(v_subctx(:),icc_subctx(:),'.')
    hold off
    
end


