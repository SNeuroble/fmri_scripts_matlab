function ploticcsurface(icc)
% does only for Gcoeff

inc=0.25;
init=1-inc;
xlimit=size(icc{2,1},1)-1; %assuming square

icc_thresh=(icc{2,1}>0.4)+(icc{2,1}>0.6)+(icc{2,1}>0.74);
icc_transitions=icc_thresh(2:17,2:17,2:17)-icc_thresh(1:xlimit,1:xlimit,1:xlimit); 

icc_thresh_smaller=icc_thresh(2:xlimit+1,2:xlimit+1,2:xlimit+1);

icc_thresh2=icc_transitions;
icc_thresh2(logical(icc_transitions))=icc_thresh_smaller(logical(icc_transitions));
icc_thresh2(icc_thresh2>1)=0;

Z=zeros(xlimit,xlimit);
for n=1:xlimit
    k=find(squeeze(icc_thresh2(n,:,:)));      
    Z(k)=n*inc+init;
end

n_end=n*inc+init;

Z(Z==0)=init+inc;
[X,Y] = meshgrid(1:inc:n_end,1:inc:n_end);
hSurface = surf(X,Y,Z);

% plot
view(135,45)
colormap(summer(100))
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha',0.25); % single color
ax = gca;
ax.GridAlpha=0.3;

end

