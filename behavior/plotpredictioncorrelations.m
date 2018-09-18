function plotpredictioncorrelations(detailedprediction,nruns,varargin)
% specify number of runs,e.g., plotpredictioncorrelations(detailed,6)
% can also pass r and p directly, e.g., plotpredictioncorrelations([r,p],6)
% can specify ICC to compare with, e.g., plotpredictioncorrelations(detailed,6,icc{2,4}{1})
% detailed is one of the outputs of the doprediction_variablenruns script, run WITHOUT loo

doreliabilityplot=0;

if (~isempty(varargin))
    reliability=varargin{1};
    doreliabilityplot=1;
end

if istable(detailedprediction) | isstruct(detailedprediction)
    r=detailedprediction.r{nruns}{1};
    p=detailedprediction.p{nruns}{1};
    if iscell(r)
        r=r{1};
        p=p{1};
    end
else
    r=detailedprediction(:,1); p=detailedprediction(:,2);
end

caxislimits=[-0.1 0.1];
% colorbarticks=[caxislimits(1):0.01:caxislimits(2)];
% colorbarlabels=cellstr(num2str(colorbarticks'))';


clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});

drawmatrix_atlas(r,'reorderimg',1,'atlascategory','subnetwork','datacategory','none')
colormap(bipolar)
caxis(caxislimits)
% colorbar('Ticks',colorbarticks,'TickLabels',colorbarlabels);
saveas(gcf,sprintf('corrplotdetailed_%s.png',timeID))
close;

domatrixsummary_avg(r,'reorderimg',1,'atlascategory','subnetwork','datacategory','none');
colormap(bipolar)
caxis(caxislimits)
% colorbar('Ticks',colorbarticks,'TickLabels',colorbarlabels);
saveas(gcf,sprintf('corrplotsummary_%s.png',timeID))
close;



% versus reliability
if doreliabilityplot
    % linear
    figure;
    h=plot(reliability,abs(r));
    
    c=colorcube; thiscolor=c([2],:);
    thesemarkers={'.'};
    markersz=22;
    set(h,'Color',thiscolor,'LineStyle','none','Marker','.')
    
    h=lsline;
    set(h,'Color','k')
    
    saveas(gcf,sprintf('behcorrvrel_%s.png',timeID))
    close;
    [r_behcorr_rel,p_behcorr_rel]=corr(reliability,abs(r));
    sprintf('Correlation: r=%0.5f,p=%0.5f',r_behcorr_rel,p_behcorr_rel)
    
    % boxplot
    figure;
    grouping=p<0.05;
    boxplot(reliability,grouping);
    saveas(gcf,sprintf('boxplotsigvrel_%s.png',timeID))
    close;
    [~,p_behcorr_rel]=ttest2(reliability,grouping);
    delmu=mean(reliability(grouping))-mean(reliability(~grouping));
    esz=delmu/std_pooled(reliability(grouping),reliability(~grouping));
    sprintf('t-test: effect size=%0.5f, p=%0.5f',esz,p_behcorr_rel)
    
end



