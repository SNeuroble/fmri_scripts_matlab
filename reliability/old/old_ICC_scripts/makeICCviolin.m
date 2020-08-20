function makeICCviolin(data)

dooutline=0;
ho=1;
dw=0.8;
mm=5;
col=parula;
col=flipud(col);

% e.g., data=icc_pccsig{2,5}{1:4:9};
figure;
if dooutline
    distributionPlot(data,'histOpt',ho,'showMM',0,'color','black','distWidth',0.9);
    hold on
    distributionPlot(data*0.97,'histOpt',ho,'showMM',mm,'color','w','distWidth',dw);
else
    distributionPlot(data,'histOpt',ho,'showMM',mm,'distWidth',dw,'color','');
    hold on
end

% axis + thresholds
axis([0 4 0 1])
plot(linspace(0,4,30),ones(1,30)*0.4,'.','color','b')
plot(linspace(0,4,30),ones(1,30)*0.6,'.','color','g')
plot(linspace(0,4,30),ones(1,30)*0.74,'.','color','y')

hold off