function plotremovals(x,thisy,theselabels,thisaxis,dovertline,vertline_xs,dohorzline,horzline_ys)

% one plot
% x=1:26;
% y=[fliplr(summary_decrease.r{i}') summary.r{i}'];
% theselabels=[fliplr(nedges_decrease)'; nedges']; % numeric
% thisaxis=[0 27 -0.1 0.3];
% x_line=13.5;

% better plot (from inv ICC edges removed, remove (1) mirrored 35778 edges; remove (2) single edge removed)
% x=1:24;
% y=[fliplr(summary_decrease.r{i}(3:end)') summary.r{i}'];
% theselabels=[fliplr(nedges_decrease(3:end))'; nedges']; % numeric
% thisaxis=[0 25 -0.1 0.25];
% x_line=12;

theselabels=num2str(theselabels);
markersz=70;
markerfilled='filled';
lw=1;

figure
scatter(x,thisy,markersz,markerfilled)
set(gca,'xtick',x)
set(gca,'XTickLabel',theselabels)
set(gca,'XTickLabelRotation',45)
axis(thisaxis)

if dovertline
    hold on
    for thisx=vertline_xs
        plot([thisx thisx],[thisaxis(3) thisaxis(4)],'r:','LineWidth',lw)
    end
    hold off
end

if dohorzline
    hold on
    for thisy=horzline_ys
        plot([thisaxis(1) thisaxis(2)],[thisy thisy],'r:','LineWidth',lw)
    end
    hold off
end


% show x in mid
hold on
xmid=floor(length(x)/2);
plot([xmid xmid],[thisaxis(3) thisaxis(4)],'k-','LineWidth',0.1)
hold off

% show y=0
hold on
plot([thisaxis(1) thisaxis(2)],[0 0],'k-','LineWidth',0.1)
hold off

