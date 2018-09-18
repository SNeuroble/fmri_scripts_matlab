% thanks to https://www.mathworks.com/matlabcentral/answers/100687-how-do-i-extract-data-from-matlab-figures?requestedDomain=www.mathworks.com
% Designed for 2D scatterplot. Nothing else can be plotted on the fig.

h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

%% If multiple subplots:
%
% dataObjs lists subplots backwards and
% lists menu first. So, for 7 subplots,
% dataObjs{1} is menu, dataObjs{2} is subplot7, and dataObjs{3} is subplot6,&c...
% Therefore, to get subplot 2 out of 7:

% xdata = get(dataObjs{7}, 'XData');
% ydata = get(dataObjs{7}, 'YData');

% And overlapping layers are listed backwards.
% Therefore, to get the first layer of two layers:

% xdata=flipud(xdata);
% ydata=flipud(ydata);

