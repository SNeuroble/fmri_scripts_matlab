function [y_predict,residuals,SSE]=fit_spline(x,y,smoothing,plot_filename)  
%% Fit function
% fit cubic spline to sliding window avg with 50% overlap
% also tried sigmoid - derived from https://stackoverflow.com/questions/16509696/fitting-sigmoid-to-data/16511754
% sigfunc = @(A, x)(A(1)+(A(2)-A(1)) ./ (1 + exp(-A(3)*x+A(4))/A(5)));
% A0 = ones(1,5);
% A_fit = nlinfit(x, y, sigfunc, A0);
% y_predict=sigfunc(A_fit,x);

% sliding window to get smooth model
window_sz=0.1;
percent_overlap=0.5; % percent overlap (0.5=50%)
overlap=window_sz*percent_overlap;
nwindows=ceil( (max(x) - min(x) - overlap) / (window_sz - overlap) );
for i=1:nwindows
    x_windowed(i)= min(x)+i*(window_sz-overlap);
    y_predict_windowed(i)=mean(y(x>(x_windowed(i)-overlap) & x<=(x_windowed(i)+overlap)));
end

% fit spline to smooth model
y_predict = csaps(double(x_windowed),double(y_predict_windowed),smoothing,x);

% descriptive stuff
residuals=y-y_predict;
SSE=sum(residuals.^2);

% OR:
% y_predict = smooth(x,y,0.5,'rloess');

