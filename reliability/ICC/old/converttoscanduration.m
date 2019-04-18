function y=converttoscanduration(x,originterval,desiredinterval)
% convert x from units to scan duration
% e.g.,
% load('iccupto6','icc')
% icc2=cell2mat(icc{2,6});
% icc3=sum((+(icc2<0.4))')';
% icc3=icc3+1; % careful here...
% icc_scdur=converttoscanduration(icc3,[1,31],[6,24]);

minunits=originterval(1);
maxunits=originterval(2);
minscandur=desiredinterval(1);
maxscandur=desiredinterval(2);

m=(maxscandur-minscandur)/(maxunits-minunits);
b=maxscandur-m*(maxunits);
y=m*x+b;

