function saveimg_squaresubplot(filenames,loadsavedimg,m,n)
% e.g.,filenames={'tmp_DStudy_21H52M22S','tmp_DStudy_22H01M15S','tmp_DStudy_22H10M09S','tmp_DStudy_22H18M32S','tmp_DStudy_22H27M23S','tmp_DStudy_22H35M37S',...
% 'tmp_DStudy_22H37M07S','tmp_DStudy_22H39M27S','tmp_DStudy_22H42M17S','tmp_DStudy_22H45M29S','tmp_DStudy_22H48M48S','tmp_DStudy_22H52M06S'};
%
% saveimg_squaresubplot(filenames,'no')

if strcmp(loadsavedimg,'no')
    loadsavedimg=0;
elseif strcmp(loadsavedimg,'yes')
    loadsavedimg=1;
else
    loadsavedimg = questdlg('Load from a saved image?');
    if strcmp(loadsavedimg,'no')
        loadsavedimg=0;
    elseif strcmp(loadsavedimg,'yes')
        loadsavedimg=1;
    else
        return;
    end
end


for filename=filenames
    
    filename=filename{1};
    
    if loadsavedimg
        figure(1)=openfig(sprintf('%s.fig',filename));
    end
    
    if ~exist('n','var')
        [m,n]=size(findall(gcf,'type','axes'));
    end
    for i=1:(m*n)
        subplot(m,n,i);
        axis('square')
    end
    
    saveas(gcf,sprintf('%s.png',filename))
    close(gcf)
end