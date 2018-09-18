function [ordered_mat]=summarize_columns(data)
% data is already filtered for p<threshold

%typical process:
% for (i=1:9) ps0001_reord{i}=double(ps{i}<0.0001).*ps{i}; end
% ps0001_reord_summary=summarize_columns(ps0001_reord);
% ps0001_toplot=reshape(cell2mat(ps0001_reord_summary),14,9);
% ps0001_toplot(:,1)=[];
% bar(ps0001_toplot)

ordering=[1,41,44,64,84,100,113,139,167,181,202,217,234,249,279];
% right: frontal-1, insula-2, parietal-3, temporal-4, occipital-5, limbic-6, subcortical+cbl-7

%ordering=[1,30,41,44,64,84,100,113,129,136,139,167,179,181,202,217,234,249,266,274,279];
% right: 1=right prefrontal, 2 = right sensorimotor, 3=right insula, 4=right par, 5=r temporal, 6 = r occipital, 7=limbic, 8=cbl, 9=thal, 10=brainstem

for(i=1:length(data))
    
    cols=sum(data{i}>0);
    
    for(j=1:(length(ordering)-1))
        
        ordered_mat{i}(j)=sum(cols(ordering(j):(ordering(j+1)-1))) / ((ordering(j+1)-1) - ordering(j));
        % these are NORMALIZED by the lobe size
        
    end
    
end