function [edge_corrs,mean_corr] = calc_edge_corrs(roi_data,varargin)
%% read inputs
% roi_data and mask assumed to be cell arrays;

% define inputs
if(size(varargin)>0)
    if(strcmp(varargin{1},'mask'))
        usemask=1;
        mask=varargin{2};
        if (size(mask,1)==1)
            mask=mask{1}(:,:);
        end
        tct=3;
    else
        usemask=0;
        tct=1;
    end
    
    if(size(varargin)==tct)
        use_two_col_tbl = varargin{tct}; % use_two_col_tbl if imported table data that has multiple columns
    else
        use_two_col_tbl = 0;
    end
else
    use_two_col_tbl = 0;
    usemask=0;
end


 


%% set up
% assuming rois order is: [scanner1_subj1_day1, scanner1_subj1_day2, scanner2_subj1_day1, &c]

if(use_two_col_tbl)
    data_column=2;
else
    data_column=1;
end



excluded_sessions=[]; % note: dimension mismatch in timeseries for [80:91]


reference_session=1;
reference_edges=roi_data{reference_session}(:,data_column);



if(use_two_col_tbl)
    reference_edges=table2array(reference_edges);
end


%edges = [1:size(roi_data{1}(:,1),1)];

ngroups=3;
sessions{1} = [1:size(roi_data,2)];
sessions{2} = [1 2 15 28 29 42 56 57 70 80 92 93];
sessions{3} = [2]; % GE Scanners = [42 70 80]

x(1:ngroups) = 1;


edge_corrs = [];


%% compute corr between sessions for all individuals


for (group=[1:ngroups])
    
    for (i=sessions{group})
        
        if (i~=reference_session && ~any(i==excluded_sessions))

            if(use_two_col_tbl)
                thisroi_edges=table2array(thisroi_edges);
            else
                thisroi_edges = roi_data{i}(:,data_column);
            end

            if(~usemask)
                edge_corrs(x(group),group)=corr(reference_edges,thisroi_edges);
            else
                edge_corrs(x(group),group)=corr(reference_edges(mask>0),thisroi_edges(mask>0));
 
            end
            
            
            x(group)=x(group)+1;
        end
        
    end
    
end


%% Draw hist
nbins=25;
cmap=hsv(8);

[series(1,:),centers] = hist(edge_corrs(1:(x(1)-1),1),nbins);

for (group=[2:(ngroups)])
    [series(group,:)] = hist(edge_corrs(1:(x(group)-1),group),centers);
end

DataSum=series(1,:); % sum(series) - ONLY USE IF using non-overlapping data
figure
width(1) = 0.5;
bar(centers,DataSum,width(1),'FaceColor',[0.2,0.2,0.5],....
                     'EdgeColor','none');

                 
for (group=[2:(ngroups)])              
    hold on
    width(group) = width(1);
    bar(centers,series(group,:),width(group),'FaceColor',cmap(group,:),...
                     'EdgeColor','none');
end


legend(num2str(transpose(1:ngroups))) % add legend
hold off




mean_corr=mean(edge_corrs(:,1));

