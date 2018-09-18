% simplescatterplotter_formyelincomparison

% load('comparison_icc_myelin.mat')

% data=[icc_hcp(~all_nans)',myelin_hcp(~all_nans)',icc_mrtool(~all_nans)',myelin_mrtool(~all_nans)'];
% datapair(1,:)=[1,2];
% datapair(2,:)=[3,4];
% datapair(3,:)=[2,4];
% datapair(4,:)=[1,3];

% data=[icc(~all_nans),my(~all_nans)];
data=[icc_parc,my_parc];
datapair(1,:)=[1,2];


for i=1:size(datapair,1)
    a=figure;
    mysc=line(data(:,datapair(i,1)),data(:,datapair(i,2)));
    set(mysc                         , ...
        'LineStyle'       , 'none'      , ...
        'Marker'          , 'o'         , ...
        'MarkerSize'      , 7           , ...
        'MarkerEdgeColor' , 'none'      , ...
        'MarkerFaceColor' , [.7 .7 1] );
    
    %         'MarkerSize'      , 7           , ...
    
    myl=lsline;
    set(myl                            , ...
        'Color'       , [.2 .7 .8]      , ...
        'LineWidth'       , 1.3         );
    saveas(a,sprintf('img%d.png',i))
end

close all
