function [summary,detailed]=predictbehavior_variabledata(data,ftbl,behavior,mask)
% doprediction for different numbers of runs and masked data; plot corrs and mse
% load('datatrt268') % data=data_new, ftbl=ftbl_new
%     data MUST be cell of 2-D matrices
% load('behaviorpackage') % behavior=stai_nansfilled.y2
% load('iccupto6','icc') % mask=icc{2,4}{1};
% ************* USER DEFINED ************

doloo=1; % also does overall plot
do_loo_plot=1;
bysessions=1;
doplot=0;
flipfigure=0;
pthresh_default=0.05; % reference p 0.05
doremovals=1;
varyps=0;
edgeremovalinterval=3200; % # edges to remove at a time
plotfitline=0;


nruns_range=[6]; % nruns=[1,3,6];
nruns_range_length=length(nruns_range);

ps_range=[0.05,0.005,0.0005]; % for varying ps

% for varying mask
nedges=size(data{1},1)*(size(data{1},1)-1)/2;
mask2{1}=logical(tril(ones(size(data{1})),-1));

if doremovals % remove certain edges
    if ~exist('mask','var')
      error('specified removals but no mask.')
    end
    
    % replace mask by rank
    [~,~,mask]=unique(mask); % third arg is rankings
    mask=structure_data(mask);
    mask=tril(mask,-1);
    
    % count number iterations
    nedges_2ndit=sum(sum(+(mask>1))); % edges left after removing all edges ranked #1
    totaliterations=floor(nedges_2ndit/edgeremovalinterval)+1; % number of times to remove edges
    %totaliterations=5; % PRESET
    nedges(totaliterations)=0;
    
    data2{length(data)}=[];
    mask2{totaliterations}=[];
    
else
    
    for i=1:length(data)
        data{i}=data{i}(mask2{1});
    end
    
    if varyps
        totaliterations=length(ps_range);
    else
        totaliterations=1;
    end
    
end

data2=data;
for it=1:totaliterations;
    
    fprintf('Iteration %0.0f\n',it)
    
    % if masking
    if doremovals
        if it~=1
            % include everything only first time
            mask2{it}=mask>(it-2)*edgeremovalinterval+1;
            nedges(it)=sum(sum(mask2{it}));
            for i=1:length(data)
                data2{i}=data{i}(mask2{it});
            end
        end
    end
    
    % if varying ps
    if varyps
        pthresh=ps_range(it);
    else
        pthresh=pthresh_default;
    end
    
    
    % call prediction
    for i=1:nruns_range_length
        
        fprintf('Averaged over %0.0f runs\n',nruns_range(i))
        
        [m1,m2,davg]=predict_behavior(data2,behavior,ftbl,...
            'nruns',nruns_range(i),'pthresh',pthresh,'doplot',doplot,...
            'bysessions',bysessions,'doloo',doloo);
        
        detailed.r{it,i}=m1.r;
        detailed.p{it,i}=m1.p;
        detailed.predictionmask{it,i}=m1.posnegmask;
        
        if doloo
            r{1}(it,i)=m1.ppredictcorr; r{2}(it,i)=m1.npredictcorr; r{3}(it,i)=m1.fullpredictcorr;
            p{1}(it,i)=m1.ppredictpval; p{2}(it,i)=m1.npredictpval; p{3}(it,i)=m1.fullpredictpval;
            mse{1}(it,i)=m1.ppredictmse; mse{2}(it,i)=m1.npredictmse; mse{3}(it,i)=m1.fullpredictmse;
            residuals{1}{it,i}=m1.presiduals; residuals{2}{it,i}=m1.nresiduals; residuals{3}{it,i}=m1.fullresiduals;
        end
    end
    
end
detailed.mask{it}=mask2;

if doloo
    summary.r=r;
    summary.p=p;
    summary.mse=mse;
    summary.residuals=residuals;
else
    summary=nan;
end



% plot
if do_loo_plot
    rtest=cell2mat(r);
if doloo && any(~isnan(rtest(:)))
    figure;
    
    if doremovals || varyps
        if flipfigure
            x=1:nruns_range_length;
        else
            x=1:totaliterations;
        end
    else
        x=1:nruns_range_length; % TODO: careful with this; non-uniform
    end
    
    y{1}=r{1}; y{2}=r{2}; y{3}=r{3};
    y{4}=mse{1}; y{5}=mse{2}; y{6}=mse{3};
    nplots=length(y);
    
    c=colorcube;
    thesecolors=c([2,4,37,8,10,13,17,21,25,32,50,57],:);  % if iterate through color
    thesemarkers={'.'};
    markersz=22;
    
    subplot_titles={'Pos. Corr.','Neg. Corr.','Full Corr.','Pos. MSE','Neg. MSE','Full MSE'};
    ylbl_top='Correlation';
    ylbl_bottom='MSE';
    
    r=cell2mat(r);
    mse=cell2mat(mse);
    
    ymin_r=min(min(r(~isnan(r))))-0.1;
    ymax_r=max(max(r(~isnan(r))))+0.1;
    ymin_mse=min(min(mse(~isnan(mse))))-2;
    ymax_mse=max(max(mse(~isnan(mse))))+2;
    
    m_subplot=2;
    n_subplot=ceil(nplots/2);
    for i=1:nplots
        subplot(m_subplot,n_subplot,i)
        hold on
        
        % do plotting based on grouped removals or ps
        if doremovals || varyps
            if flipfigure
                xlbl='# runs';
                xticklbls=num2str(nruns_range');
                xmax=nruns_range_length;
               
                if doremovals
                    lgdids=num2str(nedges');
                    lgdtitle='# edges';
                else
                    lgdids=num2str(ps_range');
                    lgdtitle='pval';
                end
                
                
            else
                if doremovals
                    xlbl='# edges';
                    xticklbls=num2str(nedges');
                    xmax=totaliterations;
                else
                    xlbl='pval';
                    xticklbls=num2str(ps_range');
                    xmax=totaliterations;
                end
                lgdids=num2str(nruns_range');
                lgdtitle='nruns';
                
                % plot fits
                if plotfitline
                for j=1:length(nruns_range)
                    pfit = polyfit(x,y{i}(:,j)',2);
                    x1 = linspace(0,xmax);
                    y1 = polyval(pfit,x1);
                    h=plot(x1,y1);
                    set(h,'LineStyle',' -','Marker','none','Color','k')
                end
                end
                
            end
            
        else
            xlbl='# runs';
            xticklbls=num2str(nruns_range');
            xmax=nruns_range_length;
            
            if plotfitline
                if size(y{1},2)==1
                    % plot fits
                    pfit = polyfit(x,y{i},2);
                    x1 = linspace(0,xmax);
                    y1 = polyval(pfit,x1);
                    h=plot(x1,y1);
                    set(h,'LineStyle',' -','Marker','none','Color','k')

                    % get R^2 of polyfit (see https://www.mathworks.com/help/matlab/data_analysis/linear-regression.html#bswj6xx)
                    yfitvals=polyval(pfit,x);
                    yresid=y{i}-yfitvals;
                    SSresid=sum(yresid.^2);
                    SStot=(length(y{i})-1)*var(y{i}); % DOF * avg variation
                    rsq_adj(i)=1-SSresid/SStot * (length(y{i})-1)/(length(y{i})-length(pfit));
                    % RSq-adjusted: for the extra DOF of a polynomial

                end
            end
        end
        
        if i<=nplots/2
            axis([0 xmax ymin_r ymax_r])
        else
            axis([0 xmax ymin_mse ymax_mse])
        end
        
        if i==1
            ylabel(ylbl_top);
        elseif i==ceil(nplots/2)+1
            ylabel(ylbl_bottom);
        end
        
        
        % plot actual data
        if flipfigure
            h2=plot(x,y{i}','markers',markersz);
        else
            h2=plot(x,y{i},'markers',markersz);
        end
        set(h2,'DefaultAxesColorOrder',thesecolors);
        set(h2,'LineStyle','none','Marker',thesemarkers{1});
        
        % set titles
        title(subplot_titles(i))
        set(gca,'Xtick',1:xmax)
        set(gca,'Xticklabel',xticklbls)
        set(gca,'XTickLabelRotation',45)
        
        if exist('rsq_adj','var')
            if i<=nplots/2
                % set textbox (report stats)
                txtboxlocation=[i*.2805-0.05 .353 .3 .3];
                txtboxstr={sprintf('R^2 Adj=%0.4f',rsq_adj(i))};
                annotation('textbox',txtboxlocation,'String',txtboxstr,'FitBoxToText','on','LineStyle','none','FontSize',8)
                hold off
            end
            
        
        end
        
           if (varyps || doremovals) && i==nplots
                lgd=legend(h2,lgdids);
                xlabel(xlbl);
            end
        
        hold off
        
    end
    
    clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
    timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
    savefig(sprintf('tmp_beh_rplot_functionof_%s.fig',timeID));
    saveimg_squaresubplot({sprintf('tmp_beh_rplot_functionof_%s',timeID)},'no',m_subplot,n_subplot) % TODO: might be in trouble with axes
    
end
end







% % %% OFF - for dynamics
% % % you've gone too far - ignore this chunk
% % 
% % dodynamics=0;
% % if dodynamics
% %     
% %     [S{i},Q(i)]=multislice_dynamic_signed(data((i-1)*24+1:i*24),4);
% %     for i=1:12
% %         for j=1:4
% %             S2{i}{j}=S{i}{j}(:,2:6)~=S{i}{j}(:,1:5);
% %         end
% %         S3((i-1)*4+1:i*4)=+S2{i};
% %     end
% %     
% %     for i=1:48
% %         S4{i}=sum(S3{i},2);
% %     end
% %     
% %     data=S3; % data=S4
% %     
% %     % from here on is the same as the above
% %     [m1,m2,davg]=predict_behavior(data,behavior,ftbl,...
% %         'nruns',1,'pthresh',0.05,'doplot',1,'bysessions',1);
% %     
% % end
% % 
% % % you've gone too far - ignore this chunk

