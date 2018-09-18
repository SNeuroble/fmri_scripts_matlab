function [mdl1,mdl2,data_avg]=predict_behavior(data,behavior_data,ftbl,varargin)
% test edgewise corr btw brain and behavior
% e.g., [m1,m2,davg]=predict_behavior(data_new,RBANS.Total,ftbl_new,'bysessions',1,'pthresh',0.05,'doloo',1,'nruns',i,'meancenter',1);
% brain organized by subs, runs within sessions
% behavior organized by people, sessions
% returns m, a structure containing the predictive edges and associated p and r vals
%   if omitting all four sessions for individual, will return 0 entries for sessions 2:4

% loo options: doloo, useloo (if doloo and bysessions, do single session each time),
% 1.1.1 predict single session for left out sub's using other sub's sessions (default): bysessions=1; useloosub=0; doloo=1;
% 1.1.2 correlation using all sessions:  bysessions=1; useloosub=0; doloo=0;
% 1.2.1 predict single session for left out sub's using other sub's sessions and loo sub's other sessions (double dipping):  bysessions=1; useloosub=1; doloo=1;
% 1.2.2 (same as 1.1.2??) correlation using all sessions: bysessions=1; useloosub=1; doloo=0;
% 2.1 predict single sub's mean for all sessions:  bysessions=0; doloo=1; (useloosub doesn't matter)
% 2.2 correlation using all subject means for all sessions: bysessions=0; doloo=0; (useloosub doesn't matter)



%% Parse input
p=inputParser;
defaultbysessions=1; % do separate sessions
defaultuseloosub=0; % use other sess in left-out sub to predict
defaultpthresh=0.05;
defaultdoloo=1;
defaultnruns=max(ftbl(:,4));
defaultispremasked=0;
defaultdosubnetworks=0;
defaultdoplot=1;

addRequired(p,'data',@iscell);
addRequired(p,'behavior_data',@isnumeric);
addRequired(p,'ftbl',@isnumeric);
addParameter(p,'bysessions',defaultbysessions,@isnumeric);
addParameter(p,'useloosub',defaultuseloosub,@isnumeric);
addParameter(p,'pthresh',defaultpthresh,@isnumeric);
addParameter(p,'doloo',defaultdoloo,@isnumeric);
addParameter(p,'nruns',defaultnruns,@isnumeric);
addParameter(p,'ispremasked',defaultispremasked,@isnumeric);
addParameter(p,'dosubnetworks',defaultdosubnetworks,@isnumeric);
addParameter(p,'doplot',defaultdoplot,@isnumeric);

parse(p,data,behavior_data,ftbl,varargin{:});

clearvars -except data behavior_data ftbl p

bysessions = p.Results.bysessions;
useloosub=p.Results.useloosub;
pthresh = p.Results.pthresh;
doloo = p.Results.doloo;
thisnruns = p.Results.nruns;
ispremasked = p.Results.ispremasked;
dosubnetworks = p.Results.dosubnetworks;
doplot = p.Results.doplot;



%% Set up data
reduced_ftbl=[ftbl(:,1),ftbl(:,4),ftbl(:,5)];
% nsites=max(ftbl(:,2)); ndays=max(ftbl(:,3));
ftbl_avg=ftbl(ftbl(:,4)==1,:);
nsubs=max(reduced_ftbl(:,1));
nruns=max(reduced_ftbl(:,2)); % thisnruns=nruns;
nsess=max(reduced_ftbl(:,3));



matdim=size(data{1});

if dosubnetworks
    summatdim=size(domatrixsummary_avg(data{1},'subnetwork',0,1,0),1);
    if ~ispremasked
        trimask=logical(tril(ones(summatdim,summatdim)));
    end
else
    if ~ispremasked
        trimask=logical(tril(ones(matdim(1),matdim(2)),-1));
    end
end

if ispremasked
    trimask=true(size(data{1}));
end


%% Average Data
if ~bysessions % btw sub - 1. avg for each sub (over all sessions and runs)
    
    meancenter=0;
    data_avg=averagecellmatrices(data,reduced_ftbl,1:nsubs,[1:thisnruns]',[1:nsess]');
    %     data_avg=averagecellmatrices(data,reduced_ftbl,1:nsubs,1,1); % single run and single sess
    
    for sub=1:nsubs
        behavior_tmp=behavior_data(ftbl_avg(:,1)==sub);
        behavior_all(sub)=nanmean(behavior_tmp);
        %         behavior_all(sub)=behavior_tmp(1); % predict single sess
    end
    behavior_all=behavior_all';
    
    
else % 2. within sub - avg for each sub/sess (over all runs)
   
    meancenter=1;
    data_avg=averagecellmatrices(data,reduced_ftbl,1:nsubs,[1:thisnruns]',1:nsess);
    %     data_avg=averagecellmatrices(data,reduced_ftbl,1:nsubs,thisnruns,1:nsess);
    behavior_all=behavior_data;
    
    % if want to, impute missing data
    nansmask=isnan(behavior_all); % remember missing data
    usemissingdata_train=0;
    if usemissingdata_train
        if any(nansmask)
            behavior_all=imputeTRTdata(behavior_all);
        end
    end
    
    % mean center
    if meancenter
        for sub=1:nsubs
            % this is to account for different subject offsets, allowing us to do a multilevel-esque model using normal correlation
            ids=(sub-1)*nsess+1:sub*nsess;
            
            % mean center data
            thisdatamean=mean(data_avg(:,ids),2);
            for idtmp=ids;
                data_avg(:,idtmp)=data_avg(:,idtmp)-thisdatamean;
            end
            
            % mean center behavior
            thisbehmean=nanmean(behavior_all(ids));
            behavior_all(ids)=behavior_all(ids)-thisbehmean;
        end
    end
    
    
end

if ~doloo
    % only need to run once bc all subs (and all sessions, if bysession) are always included
    nsubs=1;
    nsess=1;
end

clearvars behavior_tmp behavior_data data ftbl ftbl_avg tmp tmp_avg



%% Build model for each left-out sub
for sub=1:nsubs
    amt_fin_old=0;
    
    % 1. do sessions separately
    if bysessions
        
        % 1.1 predict each separate session
        % adapted for missing (NaN) data
        if ~useloosub
            
            % if doloo  1.1.1 predict single left out session
            % else      1.1.2 correlation using all sessions
            
            %             id=ftbl(:,1)==sub & ftbl(:, % TODO: should be able to do even if using selected subs/sessions
            id=((sub-1)*nsess+1):(sub*nsess);
            
            behavior_train=behavior_all;
            behavior_test=behavior_train(id); 
            testdataisnan=nansmask(id);
            
            if doloo; behavior_train(id)=[];
            end
            
            data_train=data_avg;
            data_test=data_train(:,id);
            if doloo; data_train(:,id)=[]; end
            
            mdlcounter=1;
            for i=1:nsess
               
                    id=(sub-1)*nsess+i;
                    
                    if mdlcounter==1 % build model from all other subjects; test for sess 1
                        thisnsess=nsess;
                        thismdl1=[];
                        thismdl2=[];
                    else % use above model (from all other subs) to test sess's 2-4
                        thisnsess=1;
                        thismdl1=mdl1;
                        thidmdl2=mdl2;
                    end
                    
                if ~testdataisnan(i) || ~doloo
                    [mdl1.r{sub}{i},mdl1.p{sub}{i},mdl1.posnegmask{sub}{i},...
                        mdl2.pos_mdl{sub}{i},mdl2.neg_mdl{sub}{i},mdl2.full_mdl{sub}{i},...
                        mdl2.ppredict{sub}(:,i),mdl2.npredict{sub}(:,i),mdl2.fullpredict{sub}(:,i)]...
                        =doprediction(data_train,behavior_train,data_test(:,i),...
                        pthresh,nsubs,thisnsess,doloo,thismdl1,thismdl2);
                    
                    mdlcounter=mdlcounter+1;
                else
                    mdl2.ppredict{sub}(:,i)=single(NaN);
                    mdl2.npredict{sub}(:,i)=single(NaN);
                    mdl2.fullpredict{sub}(:,i)=single(NaN);
                end
            end
            
            % 1.2 do single session each time - shouldn't really do this bc some double dipping into loo sub's data
        else
            warning('Watch out - you''re double dipping even if you mean center.')
            for sess=1:nsess
                id=(sub-1)*nsess+sess;
                
                % if doloo  1.2.1 predict single left out session
                % else      1.2.2 correlation using all sessions
                
                behavior_train=behavior_all;
                if doloo; behavior_train(id)=[];
                end
                
                data_train=data_avg;
                data_test=data_train(:,id);
                if doloo; data_train(:,id)=[]; end
                
                [mdl1.r{sub}{sess},mdl1.p{sub}{sess},mdl1.posnegmask{sub}{sess},...
                    mdl2.pos_mdl{sub}{sess},mdl2.neg_mdl{sub}{sess},mdl2.full_mdl{sub}{sess},...
                    mdl2.ppredict{sub}(:,sess),mdl2.npredict{sub}(:,sess),mdl2.fullpredict{sub}(:,sess)]...
                    =doprediction(data_train,behavior_train,data_test,...
                    pthresh,nsubs,nsess,doloo,[],[]);
                
            end; end
        
        % 2. combine over all sessions
        % adapted for missing (NaN) data - TODO set mdl NaN for full NaN
    else
        % if doloo  2.1 predict single left out session
        % else      2.2 correlation using all sessions
        
        % here, sub = leftoutsub
        
        behavior_train=behavior_all;
        behavior_test=behavior_train(sub); 
        if doloo; behavior_train(sub)=[];
        end
        
        data_train=data_avg;
        data_test=data_train(:,sub);
        if doloo; data_train(:,sub)=[];
        end
        
        if ~isnan(behavior_test)
            [mdl1.r{sub},mdl1.p{sub},mdl1.posnegmask{sub},...
                mdl2.pos_mdl{sub},mdl2.neg_mdl{sub},mdl2.full_mdl{sub},...
                mdl2.ppredict{sub},mdl2.npredict{sub},mdl2.fullpredict{sub}]...
                =doprediction(data_train,behavior_train,data_test,...
                pthresh,nsubs,1,doloo,[],[]);
        end
        
        % turn off "badly conditioned" warning (for negative vals)
        if sub==1
            w=warning('query','last');
            if ~ isempty(w); id = w.identifier;
                warning('off',id)
            end
        elseif sub==nsubs
            if exist('id','var'); warning('on',id)
            end; end
        
    end
    
    
    
    
    % report percent complete
    if (mod(round(sub*100/nsubs),ceil(nsubs/2))==0)
        amt_fin=round(sub*100/nsubs);
        if(~(amt_fin==amt_fin_old))
            amt_fin_old=amt_fin;
            fprintf('Sub%0.0f %0.0f%% finished\n',sub,amt_fin);
        end; end
    
end






% how'd we do?
if doloo
    
    beh_ppredict=cell2mat(mdl2.ppredict)';
    beh_npredict=cell2mat(mdl2.npredict)';
    beh_fullpredict=cell2mat(mdl2.fullpredict)';
    
    
    if ~exist('nansmask','var')
       nansmask=logical(zeros(size(behavior_all)));
    end

    beh_actual=behavior_all(~nansmask);
    beh_ppredict=beh_ppredict(~nansmask);
    beh_npredict=beh_npredict(~nansmask);
    beh_fullpredict=beh_fullpredict(~nansmask);
    

    [mdl1.ppredictcorr,mdl1.ppredictpval]=corr(beh_ppredict,beh_actual);
    [mdl1.npredictcorr,mdl1.npredictpval]=corr(beh_npredict,beh_actual);
    [mdl1.fullpredictcorr,mdl1.fullpredictpval]=corr(beh_fullpredict,beh_actual);
%     [mdl1.ppredictcorr,mdl1.ppredictpval]=corr(beh_ppredict,beh_actual,'type','Spearman');
%     [mdl1.npredictcorr,mdl1.npredictpval]=corr(beh_npredict,beh_actual,'type','Spearman');
%     [mdl1.fullpredictcorr,mdl1.fullpredictpval]=corr(beh_fullpredict,beh_actual,'type','Spearman');
    
    mdl1.ppredictmse=immse(double(beh_ppredict),beh_actual);
    mdl1.npredictmse=immse(double(beh_npredict),beh_actual);
    mdl1.fullpredictmse=immse(double(beh_fullpredict),beh_actual);
    
    mdl1.presiduals=abs(double(beh_ppredict-beh_actual));
    mdl1.nresiduals=abs(double(beh_npredict-beh_actual));
    mdl1.fullresiduals=abs(double(beh_fullpredict-beh_actual));
    
    plottype{1}='positive'; plottype{2}='negative'; plottype{3}='full';
    r(1)=mdl1.ppredictcorr; r(2)=mdl1.npredictcorr; r(3)=mdl1.fullpredictcorr;
    pval(1)=mdl1.ppredictpval; pval(2)=mdl1.npredictpval; pval(3)=mdl1.fullpredictpval;
    
    if doplot
        if ~bysessions
            nsess=1;
        end
        
        prediction_plot{1}=reshape(cell2mat(mdl2.ppredict),nsess,nsubs);
        prediction_plot{2}=reshape(cell2mat(mdl2.npredict),nsess,nsubs);
        prediction_plot{3}=reshape(cell2mat(mdl2.fullpredict),nsess,nsubs);
        behavior_plot=reshape(behavior_all,nsess,nsubs);
        
        figure
        
        ymin=min(behavior_all)-1;
        ymax=max(behavior_all)+1;
        if ~bysessions
            xmin=ymin+1; xmax=ymax-1;
            %             xmin=ymin+5; xmax=ymax-5;
        else
            xmin=min(min([prediction_plot{1},prediction_plot{2},prediction_plot{3}]))-0.5;
            xmax=max(max([prediction_plot{1},prediction_plot{2},prediction_plot{3}]))+0.5;
        end
        
        variablecolor=1;
        if variablecolor % vary color
            t=colorcube;
            thesecolors=t([2,4,37,8,10,13,17,21,25,32,50,57],:);  % if iterate through color
            thesemarkers={'.'};
            markersz=22;
        else % vary markers
            thesecolors=[0 0 0];
            thesemarkers={'+','o','*','.','x','s','d','^','v','>','p','h'};
            markersz=4;
        end
        
        for i=1:3
            subplot(1,3,i)
            set(0,'DefaultAxesColorOrder',thesecolors);
            set(0,'DefaultAxesLineStyleOrder',thesemarkers);
            
            % plot data
            plot(prediction_plot{i},behavior_plot,'markers',markersz);
            axis([xmin xmax ymin ymax])
            
            doindividualfitlines=0;
            if doindividualfitlines % plot automatic line for each set of datapoints:
                h=lsline;
            else % manually fit/plot LS line for all data points:
                hold on;
                nanmask=~isnan(behavior_plot);
                pfit = polyfit(prediction_plot{i}(nanmask),behavior_plot(nanmask),1);
                x1 = linspace(xmin,xmax);
                y1 = polyval(pfit,x1);
                h=plot(x1,y1);
            end
            
            axis([xmin xmax ymin ymax])
            set(h,'LineStyle',' :','Marker','none')
            axis('square')
            title(plottype(i))
            
            txtboxlocation=[i*.2805-0.148 .353 .3 .3];
            txtboxstr={sprintf('r=%0.4f',r(i)),sprintf('p=%0.4f',pval(i))};
%             txtboxstr={sprintf('r=%0.4f',r(i))};
            annotation('textbox',txtboxlocation,'String',txtboxstr,'FitBoxToText','on','LineStyle','none','FontSize',8)
            hold off
            
            %     legend('show')
            %     for i=1:nsubs
            %         plot(prediction_pplot(:,i),behavior_plot(:,i))
            %         lsline
            % %         yfit = P(1)*x+P(2);
            % %         plot(x,yfit,'r-.');
            %     end
            
            
        end
        
        % save img
        
        
        clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
        timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});
        savefig(sprintf('tmp_beh_rplot_%s.fig',timeID));
        saveas(gcf,sprintf('tmp_beh_rplot_%s.png',timeID))
        close;
        
        %         disp(sprintf('Pos corr: %0.4f, p=%0.4f',mdl1.ppredictcorr,))
        %         disp(sprintf('Neg corr: %0.4f, p=%0.4f',mdl1.npredictcorr,mdl1.npredictpval))
        %         disp(sprintf('Full model corr: %0.4f, p=%0.4f',mdl1.fullpredictcorr,mdl1.fullpredictpval))
        
        
        
    end

end

% end
