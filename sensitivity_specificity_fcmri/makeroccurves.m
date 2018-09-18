function [auc,FP_FN_edges]=makeroccurves(data_orig,ftbl,varargin)
% e.g., [auc,ps_all]=makeroccurves(data,ftbl,'esz',0.2,'niterations',100,'varnruns',1);
% ftbl = [sub ... sess]
% if 'use_same_sub', draws from different sessions
% for multiple runs of healthy controls, may want:
% sshfs to :/data_dustin/ 
% [data_hc,ftbl_hc]=load_reliability_data('/Users/stephanie/Documents/data/mnt/new_controls/results/gsr_unism/','.*matrix_[12]\.nii\.gz');
%
% takes ~11 sec to run 115 iterations - running close to ideal in connectivity sim (8 s), similar to this speed: for i=1:100; [~,t]=ttest(rand(20,30000),rand(20,30000)); [~,t]=ttest(rand(20,30000),rand(20,30000)*100); end
% for now, for HC, must use >= 3 runs

tic
%% Parse input
p = inputParser;

defaultvaresz=0;  % REMOVED bc want separate section for each esz
defaultesz=nan;
defaultvarnruns=0;
defaultvarnsubs=0;
defaultnsubs=nan;
defaultniterations=10;
defaultusesamesubs=0;
defaultdoindividual=0;
defaultmakeplots=1;

addParameter(p,'varesz',defaultvaresz,@isnumeric);
addParameter(p,'esz',defaultesz,@isnumeric);
addParameter(p,'varnruns',defaultvarnruns,@isnumeric);
addParameter(p,'varnsubs',defaultvarnsubs,@isnumeric);
addParameter(p,'nsubs',defaultnsubs,@isnumeric);
addParameter(p,'niterations',defaultniterations,@isnumeric);
addParameter(p,'usesamesubs',defaultusesamesubs,@isnumeric);
addParameter(p,'doindividual',defaultdoindividual,@isnumeric);
addParameter(p,'makeplots',defaultmakeplots,@isnumeric);

parse(p,varargin{:});

varesz = p.Results.varesz;
esz_range = p.Results.esz;
varnruns = p.Results.varnruns;
varnsubs = p.Results.varnsubs;
nsubs_range = p.Results.nsubs;
niterations = p.Results.niterations;
usesamesubs = p.Results.usesamesubs;
doindividual=p.Results.doindividual;
makeplots=p.Results.makeplots;

clearvars p varargin

%% Setup factor names
% everything should be: Subject, Run, Session
if size(ftbl,2)==2 % assume "HC" but maybe make this explicit in future?
    ftbl=table(ftbl(:,1), ftbl(:,2));
    ftbl.Properties.VariableNames={'Subject','Run'};
    if sum(strcmp('Session', ftbl.Properties.VariableNames)) == 0
        ftbl.Session=ones(length(ftbl.Subject),1);
    end
    
elseif size(ftbl,2)==5 % assume "TRT
    ftbl=table(ftbl(:,1),ftbl(:,2),ftbl(:,3),ftbl(:,4),ftbl(:,5));
    ftbl.Properties.VariableNames={'Subject','Site','Day','Run','Session'};
    if ~usesamesubs
        warning('Doing test-retest data but not using same subjects. May result in limited numbers of subjects (n=6).')
    end
else
    error('Factor table structure is not recognized.')
end


%% Setup data - trimask data

nscans=length(data_orig);
matsz=size(data_orig{1});
trimask=logical(tril(ones(matsz),-1));
data=nan(nscans,sum(+trimask(:)));

for i=1:nscans
    data(i,:)=data_orig{i}(trimask)'; % This is flipped from normal since connsim was writted flipped
end
nedges=length(data(1,:));


%% Setup variables

% vary numbers of runs?
maxnruns=max(ftbl.Run);
if varnruns
    if maxnruns==1
        error('Cannot vary number of runs if only one run provided.')
    end
    
    if maxnruns/2>1
        nruns_range=[1 ceil(maxnruns/2) maxnruns];
    elseif maxnruns==2
        nruns_range=[1 2];
    end
else
    nruns_range=ceil(maxnruns/2);
end

clearvars data_orig trimask



% vary numbers of subjects?
maxnsubs=max(ftbl.Subject);
if varnsubs
    if usesamesubs
        nsubs_range=[ceil(maxnsubs/4) ceil(maxnsubs/2) maxnsubs];
    else
        nsubs_range=[ceil(maxnsubs/8) ceil(maxnsubs/4) ceil(maxnsubs/2)];
    end
else
    if isnan(nsubs_range) % check whether predefined by user
        if usesamesubs
            nsubs_range=ceil(maxnsubs/2);
        else
            nsubs_range=ceil(maxnsubs/4);
        end
    end
end

if doindividual
	warning('Individual analysis uses different nsubjects than specified. Pre/post: subjects_actual=subjects_specified-1 (e.g., 12 -> 11). Group contrast: subjects_actual=(subjects_specified*2)-1 (e.g., 20 -> 39).')
end

% vary effect size?
if varesz;
    esz_range=[0.2 0.5 0.8]; % 0.2 (sm), 0.5 (med), 0.8 (lg)
else
    if isnan(esz_range) % check whether predefined by user
        esz_range=0.5;
    end
end

% vary p thresholds
pthresh_range=[0.02:0.02:0.1 0.2 0.4 0.8]; % for now, just drop some points to make faster


%% Simulation

%preallocate
FP_edges=nan(nedges,length(pthresh_range));
FN_edges=nan(nedges,length(pthresh_range));
FP_FN_overalledges{length(nsubs_range)*length(nruns_range)*length(esz_range)}=[];
FP_FN_edges{length(nsubs_range)*length(nruns_range)*length(esz_range)}=[];
% ps_all{length(nsubs_range)*length(nruns_range)*length(esz_range)}=[];

it=1;

for nsubs=nsubs_range
    for nruns=nruns_range
        for esz=esz_range
            
            % ***** the workhorse ****
            ps=connectivitysim(data,ftbl,'iterations',niterations,'effectsize',esz,'nruns',nruns,'nsubs',nsubs,'dofalsepositive',1,'dofalsenegative',1,'usesamesubs',usesamesubs,'doindividual',doindividual);
            
            it_p=2;
            
            for pthresh=pthresh_range
                FP=ps{1} < pthresh;
                FN=ps{2} > pthresh;
                
                FP_edges(:,it_p)=mean(FP,2);
                FN_edges(:,it_p)=mean(FN,2);
                
                FP_overalledges(it_p)=mean(FP_edges(:,it_p));
                FN_overalledges(it_p)=mean(FN_edges(:,it_p));

                it_p=it_p+1;  
            end
            
            % no need to compute the first or last cols (p=0: FP=0, FN=1; p=1: FP=1,FN=0)
            FP_edges(:,1)=0; FP_overalledges(1)=0;
            FN_edges(:,1)=1; FN_overalledges(1)=1;

            FP_edges(:,it_p)=1; FP_overalledges(it_p)=1;
            FN_edges(:,it_p)=0; FN_overalledges(it_p)=0;
            
            FP_FN_overalledges{it}=[FP_overalledges' FN_overalledges'];
            FP_FN_edges{it}=[FP_edges FN_edges];
            it=it+1;
        end
    end
end



%% Calculate AUCs

% AUCs via FP/TP avgd across edges
for i=1:length(FP_FN_overalledges)
    FP_overall=FP_FN_overalledges{i}(:,1);   % x = FP = 1-TN = 1-specificity
    TP_overall=1-FP_FN_overalledges{i}(:,2); % y = TP = 1-FN = sensitivity
    auc{1}(i)=trapz(FP_overall,TP_overall);
end

% AUCs via FP/TP @ individual edges (note: FP_edges already scaled above)
% THIS ONLY DOES THE LAST PLOTTED ROC CURVE
FP_edges=FP_edges;
TP_edges=1-FN_edges;
for i=1:size(FP_edges,1)
    auc{2}(i)=trapz(FP_edges(i,:),TP_edges(i,:));
end
auc{2}(auc{2}>1)=1; % some of this stuff is >1 bc taking graph sometimes doubles back (so integral double counts)
    


%% Plot ROCs and AUC spatial map

if makeplots % TODO: change so auc's are calculated even when not plotted
    % Plot ROCs
    figure;
    hold on
    ytextbump=0.008;
    for i=1:length(FP_FN_overalledges)
        FP_overall=FP_FN_overalledges{i}(:,1);   % x = FP = 1-TN = 1-specificity
        TP_overall=1-FP_FN_overalledges{i}(:,2); % y = TP = 1-FN = sensitivity
        plot(FP_overall,TP_overall,'blue')
        if i==1
            h=text(FP_overall,TP_overall+ytextbump,num2str([0,pthresh_range,1]'),'HorizontalAlignment','right','FontSize',8);
            set(h, 'rotation', -45)
        end
    end
    plot([0 1],[0 1],'black')
    plot([0.05 0.05],[0 1],'red')
    hold off
    
    xlabel('False Positive Rate (1-Specificity)')
    ylabel('True Positive Rate (Sensitivity)')
    
    % Plot AUCs and FPR spatial maps for individual edges (does last plotted line (?))
    % note: FPs can be pretty sparse if only few iterations
    dospatialplots=1;
    thisp=4; % p column for plotting FPR; 4 corresponds with p=0.06
    if dospatialplots
        % set plotting parameters
        offset=min(auc{2})-std(auc{2})/4;
        scaling=1/(max(auc{2})-offset);
        
        offsetFP=min(+FP_edges(:,thisp))-std(+FP_edges(:,thisp))/4;
        scalingFP=1/(max(+FP_edges(:,thisp))-offsetFP);
        
        %drawmatrix_atlas((auc{2}-offset)*scaling,'reorderimg',1,'atlascategory','subnetwork','datacategory','ICC');
        domatrixsummary_avg((auc{2}-offset)*scaling,'reorderimg',1,'atlascategory','subnetwork','datacategory','ICC');
        domatrixsummary_avg((+FP_edges(:,thisp)-offsetFP)*scalingFP,'reorderimg',1,'atlascategory','subnetwork','datacategory','ICC');
    end
end



%% Plot #s of errors
% % plot at line 162 (just got ps)
% for i=1:100
% [a,b]=histcounts(ps{1}(:,i));
% a2(i,:)=a;
% end
% x=b(1:34)+0.03;
% y=mean(a2);
% err=std(a2);
% errorbar(x,y,err)
% %FDR-corrected
% theseedges=0:0.03:1;
% for i=1:100
% [fdr,q]=mafdr(ps{1}(:,i));
% [a,b]=histcounts(q,theseedges);
% a2(i,:)=a;
% end
% x=b(1:33)+0.03;
% y=mean(a2);
% err=std(a2);

toc
