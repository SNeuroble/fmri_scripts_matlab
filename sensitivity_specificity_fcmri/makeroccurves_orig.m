function [auc,ps_allcounts]=makeroccurves(data_orig,ftbl,varargin)
% e.g., [auc,ps_all]=makeroccurves(data,ftbl,'esz',0.2,'niterations',100,'varnruns',1);
% ftbl = [sub ... sess]
% if 'use_same_sub', draws from different sessions
% for multiple runs of healthy controls, may want:
% [data_hc,ftbl_hc]=load_reliability_data('/Users/stephanie/Documents/data/mnt/new_controls/results/gsr_unism/','.*matrix_[12]\.nii\.gz');

% takes ~15 sec to run 115 iterations. Now 20 sec for 135 iterations.
% for now, for HC, must use >= 3 runs

% TODO: probably just need to simulate once and sweep through

tic
%% Parse input
p = inputParser;

defaultvaresz=0;  % REMOVED bc want separate section for each esz
defaultesz=nan;
defaultvarnruns=0;
defaultvarnsubs=0;
defaultnsubs=nan;
defaultniterations=5; % IMPORTANT: pvals>=0.1 use only 2 iterations (we are less interested in these)
% 7 iterations works kinda well at least to remove jagged parts
defaultusesamesubs=0;
defaultdoindividual=0;


addParameter(p,'varesz',defaultvaresz,@isnumeric);
addParameter(p,'esz',defaultesz,@isnumeric);
addParameter(p,'varnruns',defaultvarnruns,@isnumeric);
addParameter(p,'varnsubs',defaultvarnsubs,@isnumeric);
addParameter(p,'nsubs',defaultnsubs,@isnumeric);
addParameter(p,'niterations',defaultniterations,@isnumeric);
addParameter(p,'usesamesubs',defaultusesamesubs,@isnumeric);
addParameter(p,'doindividual',defaultdoindividual,@isnumeric);


parse(p,varargin{:});

varesz = p.Results.varesz;
esz_range = p.Results.esz;
varnruns = p.Results.varnruns;
varnsubs = p.Results.varnsubs;
nsubs_range = p.Results.nsubs;
niterations = p.Results.niterations;
usesamesubs = p.Results.usesamesubs;
doindividual=p.Results.doindividual;

clearvars p varargin


%% Setup factor names
% everything should be: Subject, Run, Session
if size(ftbl,2)==2 % assume "HC" but maybe make this explicit in future?
    ftbl=table(ftbl(:,1), ftbl(:,2));
    ftbl.Properties.VariableNames={'Subject','Run'};
    
    if sum(strcmp('Session', ftbl.Properties.VariableNames)) == 0
        nsess=1;
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


%% Setup data

% trimask data
nscans=length(data_orig);
matsz=size(data_orig{1});
trimask=logical(tril(ones(matsz),-1));
for i=1:nscans
    data(i,:)=data_orig{i}(trimask)'; % This is flipped from normal since connsim was writted flipped
end
nedges=length(data(1,:));

% create new ftbl
% ftbl_new=ftbl;
% ftbl_new=[ftbl.Subject ftbl.Run ftbl.Session];


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
if varnsubs;
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

% vary effect size?
if varesz;
    esz_range=[0.2 0.5 0.8]; % 0.2 (sm), 0.5 (med), 0.8 (lg)
else
    if isnan(esz_range) % check whether predefined by user
        esz_range=0.5;
    end
end

pthresh_range=[0.02:0.02:0.1 0.2 0.4 0.8]; % for now, just drop some points to make faster
pthresh_iterationslimit=0.1;

%% Simulation

%preallocate
f2{length(nsubs_range)*length(nruns_range)*length(esz_range)}=[];
% ps_all{length(nsubs_range)*length(nruns_range)*length(esz_range)}=[];

it=1;

for nsubs=nsubs_range
    for nruns=nruns_range
        for esz=esz_range
            
            counts=NaN(niterations,2*length(pthresh_range));
            % don't need to compute the first two cols (FP=0,FN=1)
            counts=[zeros(niterations,1),ones(niterations,1)];
            ps_allcounts=[zeros(nedges,1),ones(nedges,1)];
            it_ct=3;
            
            for pthresh=pthresh_range
                % to speed things up, do fewer iterations for larger (less interesting) pvals
                if pthresh<=pthresh_iterationslimit; niterations_mod=niterations;
                else niterations_mod=10;
                end
                
                % ***** the workhorse ****
                [ct,ps]=connectivitysim(data,ftbl,'iterations',niterations_mod,'effectsize',esz,'pthresh',pthresh,'nruns',nruns,'nsubs',nsubs,'dofalsepositive',1,'dofalsenegative',1,'usesamesubs',usesamesubs,'doindividual',doindividual);
                % consider adding: 'names',{'Subject' 'Run' 'Session'}
                
                
                if pthresh>pthresh_iterationslimit % if doing fewer iterations
                    ct=[ct;repmat(mean(ct),niterations-niterations_mod,1)];
                end
                
                counts(:,it_ct:it_ct+1)=ct/nedges; % [FP FN FP ...]
                ps_allcounts(:,it_ct)  = mean((ps{1} < pthresh),2); % [FP FN FP ...] % FP
                ps_allcounts(:,it_ct+1)= mean((ps{2} > pthresh),2); % FN
                it_ct=it_ct+2;
            end
            
            % no need to compute the last two cols either (FP=1,FN=0)
            counts(:,it_ct:it_ct+1)=[ones(niterations,1),zeros(niterations,1)];
            ps_allcounts(:,it_ct:it_ct+1)=[ones(nedges,1),zeros(nedges,1)];
            
            f2{it}=reshape(mean(counts),2,size(counts,2)/2)'; % [FP FN]
            it=it+1;
        end
    end
    
end


%% Plot ROCs

figure;
hold on
ytextbump=0.008;
fsz=8;
for i=1:length(f2)
    x=f2{i}(:,1);   % x = FP = 1-TN = 1-specificity
    y=1-f2{i}(:,2); % y = TP = 1-FN = sensitivity
    auc{1}(i)=trapz(x,y);
    plot(x,y,'blue') %  (special for this case since 100% should be true pos)
    if i==1
        h=text(x,y+ytextbump,num2str([0,pthresh_range,1]'),'HorizontalAlignment','right','FontSize',8);
        set(h, 'rotation', -45)
    end
end
plot([0 1],[0 1],'black')
plot([0.05 0.05],[0 1],'red')
hold off

xlabel('False Positive Rate (1-Specificity)')
ylabel('True Positive Rate (Sensitivity)')


%% Spatial AUCs and FPR
% does last plotted line
dospatialplots=1;
thisp=4; % p column for plotting FPR; 4 corresponds with p=0.06



if dospatialplots
     x=ps_allcounts(:,1:2:end); % already scaled by nedges above
     y=1-ps_allcounts(:,2:2:end);
    for i=1:size(ps_allcounts,1)
        auc{2}(i)=trapz(x(i,:),y(i,:));
    end
    
    auc{2}(auc{2}>1)=1; % some of this stuff is >1 bc taking graph sometimes doubles back (so integral double counts)
    
    % set plotting parameters
    offset=min(auc{2})-0.05;
    offset2=min(+x(:,thisp))-0.05;
    scaling=1/(max(auc{2})-offset);
    scaling2=1/(max(+x(:,thisp))-offset2);
    
    %drawmatrix_atlas(auc,'reorderimg',1,'atlascategory','subnetwork','datacategory','ICC');
    domatrixsummary_avg((auc{2}-offset)*scaling,'reorderimg',1,'atlascategory','subnetwork','datacategory','ICC');
    domatrixsummary_avg((+x(:,thisp)-offset2)*scaling2,'reorderimg',1,'atlascategory','subnetwork','datacategory','ICC');
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
