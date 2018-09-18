function [icc_summary,var,stats,sigmask] = run_reliability(varargin)
% does sig masking and computes 1, 2- or 3-factor ICC (detects from nfactors in factors) 
% note: 2- and 3-factor ICC via G-Theory (Webb and Shavelson, 2005)
% provide data and ftbl or procedure (e.g., 'PCC')
    % data is cell array (nscans) of 1-, 2-, or 3-D matrices (any type, e.g., raw data, edges, masked data &c)
    % procedures can be: {'PCC','gbc1_unnorm','gbc1','fspectra_norm' ...}
% optionals: datadir, thismodel, correctiontype
%
% ** Ex1: [icc,var,stats,sigmask]=run_reliability('data',data,'factors',ftbl) % if loading data manually
% ** Ex2: ignore - this hasn't been updated: [icc,var,stats]=run_reliability('correctiontype','Bonf','procedure','PCC') % to get/load files locally

%% Parse input
p = inputParser;

defaultdata{1}=0; % no default
defaultfactors=0; % no default
defaultprocedure='NA'; % no default - procedure options: {'PCC', 'RMC',...} to load from directory
defaultdatadir='/Users/stephanie/Documents/data/traveling_subs/'; % automatically adds results_%procedure%
defaultthismodel='full'; % anova model options: {'full', 'none'}
defaultcorrectiontype='none'; % multiple comparison corr for sigmask options: {'none','Bonf','FDR'}

addParameter(p,'data',defaultdata,@iscell);
addParameter(p,'factors',defaultfactors,@isnumeric);
addParameter(p,'procedure',defaultprocedure,@isnumeric);
addParameter(p,'datadir',defaultdatadir,@ischar);
addParameter(p,'thismodel',defaultthismodel,@ischar);
addParameter(p,'correctiontype',defaultcorrectiontype,@ischar);

parse(p,varargin{:});

data = p.Results.data;
ftbl = p.Results.factors;
procedure = p.Results.procedure;
datadir=p.Results.datadir;
thismodel=p.Results.thismodel;
correctiontype=p.Results.correctiontype;

clearvars p varargin

%% setup data
if ~strcmp(procedure,'NA')
    warning('Loading in via file here might have to be updated.')
    % TODO: either manually replace or use pwd (after cd to datadir)
    datadir=strcat(datadir,sprintf('results_%s/',procedure));
    [data,ftbl]=load_reliability_data(datadir,'*nii*');     % TODO: update input to load_reliability_data based on new version

else
    if sum(sum(ftbl))==0 | sum(sum(data{1}))==0;
        error('Data unspecified. Either pass data/ftbl or specify procedure.')
    end
    procedure='tmp';
    
    % check for missing data (per sub)
    for i=1:size(ftbl,1);
        t(i)=sum(ftbl(:,1)==i);
    end
    if length(unique(t))>2  % 2 bc unique includes "0"
        disp('Runs per subject:')
        disp(sprintf('%d ',unique(t)))
        error('Missing data. Please remove partial data.')
    end
end

%% load data - old
% 
% thismodel='full';
% 
% if isempty(varargin)
%     error('Please either pass data/ftbl or specify procedure.')
% 
% elseif size(varargin,2)==1
%     % TODO: either manually replace or use pwd (after cd to datadir)
%     % TODO: update input to load_trt_files based on new version
%     procedure=varargin{1};
%     datadir=sprintf('/Users/stephanie/Documents/data/traveling_subs/results_%s/',procedure);
%     [data,ftbl]=load_reliability_data(datadir,'*nii*');
% 
% elseif size(varargin,2)==2
%     if ~(size(varargin{1},2)==size(varargin{2},1))
%         error('input should be: correctiontype, data, ftbl')
%     end
%     procedure='tmp';
%     data=varargin{1};
%     ftbl=varargin{2};
%     
% elseif size(varargin,2)==3
%     if ~(size(varargin{1},2)==size(varargin{2},1))
%         error('input should be: correctiontype, data, ftbl')
%     end
%     procedure='tmp';
%     data=varargin{1};
%     ftbl=varargin{2};
%     thismodel=varargin{3};
% 
% % check number of factors accounted for
% 
%     for i=1:size(ftbl,1);
% 	t(i)=sum(ftbl(:,1)==i);
%     end
% 
%     if length(unique(t))>2  % 2 bc unique includes "0"
%        
% 	disp('Runs per subject:')
% 	disp(sprintf('%d ',unique(t)))
% 	error('Missing data. Please remove partial data.')
%     end
% 
% 
% 
% %    predictedsize=1;
% %    for i=1:size(ftbl,2)
% %	  predictedsize=predictedsize*max(unique(ftbl(:,i)))
% %	%predictedsize=max(unique(ftbl(:,1)))*max(unique(ftbl(:,4)))*max(unique(ftbl(:,5)));
% %    end	
% %    if size(ftbl,1) ~= predictedsize
% %	  error(sprintf('%f runs missing. Please remove partial factors and continue.', predictedsize-size(ftbl,1)))
% %    end
% 	  %nsubj=unique(ftbl(:,1));
% 	  %for j=2:size(ftbl,2) 
% 	  %end
% 
% end
% 
% 
% 
% 
% clearvars varargin 

%% main

% determine whether doing voxelwise image or matrix

if(ndims(data{1})==3)
    doing_voxelwise=1;
elseif (ndims(data{1})==2)
    doing_voxelwise=0;
else 
    error('weird dimensions')
end

if(doing_voxelwise) % mask bc voxelwise
%    gmmask = load_untouch_nii('5thpass_template_GM_WM_CSF_resl_crop_resampled.nii.gz');
%     (if mask is in uint8, do: im2single(gmmask) )
%    gmmask=gmmask.img;
%    gmmask(gmmask<3)=0;
%    gmmask(gmmask==3)=1;

	gmmask=ones(size(data{1})); % NEW for no masking
    [masked_data,gmmask_composite]=get_masked_data_special_trav(data,gmmask);  % get GM; remove zero entries
else % no GM mask bc matrix or other
    s=size(data{1},1);
    % use only lower triangle of matrix if:
    % 1) full matrix, 2) not a single entry
    if size(data{1},2)==s && size(data{1},2)>2
        for i=1:length(data)
            masked_data{i}=data{i}(logical(tril(ones(s,s),-1)));
        end
    else % only subset of edges or other
        initmask=ones(size(data{1}));
        [masked_data,zeromask]=get_masked_data_special_trav(data,initmask); % may contain redundant edges - dof
%       masked_data=data; % SMN
    end
    
end

sigmask = create_sigedge_mask(masked_data,0.05,correctiontype); % "none"->all ones
data_sig = get_masked_data(masked_data,sigmask);
[icc_summary,var,stats]=calc_roi_iccs(data_sig,ftbl,'all',thismodel);

% if using no sigmask, still return a real sigmask in case it's needed for future analyses
if strcmp(correctiontype,'none')
    sigmask = create_sigedge_mask(masked_data,0.05,'Bonf'); % return Bonf mask anyways
end


%% save data
clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');


if(doing_voxelwise) % save mask bc voxelwise
%     save(sprintf('%s_%s.mat',procedure,correctiontype),'sigmask','stats','icc_summary','var','gmmask_composite')
    save(sprintf('%s_icc_%sH%sM%sS.mat',procedure,clockinfo{1},clockinfo{2},clockinfo{3}),'sigmask','stats','icc_summary','var','gmmask_composite')
else
%     save(sprintf('%s_%s.mat',procedure,correctiontype),'sigmask','stats','icc_summary','var')
    save(sprintf('%s_icc_%sH%sM%sS.mat',procedure,clockinfo{1},clockinfo{2},clockinfo{3}),'sigmask','stats','icc_summary','var')
end

soundthealarm=0;
if soundthealarm
    load handel
    sound(y,Fs)
end

