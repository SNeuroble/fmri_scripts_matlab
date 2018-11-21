function [data,ftbl] = load_reliability_data(thisFolder, thisPattern)
% built to parse metadata from traveling_subs and test_retest
% supports regex
% e.g., [data,ftbl] = load_trt_files('/Users/stephanie/Documents/data/mnt/gsr_unism/','.*matrix_.\.nii\.gz');
% out:
%   ftbl(subj,scanner,day)
%   data = cell matrix; each cell contains an image vector corresponding with a file
% remember! for trav seeds, use /more_results/pcc/*crop_resampled.nii.gz

flag_studytype=0;

if ~isdir(thisFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', thisFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% filePattern = fullfile(thisFolder, thisPattern);
% filePattern = strcat(thisFolder, thisPattern);
% theseFiles = dir(filePattern);

theseFiles = regexpdir(thisFolder, thisPattern);

if isempty(theseFiles)
  errorMessage = sprintf('Error: Nothing found matching\n%s', thisPattern);
  uiwait(warndlg(errorMessage));
  return;
end


for k = 1:length(theseFiles)
%     baseFileName = theseFiles(k).name;
%     fullFileName = fullfile(thisFolder, baseFileName);
      fullFileName = theseFiles{k};

          % report percent completion
        if ~exist('amt_fin_old')
            amt_fin_old=0;
        end
        if (mod(round(k*100/length(theseFiles)),5)==0)
            amt_fin=round(k*100/length(theseFiles));
            if(~(amt_fin==amt_fin_old))
                amt_fin_old=amt_fin;
                fprintf('%0.0f%% finished\n',amt_fin);
            end
        end

%          fprintf(1, 'Now reading %s\n', fullFileName);
      
    [pathstr,name,ext] = fileparts(fullFileName);

    baseFileName=[name,ext];

    % check correct filetype and assign study_type
    if ~flag_studytype
        if(~isempty(strfind(baseFileName,'000_')))
            study_type='traveling_subs';
        elseif (~isempty(strfind(baseFileName,'TRT')))
            study_type='test_retest';
        else
            study_type='undefined';
        end
        if ~(strfind(ext,'.nii')) & ~(strfind(name,'.nii'))
            errorMessage = sprintf('Error: The file is not .nii or .nii.gz:\n%s', fullFileName);
            uiwait(warndlg(errorMessage));
            return;
        end
        flag_studytype=1;
    end
  

    % if(strfind(ext,'.txt'))
    % don't use this. scripts built to use output of nii.
    % see final_project for info on squeezing out singletons
    %       a=readtable(fullFileName,'delimiter', 'tab');
    %       a(:,3) = [];
    %       output{k} = a; 

    a = load_untouch_nii(fullFileName);
    data{k} = a.img;
  
    if strcmp(study_type,'traveling_subs')
        % file prototype: 01_V3000_2_bis_matrix.nii
        % ftbl: subj scanner day (within a particular scanner)
        identifier=strfind(baseFileName,'000_');
        this_name=baseFileName((identifier-5):(identifier+4)); % returns "01_V1000_1"
        this_subj=this_name(5);
        this_scanner=this_name(2);
        this_day=this_name(10);

        ftbl(k,1)=str2num(this_subj);
        ftbl(k,2)=str2num(this_scanner);
        ftbl(k,3)=str2num(this_day);

    elseif strcmp(study_type,'test_retest')
        % file prototype: TRT001_1_TA_S001_bis_matrix.nii
        % ftbl: subj scanner day (within a particular scanner) run session
        identifier=strfind(baseFileName,'TRT');
        this_name=baseFileName((identifier+4):(identifier+15)); % returns "01_1_TA_S001"
        this_subj=this_name(1:2);
        this_session=this_name(4);
        this_run=this_name(12);
        
        % coding scanner 1=TA, 2=TB
        if(strcmp(this_name(6:7),'TA')); this_scanner='1';
        else this_scanner='2';
        end

        this_day=this_name(10);

        ftbl(k,1)=str2num(this_subj);
        ftbl(k,2)=str2num(this_scanner);
        % day/occasion is added after reading all
        if ~(strcmp(this_run,'_'))
                ftbl(k,4)=str2num(this_run);
            else ftbl(k,4)=1;
        end
        ftbl(k,5)=str2num(this_session);

        else ftbl=[];
    end
    
end

% assign day/occasion
if strcmp(study_type,'test_retest')
    for subj=unique(ftbl(:,1))'
        for scanner=unique(ftbl(:,2))'
            
            t=(find(ftbl(:,1)==subj & ftbl(:,2)==scanner));
            t(:,2)=ftbl(t,5);
            
            l=unique(t(:,2));
            l=sortrows(l);
            for i=1:length(l);
                t(t(:,2)==l(i),2)=i;
            end
            
            ftbl(t(:,1),3)=t(:,2);
            
        end
    end
end
  
end
