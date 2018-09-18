function multiseed_result = multiseed_GLM(fhandle,varargin)
% IMPORTANT: I AM DEPENDENT--USE RUN_GLM INSTEAD!
% e.g., ex = multiseed_analysis(@calc_roi_iccs,test_data,test_factortbl,'all','site');

% TODO: THIS IS INEFFICIENT and needs to be updated

nsessions=prod(max(varargin{2}{1}));
if nsessions ~=length(varargin{1}{1})
    warning('nsess does not match length data. Will convert to length of data.')
end
    
if(strcmp(func2str(fhandle),'fit_trav_glm'))
    % turn into cell matrix with one row per seed
    %varargin{1} = rearrange_trav_cell_matrix(varargin{1}); % too much time
    %[varargin{2}{1:size(varargin{1},2)}]=deal(varargin{2}); % repeat ftbl
    
    for(i=1:1:size(varargin{1},2)) %for each data
        [multiseed_result{i}] = fhandle(varargin{1}{i},varargin{2}{i},varargin(3),varargin(4));
        
        % report percent completion
        if ~exist('amt_fin_old')
            amt_fin_old=0;
        end
        if (mod(round(i*100/size(varargin{1},2)),5)==0)
            amt_fin=round(i*100/size(varargin{1},2));
            if(~(amt_fin==amt_fin_old))
                amt_fin_old=amt_fin;
                fprintf('%0.0f%% finished\n',amt_fin);
            end
        end

    end
    
    
else  % normal 
    for(i=1:1:size(varargin{1},2))
        [multiseed_result{i}] = fhandle(varargin{:}{i});
    end
end



