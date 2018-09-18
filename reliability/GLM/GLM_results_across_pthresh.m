%{
clearvars counts_master
% TODO: switch matrix and icd next time this runs
% procedures={'pcc','rmc','lt','matrix','icd'};
procedures={'icd'};
factors={'subj','site','day'};

it=1;
for this_procedure=procedures
    for this_factor=factors
        
         file=sprintf('%s_glm_%s_all.mat',this_procedure{1},this_factor{1});
         load(file,'glm_data')
         ps = get_multiseed_ps(glm_data);
        [fdr_ps,qs]=fdr_correct_multiseed_pvals(ps);
        
        
        ps=cell2mat(ps(2:end)); % first col is intercept
        thresh_ps=threshold_multiseed_pvals(ps);
        counts=0;
        for i=2:size(thresh_ps,2)
            for j=2:size(thresh_ps,1) % p-thresh = 0.1, 0.05, 0.01, 0.005, 0.001
            counts(i-1,j-1)=sum(thresh_ps{j,i});
            end
        end
        
        
        dim=size(counts,2);
        thresh_ps=threshold_multiseed_pvals(qs);
        for i=2:size(thresh_ps,2)
            for j=2:size(thresh_ps,1) % q-thresh
            counts(i-1,j-1+dim)=sum(thresh_ps{j,i});
            end
        end
        
        counts_master{1,it}=sprintf('%s_glm_%s',this_procedure{1},this_factor{1});
        counts_master{2,it}=counts;
        
        it=it+1;
        
        clearvars -except ps qs counts_master it this_procedure procedures this_factor factors
        
    end
end 

clearvars it this_procedure procedures this_factor factors
%}

% can just load counts_master and run from here
tmean=[];
terr=[];

n=[42784,42859,42784,38503,42733];

figure
for i=3
% for i=1:3 % subj site day
    tmean=[];
    terr=[];
    
    
    for j=i:3:length(counts_master)
        tmean=[tmean;mean(counts_master{2,j})];
        terr=[terr;std(counts_master{2,j})];
    end
    
    for k=1:size(tmean,1)
        tmean2(k,:)=tmean(k,:)./n(k)*100;
        terr2(k,:)=terr(k,:)./n(k)*100;
    end
    
    tmean_master{1,i}=tmean2;
    tmean_master{2,i}=terr2;
    
    for m=1:2
        subplot(3,2,2*(i-1)+m)
        id=5*(m-1)+1;
        id2=id+4;
        
        errorbar(tmean2(:,id:id2)',terr2(:,id:id2)');
%         errorbar(tmean',terr');

        xmin=0; xmax=6; ymin=-5; ymax=50;
        axis([xmin xmax ymin ymax])
        if m==1
            xlabel('p-val')
            set(gca,'XTickLabel',{'','p<0.1','p<0.05','p<0.01','p<0.005','p<0.001',''}) %shows 1 to 11
        elseif m==2
            xlabel('q-val')
            set(gca,'XTickLabel',{'','q<0.1','q<0.05','q<0.01','q<0.005','q<0.001', ''}) %shows 1 to 11
        end
        
        
        ylabel('Percent')
        if i==1
            title('Subj')  
        elseif i==2
            title('Site') 
        elseif i==3
            title('Day') 
        end

        % blue PCC red RMC orange LT purple matrix green ICD
    end
end


