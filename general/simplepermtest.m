scenario=2; % 1: tasks (2 of 3), 2: tasks (1 of 2), 3: runs (1 of 2 tasks)


switch scenario
%% perm test choosing two of three

case 1
nraters=2;
nchoices=2;
noptions=3;

nperm=10000;
nsubj=22;

for i=1:nperm
    for j=1:nsubj
    
        for rater=1:nraters
            t(rater,:)=randperm(noptions,nchoices);
        end

        % only two raters
        overlap(i,j)=isequal(t(1,:),t(2,:));
    end    
end

probe=12;

tot_overlap=sum(overlap,2);
histogram(tot_overlap)
hold on
bar(probe,100,'r')
hold off

overlaps_sorted=sort(tot_overlap,'descend');
pval=sum(overlaps_sorted>=probe)/nperm; % one-sided t-test

break;

%% perm test choosing one of two (third is fixed to be ARN)
    case 2
        
nraters=2;
nchoices=1;
noptions=2;

nperm=10000;
nsubj=22;

for i=1:nperm
    for j=1:nsubj
    
        for rater=1:nraters
            t(rater,:)=randperm(noptions,nchoices);
        end

        % only two raters
        overlap(i,j)=isequal(t(1,:),t(2,:));
    end    
end

probe=12;

tot_overlap=sum(overlap,2);
histogram(tot_overlap)
hold on
bar(probe,100,'r')
hold off

overlaps_sorted=sort(tot_overlap,'descend');
pval=sum(overlaps_sorted>=probe)/nperm; % one-sided t-test




break;

%% perm test including prob of choosing same run
    case 3
nraters=2;
nchoices=1;
noptions=5;

nchoices2=nchoices+1; %add 1 to account for two runs in ARN task

nperm=10000;
nsubj=22;

for i=1:nperm
    for j=1:nsubj
    
        for rater=1:nraters
            t(rater,:)=[randperm(noptions,nchoices) randi(2,1,(nchoices2))];
            % first two are the tasks chosen, next two are the runs
        end

    
        % only two raters
        overlap(i,j)=isequal(t(1,:),t(2,:));
    end    
end

probe=9;

tot_overlap=sum(overlap,2);
figure()
histogram(tot_overlap)
hold on
bar(probe,100,'r')
hold off

overlaps_sorted=sort(tot_overlap,'descend');
pval=sum(overlaps_sorted>=probe)/nperm; % one-sided t-test

end
break;