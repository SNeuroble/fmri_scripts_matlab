function plot_motion(mot)
% plots motion by run (1:6)
% load motion from ../workspace/test_retest

subplot(2,2,1)
histogram(mot)
ids_highmot=find(mot>0.3);

mot_bysess=reshape(mot,6,48);
mot_bysub=reshape(mot,24,12);

% mean-centered for each sub and session
for i=1:48
mot_bysess_meanctr(:,i)=mot_bysess(:,i)-mean(mot_bysess(:,i));
end

% mean motion by sub
meansubmot(1,:)=mean(mot_bysub);
meansubmot(2,:)=std(mot_bysub);

%% plots
% nice jittered x for runs
x=repmat([1:6]',48,1)+rand(288,1)/2-0.25;

subplot(2,2,2)
plot(reshape(mot_bysess,24,12),'-') % subjects are lines, sessxrun are x axis
% scatter(x,reshape(mot2_2,288,1),'.')
subplot(2,2,3)
scatter(x,reshape(mot_bysess_meanctr,288,1),'.') % by run
subplot(2,2,4)
boxplot(mot_bysess_meanctr') % by run

% plot corr btw sub mot and sub std
figure
scatter(meanmot(1,:),meanmot(2,:))
[r,p]=corr(meanmot(1,:)',meanmot(2,:)')

