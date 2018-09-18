function [r,p]=imgcorr_v_interscaninterval(cmat,doplot,dobeh)
% function slopes=imgcorr_v_interscaninterval(isi,cmat)

% dobeh=1; % HARD CODED
nsubs=12; % HARD CODED

% load isi (36x1 vector of inter-scan intervals (within sub; meaningless between)
% IMPORTANT: expecting that t=0 is missing for all subs in isi
fid=fopen('../workspaces/test_retest/intersession intervals/ISI.txt');  % workspaces/test_retest/ISI.txt
isi=textscan(fid,'%d');
isi=isi{1};
fclose(fid);

% load cmat (matrix of all inter-scan correlations (48x48)
if cmat==0
    load('sandbox_cmat.mat') % workspaces/playground/sandbox_cmat.mat
end

% load behavior if applicable
if dobeh

%     load('behaviorpackagenew') % workspaces/test_retest/behaviorpackagenew.mat
%     beh=RBANS.Total;
%     beh=STAI.Y2;
    
%     % for mot
    fid=fopen('motioncorrection_allfiles_ffdonly.txt'); % workspaces/test_retest/motioncorrection_allfiles_ffdonly.txt
    beh=textscan(fid,'%f');
    beh=beh{1};
    beh=reshape(beh,6,48);
    beh=mean(beh)';
    
    beh_ref=repmat(beh,1,48);
    beh_cmat=abs(beh_ref-beh_ref');
    beh_ref=repmat(beh,1,48);
    beh_cmat=beh_ref-beh_ref';
    beh_cmat=abs(beh_cmat);
end

nsess=length(isi)/nsubs+1;
nscans=nsubs*nsess;

% fill in 0's
for i=1:nsubs
    isi=[isi(1:(i-1)*nsess); 0; isi((i-1)*nsess+1:end)];
end

% work backwards to assign scan days
scanday=isi;
for i=1:nsubs
    for j=2:nsess
        scanday((i-1)*nsess+j)=scanday((i-1)*nsess+j)+scanday((i-1)*nsess+j-1);
    end
end

% matrix of inter-scan intervals (only relevant within-sub)

isimat=repmat(scanday,1,nscans)-repmat(scanday,1,nscans)';
isimat=abs(isimat);


% mask only within-sub
mask=zeros(size(isimat));
for i=1:nsubs
    mask((i-1)*nsess+1:i*nsess,(i-1)*nsess+1:i*nsess)=1;
end
mask=logical(tril(mask,-1));
nisi=sum(sum(mask))/nsubs;

% get data in mask
X=isimat(mask);
Y=cmat(mask);

X2=reshape(X,nisi,nsubs);
Y2=reshape(Y,nisi,nsubs);

if doplot
    %similarity matrix
    m2=double(mask(1:nsess,1:nsess));
    ids=find(m2);
    Y_simmat=mean(Y2,2);
    m2(ids)=Y_simmat;
    % to flip
    m2=m2+m2'+diag(ones(size(m2,1),1));
    
    datacategory='icc'; setcolors
    figure;
    image((m2-minthresh)*80)

    cb=colorbar;
    set(get(cb,'Title'),'String',ctitle)
    set(cb,'YTick',cticks); set(cb,'YTickLabel',cticklbls)
    set(cb, 'ylim', cylim)
end

order=[1,4,6,2,5,3]; % HARD CODED to reflect reasonable pairs; last is sess1-sess4
npairs=length(order);
X2=single(X2(order,:));
Y2=Y2(order,:);

% ColorSet = varycolor(50);
% set(gca, 'ColorOrder', ColorSet);

[r,p]=corr(X2(:),Y2(:));
if doplot
    % figure
    % To add line plots - note that this messes up the scatterplot legend colors
    figure
    %set_better_scatter_colors

    for i=1:nsubs
        %     plot(X2([1,find(order==3)],i),Y2([1,find(order==3)],i),'-k','LineWidth',0.1)
        fit=polyfit(X2(:,i),Y2(:,i),1);
        slopes(i)=fit(1);
        plot(X2(:,i),polyval(fit,X2(:,i)),'-k','LineWidth',0.1)
        hold on
    end

    X2=X2'; Y2=Y2';
    for i=1:nisi
        scatter(X2(:,i),Y2(:,i),50,'filled')
        hold on
    end

    
    % should first plot lines connecting ppl

    str={'sess1-sess2','sess1-sess3','sess1-sess4','sess2-sess3','sess2-sess4','sess3-sess4'};
    str=str(order);
    str=horzcat({'s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11','s12'},str);
    legend(str,'Location','SouthEast');
    
    
    fprintf('Connectivity sim v. interval: r=%0.4f, p=%0.4f\n',r,p);
end


% corr of mat similarity and behavior similarity
if dobeh
    
    mask=zeros(48);
    for i=1:12
        mask((i-1)*4+1:i*4,(i-1)*4+1:i*4)=1;
    end
    mask=logical(tril(mask,-1));
    
    mask2=isnan(beh_cmat);
    mask=mask&(~mask2);
    
    cmat2=cmat(mask);
    beh_cmat2=beh_cmat(mask);
    
    [r2,p2]=corr(cmat2,beh_cmat2);
        
    if doplot
        set_better_scatter_colors
        for i=1:12
            scatter(cmat2((i-1)*6+1:i*6),beh_cmat2((i-1)*6+1:i*6))
        end
        
         fprintf('Behavior sim v. connectivity sim: r=%0.4f, p=%0.4f\n',r2,p2);
    end
    

    % figure; scatter(cmat(mask),beh_cmat((mask)))
    % multiple regression - try MLM
    isimat2=single(isimat(mask));
    X=[isimat2 beh_cmat2];
    fitlm(X,cmat2)
%     fitlme(data,'behav ~ 1 + sess * connectivity + (1|subj)');
    
% get minicorrs
%     littleprod=minicorrs(isimat2,cmat2);
end

