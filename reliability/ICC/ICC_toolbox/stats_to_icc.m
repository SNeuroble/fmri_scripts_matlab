function [icc,var]=stats_to_icc(stats,factor_tbl,DSrange,Dstepsz)
% computes 1, 2- or 3-factor ICC (detects from nfactors in ftbl) 
% note: 2- and 3-factor ICC are G-Theory ICC (Webb and Shavelson, 2005)
% note: ftbl structure doesn't affect G-Theory ICC. Can input factor_tbl=0.
% uses subj as reference by default; can optionally do for site/day factors

% to do sigmask-ed:
% 1) load('sigmask','stats'), 2) mask the stats, and 3) use this script

tempsave=0; % turn on to save ICC if running directly, not within ICC script

if mod(DSrange, Dstepsz) ~= 0 || mod(DSrange, Dstepsz) ~= 0
    error('D-Study is misspecified.')
end

set_neg_var_zero=1;


nvarcomp=size(stats{1},1)-2;
nreffactors=1; % 1=subj, 2=[subj,scanner], 3=[subj,scanner,day]
% since it's not intended to be done over subjs, have to do some work with
% saving vars if want to add scanner and day
% for i=1:size(factor_tbl,2)
%    factor_size(i)=size(unique(factor_tbl(:,i)),1);
% end

if factor_tbl ~= 0
    f1size=length(unique(factor_tbl(:,1)));
else
    f1size=0;
end


var_comp{nvarcomp}=[];
for i=1:nvarcomp
    
    var_comp{i}=getvalsfromcellarray(stats,i+1,5,1);
    
    if set_neg_var_zero % set negative variances to zero
        var_comp{i}(var_comp{i}<0)=0;
    end
    
    var(1,i)=nanmean(nanmean(var_comp{i}));
end

% get mean % variance
for i=1:length(var_comp{1})
    var_tmp=getvalsfromcellarray(var_comp,i,1,1);
    varperc(i,:)=var_tmp/sum(var_tmp)*100;
end

var(2,:)=mean(varperc,1);
clearvars varperc


if nvarcomp==2
    % classical ICC based on mean squares; no D-Study
    % NOTE: output "var" is actually mean of mean squares
    
    testsubj=min(factor_tbl(:,1));
    k=sum(factor_tbl(:,testsubj)==1); % get number of observations for a single subj
    
    for i=1:nvarcomp
        MS{i}=getvalsfromicccellarray(stats,i+1,2);
        if set_neg_var_zero % set negative variances to zero
            MS{i}(MS{i}<0)=0;
        end
        var(i)=nanmean(nanmean(MS{i}));
    end
    
    
    for i=1:length(stats)
        icc(i)=(MS{1}(i)-MS{2}(i))/(MS{1}(i)+(k-1)*MS{2}(i));
    end
    
    
    return;
    
    
    
elseif nvarcomp==3
    % from 2-way ANOVA, pairwise interaction+residual
    % composition:  1: var_su         2: var_sc       3: var_su_sc_e
    
    %     error('Have not designed D-Study for only subj and sc')
    
    ref_factor=1;
    n1_all=[1:Dstepsz:DSrange];
    
    Dstudy_G=zeros(1,length(n1_all));
    Dstudy_D=zeros(1,length(n1_all));
    icc_G{ref_factor}=zeros(size(var_comp{1},1),size(var_comp{1},2));
    icc_D{ref_factor}=zeros(size(var_comp{1},1),size(var_comp{1},2));
    
    it1=1;
    for n1=n1_all % factor 2
        factor_size=[f1size n1];
        
        for i=1:size(var_comp{1},1)
            for j=1:size(var_comp{1},2)
                
                v(1)=var_comp{1}(i,j); %
                v(2)=var_comp{2}(i,j)/factor_size(2);
                v(3)=var_comp{3}(i,j)/factor_size(2);

                icc_G{ref_factor}(i,j)=v(1)/(v(1)+v(3));
                icc_D{ref_factor}(i,j)=v(1)/(sum(v));
                
                % OLD MANUAL
%                 icc_G{ref_factor}(i,j)=var_comp{1}(i,j)/(var_comp{1}(i,j)+var_comp{3}(i,j)/factor_size(2));
%                 icc_D{ref_factor}(i,j)=var_comp{1}(i,j)/(var_comp{1}(i,j)+var_comp{2}(i,j)/factor_size(2)+var_comp{3}(i,j)/factor_size(2));
            end
        end
        
        %         icc_G_overall=var(1)/(var(1)+var(3)/factor_size(2));
        %         icc_D_overall=var(1)/(var(1)+var(2)/factor_size(2)+var(3)/factor_size(2));
        
        Dstudy_G(it1)=nanmean(nanmean(icc_G{ref_factor})); % half new UPDATE
        Dstudy_D(it1)=nanmean(nanmean(icc_D{ref_factor})); % half new UPDATE
        
        if it1==1 %new
            first_iccG=icc_G;
            first_iccD=icc_D;
        end
        
        if it1 < 10 %new
            detailed_icc_G(it1)=icc_G;
            detailed_icc_D(it1)=icc_D;
        end
        
        it1=it1+1;
    end %new
    
    
    
    
elseif nvarcomp==7
    % from 3-way ANOVA, pairwise interactions, triple interaction+residual
    % composition:  1: var_su         2: var_sc       3: var_d
    %               4: var_su_sc      5: var_su_d     6: var_sc_d   7: var_su_sc_d_e
    
    %     icc_G=[{zeros(size(var_comp{1},1),size(var_comp{1},2))} {zeros(size(var_comp{1},1),size(var_comp{1},2))} {zeros(size(var_comp{1},1),size(var_comp{1},2))}];
    %     icc_D=[{zeros(size(var_comp{1},1),size(var_comp{1},2))} {zeros(size(var_comp{1},1),size(var_comp{1},2))} {zeros(size(var_comp{1},1),size(var_comp{1},2))}];
    
    for ref_factor=1:nreffactors % for each factor: subj, scanner, day
        if ref_factor==1
            a=1;
            b=2;
            c=3;
            d=4;
            e=5;
            f=6;
        elseif ref_factor==2
            a=2;
            b=3;
            c=1;
            d=6;
            e=4;
            f=5;
        else
            a=3;
            b=1;
            c=2;
            d=5;
            e=6;
            f=4;
        end
        g=7;
        
        n1_all=[1:Dstepsz:DSrange];
        n2_all=[1:Dstepsz:DSrange];
        
        Dstudy_G=zeros(length(n1_all),length(n2_all));
        Dstudy_D=zeros(length(n1_all),length(n2_all));
        icc_G{ref_factor}=zeros(size(var_comp{1},1),size(var_comp{1},2));
        icc_D{ref_factor}=zeros(size(var_comp{1},1),size(var_comp{1},2));
        
        it1=1;
        for n1=n1_all % factor 2
            
            it2=1;
            for n2=n2_all % factor 3
                factor_size=[f1size n1 n2];

                
                for i=1:size(var_comp{1},1)
                    for j=1:size(var_comp{1},2)
                        
                        v(1)=var_comp{a}(i,j); %
                        v(2)=var_comp{b}(i,j)/factor_size(b);
                        v(3)=var_comp{c}(i,j)/factor_size(c);
                        v(4)=var_comp{d}(i,j)/factor_size(b); %
                        v(5)=var_comp{e}(i,j)/factor_size(c); %
                        v(6)=var_comp{f}(i,j)/(factor_size(b)*factor_size(c));
                        v(7)=var_comp{g}(i,j)/(factor_size(b)*factor_size(c)); %
                        
                        icc_G{ref_factor}(i,j)=v(1)/(v(1)+v(4)+v(5)+v(7));
                        icc_D{ref_factor}(i,j)=v(1)/(sum(v));
                        
                        % OLD MANUAL
%                         icc_G{ref_factor}(i,j)=var_comp{a}(i,j)/(var_comp{a}(i,j)+var_comp{d}(i,j)/factor_size(b)+var_comp{e}(i,j)/factor_size(c)+var_comp{g}(i,j)/(factor_size(b)*factor_size(c)));
%                         icc_D{ref_factor}(i,j)=var_comp{a}(i,j)/(var_comp{a}(i,j)+var_comp{b}(i,j)/factor_size(b)+var_comp{c}(i,j)/factor_size(c)+var_comp{d}(i,j)/factor_size(b)+var_comp{e}(i,j)/factor_size(c)+var_comp{f}(i,j)/(factor_size(b)*factor_size(c))+var_comp{g}(i,j)/(factor_size(b)*factor_size(c)));
                        
                    end
                end
                
                % this does ICC using mean var components, and is unnecess.
                %             icc_G_overall{ref_factor}=var_comp_mean(a)/(var_comp_mean(a)+var_comp_mean(d)/factor_size(b)+var_comp_mean(e)/factor_size(c)+var_comp_mean(g)/(factor_size(b)*factor_size(c)));
                %             icc_D_overall{ref_factor}=var_comp_mean(a)/(var_comp_mean(a)+var_comp_mean(b)/factor_size(b)+var_comp_mean(c)/factor_size(c)+var_comp_mean(d)/factor_size(b)+var_comp_mean(e)/factor_size(c)+var_comp_mean(f)/(factor_size(b)*factor_size(c))+var_comp_mean(g)/(factor_size(b)*factor_size(c)));
                %             icc_G_mean{ref_factor} = nanmean(nanmean(icc_G{ref_factor}));
                %             icc_D_mean{ref_factor} = nanmean(nanmean(icc_D{ref_factor}));
                % >> end old stuff
                
                Dstudy_G(it1,it2)=nanmean(nanmean(icc_G{ref_factor})); % half new
                Dstudy_D(it1,it2)=nanmean(nanmean(icc_D{ref_factor})); % half new
                
                if it1==1 && it2==1 %new
                    first_iccG=icc_G;
                    first_iccD=icc_D;
                end
                
                detailedlim=10;
%                 if it1<=detailedlim/Dstepsz && it2==1 %new
                if it1<= DSrange/Dstepsz && it2==1 %new
                    detailed_icc_G(it1,it2)=icc_G;
                    detailed_icc_D(it1,it2)=icc_D;
                end
                
                it2=it2+1;
            end %new
            
            it1=it1+1;
        end %new
        
        
    end
else
    
    % assume subj x run x scanner x day
    %                 var_comp=tbl(2:15,13); % variance components
    %                 % composition:  1: var_su         2: var_r          3: var_sc        4: var_d
    %                 %               5: var_su_r       6: var_su_sc      7: var_su_d      8: var_r_sc    9: var_r_d	10: var_sc_d       
    %                 %               11: var_su_r_sc   12: var_su_r_d    13:var_su_sc_d   14: var_r_sc_d
    %                 %               15: var_su_r_sc_d
    
    ref_factor=1;
    
    n1_all=[1:Dstepsz:DSrange];
    n2_all=[1:Dstepsz:DSrange];
    n3_all=[1:Dstepsz:DSrange];
    
    Dstudy_G=zeros(length(n1_all),length(n2_all),length(n3_all));
    Dstudy_D=zeros(length(n1_all),length(n2_all),length(n3_all));
    icc_G{ref_factor}=zeros(size(var_comp{1},1),size(var_comp{1},2));
    icc_D{ref_factor}=zeros(size(var_comp{1},1),size(var_comp{1},2));
    
    it1=1;
    for n1=n1_all % factor 2
        it2=1;
        for n2=n2_all % factor 3
            it3=1;
            for n3=n3_all % factor 3
                
                factor_size=[f1size n1 n2 n3];
                
                for i=1:size(var_comp{1},1)
                    for j=1:size(var_comp{1},2)
                        
                        
                        if size(var_comp,2)==12
                            %12 components, from model
                            v(1)=var_comp{1}(i,j); %
                            v(2)=var_comp{2}(i,j)/factor_size(2);
                            v(3)=var_comp{3}(i,j)/factor_size(3);
                            v(4)=var_comp{4}(i,j)/factor_size(4);
                            v(5)=var_comp{5}(i,j)/factor_size(2); %
                            v(6)=var_comp{6}(i,j)/factor_size(3); %
                            v(7)=var_comp{7}(i,j)/factor_size(4); %
                            v(8)=var_comp{8}(i,j)/(factor_size(2)*factor_size(3));
                            v(9)=var_comp{9}(i,j)/(factor_size(3)*factor_size(4));
                            v(10)=var_comp{10}(i,j)/(factor_size(2)*factor_size(3)); %
                            v(11)=var_comp{11}(i,j)/(factor_size(3)*factor_size(4)); %
                            v(12)=var_comp{12}(i,j)/(factor_size(2)*factor_size(3)*factor_size(4)); %
                            
                            icc_G{ref_factor}(i,j)=v(1)/(v(1)+v(5)+v(6)+v(7)+v(10)+v(11)+v(12));
                            icc_D{ref_factor}(i,j)=v(1)/(sum(v));
                            
                        elseif size(var_comp,2)==15
                            %15 components, automatically crossed
                            v(1)=var_comp{1}(i,j); %
                            v(2)=var_comp{2}(i,j)/factor_size(2);
                            v(3)=var_comp{3}(i,j)/factor_size(3);
                            v(4)=var_comp{4}(i,j)/factor_size(4);
                            v(5)=var_comp{5}(i,j)/factor_size(2); %
                            v(6)=var_comp{6}(i,j)/factor_size(3); %
                            v(7)=var_comp{7}(i,j)/factor_size(4); %
                            v(8)=var_comp{8}(i,j)/(factor_size(2)*factor_size(3));
                            v(9)=var_comp{9}(i,j)/(factor_size(2)*factor_size(4));
                            v(10)=var_comp{10}(i,j)/(factor_size(3)*factor_size(4));
                            v(11)=var_comp{11}(i,j)/(factor_size(2)*factor_size(3)); %
                            v(12)=var_comp{12}(i,j)/(factor_size(2)*factor_size(4)); %
                            v(13)=var_comp{13}(i,j)/(factor_size(3)*factor_size(4)); %
                            v(14)=var_comp{14}(i,j)/(factor_size(2)*factor_size(3)*factor_size(4));
                            v(15)=var_comp{15}(i,j)/(factor_size(2)*factor_size(3)*factor_size(4)); %
                            
                            icc_G{ref_factor}(i,j)=v(1)/(v(1)+v(5)+v(6)+v(7)+v(11)+v(12)+v(13)+v(15));
                            icc_D{ref_factor}(i,j)=v(1)/(sum(v));
                            
                        end
                        
                    end
                end
                
                
                Dstudy_G(it1,it2,it3)=nanmean(nanmean(nanmean(icc_G{ref_factor}))); % half new
                Dstudy_D(it1,it2,it3)=nanmean(nanmean(nanmean(icc_D{ref_factor}))); % half new
                
                if it1==1 && it2==1 && it3==1 %new
                    first_iccG=icc_G;
                    first_iccD=icc_D;
                end
                
                if it1==1 && it2 < 3 && it3 < 3 %new
                    detailed_icc_G(it1,it2)=icc_G;
                    detailed_icc_D(it1,it2)=icc_D;
                end
                
                it3=it3+1;
            end % new
            
            it2=it2+1;
        end %new
        
        it1=it1+1;
    end %new
    
    
    
    
end




icc = {'DStudy_Gmean','DStudy_Dmean','first_Gmap (n_s=1,n_d=1)','first_Dmap (n_s=1,n_d=1)','detailed G','detailed D';Dstudy_G,Dstudy_D,first_iccG,first_iccD,detailed_icc_G,detailed_icc_D};

clockinfo{1}=datestr(now, 'HH'); clockinfo{2}=datestr(now, 'MM'); clockinfo{3}=datestr(now, 'SS');
timeID=sprintf('%sH%sM%sS',clockinfo{1},clockinfo{2},clockinfo{3});

if nvarcomp<=7
    
    exc_scale=18; % simple tweak to get desired colorscale
    figure;
    
    for i=1:2 % G-coeff plot, D-coeff plot
        subplot(1,2,i)
        image(((icc{2,i}>0.4)+(icc{2,i}>0.6)+(icc{2,i}>0.74))*exc_scale); %[0:3]
        ax=gca;
        
        if nvarcomp==3
            tickinc=find(n1_all==2)-1;
            
            ax.XTick=1:tickinc:length(n1_all);
            ax.XTickLabel=1:1:DSrange;
        elseif nvarcomp==7
            %         ax.XTick=linspace(1,size(icc{2,1},2),length(n1_all));
            %         ax.YTick=linspace(1,size(icc{2,1},1),length(n2_all));
            %         ax.XTickLabel=n1_all;
            %         ax.YTickLabel=n2_all;
            
            tickincX=find(n1_all==2)-1;
            tickincY=find(n2_all==2)-1;
            
            ax.XTick=1:tickincX:length(n1_all);
            ax.XTickLabel=1:1:DSrange;
            ax.YTick=1:tickincY:length(n2_all);
            ax.YTickLabel=1:1:DSrange;
        end
    end
    
    
    
    savefig(sprintf('tmp_DStudy_%s.fig',timeID));
    saveimg_squaresubplot({sprintf('tmp_DStudy_%s',timeID)},'no')
    close;
end

disp(sprintf('G: %0.2f+/-%0.2f',mean(first_iccG{1}),std(first_iccG{1})))
disp(sprintf('D: %0.2f+/-%0.2f',mean(first_iccD{1}),std(first_iccD{1})))

if(tempsave)
    save(sprintf('tmp_icc_%s.mat',timeID),'icc','var')
end


end
