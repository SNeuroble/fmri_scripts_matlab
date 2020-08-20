function [icc_summary,var_comp_mean] = calculate_masked_ICCs(icc_summary,stats,factor_tbl,mask)
% data is cell array, mask is logical
% stats from ICC calculation

mask=logical(mask);

iccs=icc_summary(3,1:2);
for i=1:length(iccs)
    for j=1:length(iccs{i});
        mean_iccs{i}{j}=nanmean(iccs{i}{j}(mask));
    end
end


do_stats=1;
if do_stats  % assuming 3-way ANOVA w pairwise interactions
    for i=1:7 % refvar, 2var, 3var, refx2, refx3, 2x3, residual
        var{i}=getvalsfromicccellarray(stats,i+1,5);
        var{i}(var{i}<0)=0;
    end
end

for i=1:7
    var_comp_mean(i)=nanmean(nanmean(var{i}(mask)));
end


factor_size=max(factor_tbl);
for ref_factor=1:3 % for each factor: subj, scanner, day
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
    
    icc_G_overall{ref_factor}=var_comp_mean(a)/(var_comp_mean(a)+var_comp_mean(d)/factor_size(b)+var_comp_mean(e)/factor_size(c)+var_comp_mean(g)/(factor_size(b)*factor_size(c)));
    icc_D_overall{ref_factor}=var_comp_mean(a)/(var_comp_mean(a)+var_comp_mean(b)/factor_size(b)+var_comp_mean(c)/factor_size(c)+var_comp_mean(d)/factor_size(b)+var_comp_mean(e)/factor_size(c)+var_comp_mean(f)/(factor_size(b)*factor_size(c))+var_comp_mean(g)/(factor_size(b)*factor_size(c)));
    
end

icc_summary = {'icc_G_mean','icc_G_overall','icc_D_mean','icc_D_overall'; mean_iccs{1},icc_G_overall,mean_iccs{2},icc_D_overall};

end