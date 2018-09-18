function ordered_betas = categorize_betas(betas)

step=1/size(betas,1);
margin=0.01;

threshold=0.02;



for(is_pos=0:1:1) % is_pos==0 is negative

    ind=is_pos*2+1;
    
    if(is_pos==1)
        ordered_betas{1,ind}='Positive';
        ordered_betas{1,ind+1}='# Pos';
        temp=find(betas>threshold);
    else
        ordered_betas{1,ind}='Negative';
        ordered_betas{1,ind+1}='# Neg';
        temp=find(betas<(-threshold));
    end
    
    temp=temp/size(betas,1);
    temp2=temp-floor(temp); % could have just coded this by row... :(
    
    
    for(i=1:1:(size(betas,1)-1))
        ordered_betas{i+1,ind}=find(temp2<=((i)*step) & temp2>((i-1)*step));
        ordered_betas{i+1,ind}=temp(ordered_betas{i+1,ind});
        ordered_betas{i+1,ind+1}=size(ordered_betas{i+1,ind},1);
    end

    ordered_betas{((size(betas,1)+1)),ind}=find(temp2==0);
    ordered_betas{((size(betas,1)+1)),ind}=temp(ordered_betas{((size(betas,1)+1)),ind});
    ordered_betas{((size(betas,1)+1)),ind+1}=size(ordered_betas{((size(betas,1)+1)),ind},1);
end