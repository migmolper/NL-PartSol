function [p_lme,dp_lme,near_s_lme]=lme_wrapper...
    (target_zero,TolLag,gamma,blk,total,xsample,ndim,elem,Nelder,wrap,near,p)

    n_a=size(total,1);
    n_s=size(xsample,1);
    all_=[1:size(total,1)];

    beta_=zeros(n_s,1);
    for i=1:n_s
        beta_(i)=gamma/blk(i)^2;
    end

    %MAKE NEAR
    if target_zero<1e-11
        %t=cputime;
        for i=1:n_s
           near_s_lme(i)={all_};
        end
    else
        range1=zeros(n_s,1);
        range_lme=zeros(n_s,1);
        for i=1:n_s
            range1(i)=blk(i)*sqrt(-1/gamma*log(target_zero));
            range_lme(i)=max(range1(i)*1, 2*blk(i));
        end
        %t=cputime;
        [near_s_lme]=make_near(total,xsample,range_lme,all_,ndim,n_a,elem,wrap,near);
    end

    %MAKE SHAPE FUNCTION
    if ndim==1
        [p_lme,dp_lme]=shapef_1D...
            (ndim,total,beta_,xsample,near_s_lme,TolLag,elem,Nelder,wrap,p);
    elseif ndim==2
        [p_lme,dp_lme]=shapef...
            (ndim,total,beta_,xsample,near_s_lme,TolLag,elem,Nelder,wrap,p);
    end
    
end

