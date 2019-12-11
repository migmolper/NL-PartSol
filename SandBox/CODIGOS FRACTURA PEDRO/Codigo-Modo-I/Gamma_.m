function [gam,dgam,hgam,p_a]=Gamma_(ndim,x_a,x,beta,lam,near)

if ndim==1
    temp=exp(-beta*((x-x_a(near)).^2) ...
        + lam*(x-x_a(near)) );
    Z=sum(temp);
    p_a=temp/Z;
    gam=log(Z);
    temp1(1:length(near))=(x-x_a(near)).*p_a;
    dgam=sum(temp1(:));
    hgam = sum( p_a.*(x-x_a(near)).^2   ) - dgam*dgam;
else
    sum1=0;
    sum2=0;
    for id=1:ndim
        sum1=sum1 + (x(id)-x_a(near,id)).^2;
        sum2=sum2 + lam(id)*(x(id)-x_a(near,id));
    end

    temp=exp(-beta*sum1 + sum2);
    Z=sum(temp);
    p_a=temp/Z;
    gam=log(Z);

    for id=1:ndim
        dgam(id)=sum( (x(id)-x_a(near,id)).*p_a(:) );
    end
    for id=1:ndim
        for jd=1:ndim
            hgam(id,jd)=sum( p_a(:).* ...
                (x(id)-x_a(near,id)).*(x(jd)-x_a(near,jd)) )  ...
                - dgam(id)*dgam(jd);
        end
    end
end