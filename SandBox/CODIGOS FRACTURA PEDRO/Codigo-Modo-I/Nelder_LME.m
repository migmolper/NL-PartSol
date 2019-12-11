function [p_a,J,cnt,lam]=Nelder_LME(x,x_a,beta,near,TolLag,lam,ndim)
    %
    max1=200;
    min1=5;
    %
    rho=1;
    xi=2;
    gam=0.5;
    sig=0.5;

    pts=ndim+1;

    %Lambda

    L(1,:)=[0,0];
    L(2,:)=[lam(1),lam(2)];
    L(3,:)=[10+2*lam(1),10+2*lam(2)];

    for i=1:pts
        lam=L(i,:);
        [f(i),~,~,~]=Gamma_(ndim,x_a,x,beta,lam,near);
    end


    % 1- Ordenar

    [LL,f_m,f_mm,f_w]=order(L,f);

    cnt=0;
    tol=10;
    while (tol>TolLag && cnt<max1) || cnt<min1

        %2 - Reflexion
        sum=zeros(1,ndim);
        for i=1:ndim
            sum=sum+LL(i,:);
        end
        x_bar=sum/ndim;
        x_r=x_bar+rho*(x_bar-LL(3,:));

        [f_r,~,~,~]=Gamma_(ndim,x_a,x,beta,x_r,near);

        if f_r<f_m   
            %3-Expansion
            x_e=x_r+xi*(x_r-x_bar);
            [f_e,~,~,~]=Gamma_(ndim,x_a,x,beta,x_e,near);

            if f_e<f_m
                LL(3,:)=LL(2,:);
                LL(2,:)=LL(1,:);
                LL(1,:)=x_e;
                f_w=f_mm;
                f_mm=f_m;
                f_m=f_e;
            else
                LL(3,:)=LL(2,:);
                LL(2,:)=LL(1,:);
                LL(1,:)=x_r;
                f_w=f_mm;
                f_mm=f_m;
                f_m=f_r; 
            end

        elseif f_r<f_mm
            LL(3,:)=LL(2,:);
            LL(2,:)=x_r;
            f_w=f_mm;
            f_mm=f_r;

        else
            cont=0;
            if f_r<f_w
                %4- Contraccion fuera
                x_c=x_bar+gam*(x_r-x_bar);
                [f_c,~,~,~]=Gamma_(ndim,x_a,x,beta,x_c,near);
                if f_c<f_r
                    cont=1;
                end

            else
                %4- Contraccion dentro
                x_c=x_bar-gam*(x_bar-LL(3,:));
                [f_c,~,~,~]=Gamma_(ndim,x_a,x,beta,x_c,near);
                if f_c<f_w
                    cont=1;
                end
            end

            if cont
                if f_c<f_m
                    LL(3,:)=LL(2,:);
                    LL(2,:)=LL(1,:);
                    LL(1,:)=x_c;
                    f_w=f_mm;
                    f_mm=f_m;
                    f_m=f_c;

                elseif f_c<f_mm
                    LL(3,:)=LL(2,:);
                    LL(2,:)=x_c;
                    f_w=f_mm;
                    f_mm=f_c;
                else
                    LL(3,:)=x_c;
                    f_w=f_c;
                end

            else
                %5- Encogimiento

                for i=1:pts
                    LL(i,:)=LL(1,:)+sig*(LL(i,:)-LL(1,:));
                    lam=LL(i,:);
                    [f(i),~,~,~]=Gamma_(ndim,x_a,x,beta,lam,near);
                end
                f(1)=f_m;
                [LL,f_m,f_mm,f_w]=order(LL,f);
            end
        end
        cnt=cnt+1;
        tol=abs(f_w-f_m);

    end
    lam=LL(1,:);
    [~,~,J,p_a]=Gamma_(ndim,x_a,x,beta,lam,near);
end


function [LL,f_m,f_mm,f_w]=order(L,f)
    [pts,sp]=size(L);

    [f_m, lo] = min(f);            
    [f_w, hi] = max(f);

    for i=1:pts
        if i~=lo && i~=hi
            mm=i;
        end
    end
    f_mm=f(mm);
    LL(1,:)=L(lo,:);
    LL(2,:)=L(mm,:);
    LL(3,:)=L(hi,:);
end
