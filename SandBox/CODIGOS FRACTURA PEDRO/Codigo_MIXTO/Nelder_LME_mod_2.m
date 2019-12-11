function [p_a,J,cnt,lam,TOL]=Nelder_LME_mod_2(e,x,x_a,beta,near,TolLag,lam,ndim)
    %
    max1=250;
    min1=5;
    %
    rho=1;
    xi=2;
    gam=0.5;
    sig=0.5;
    epsilon=1e-3;

    pts=ndim+1;
    f(pts)=0;

    %Lambda


    L(1,:)=[lam(1),0];
    L(2,:)=[lam(1),lam(2)];
    L(3,:)=[0,lam(2)];

    for i=1:pts
        lam=L(i,:);
        [~,val,~,~]=Gamma_(ndim,x_a,x,beta,lam,near);
        f(i)=norm(val);
    end


    % 1- Ordenar

    [LL,f_m,f_mm,f_w]=order(L,f);
    
    %%%%% ROTATION ?? %%%%%
    ROT=1;
    iter=1;
    while ROT==1 && iter<100
        [Vn,diam]=von(LL);
        if Vn < epsilon
            ROT=1;
            [LL,f_m,f_mm,f_w]=rota2D(LL,near,x_a,x,beta,diam);
        else
            ROT=0;
        end
        iter=iter+1;
    end
    %%%%%%%%%%%%%%%
    
    cnt=1;
    tol(cnt)=10;
    while (tol(cnt)>TolLag && cnt<max1) || cnt<min1

        %2 - Reflexion
        sum=zeros(1,ndim);
        for i=1:ndim
            sum=sum+LL(i,:);
        end
        x_bar=sum/ndim;
        x_r=x_bar+rho*(x_bar-LL(3,:));

        [~,val,~,~]=Gamma_(ndim,x_a,x,beta,x_r,near);
        f_r=norm(val);

        if f_r<f_m   
            %3-Expansion
            x_e=x_bar+xi*(x_r-x_bar);
            [~,val,~,~]=Gamma_(ndim,x_a,x,beta,x_e,near);
            f_e=norm(val);

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
                [~,val,~,~]=Gamma_(ndim,x_a,x,beta,x_c,near);
                f_c=norm(val);
                if f_c<f_r
                    cont=1;
                end

            else
                %4- Contraccion dentro
                x_c=x_bar-gam*(x_bar-LL(3,:));
                [~,val,~,~]=Gamma_(ndim,x_a,x,beta,x_c,near);
                f_c=norm(val);
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
                    [~,val,~,~]=Gamma_(ndim,x_a,x,beta,lam,near);
                    f(i)=norm(val);
                end
                f(1)=f_m;
                [LL,f_m,f_mm,f_w]=order(LL,f);
            end
        end
        cnt=cnt+1;
        %%%%% ROTATION ?? %%%%%
        ROT=1;
        iter=1;
        while ROT==1 && iter<100
            [Vn,diam]=von(LL);
            if Vn < epsilon
                ROT=1;
                [LL,f_m,f_mm,f_w]=rota2D(LL,near,x_a,x,beta,diam);
            else
                ROT=0;
            end
            iter=iter+1;
        end
        %%%%%%%%%%%%%%%
        tol(cnt)=abs(f_w-f_m);
        if cnt>=230
            cnt;
        end   
    end    
    TOL=tol(cnt);
    lam=LL(1,:);
    [~,~,J,p_a]=Gamma_(ndim,x_a,x,beta,lam,near);
end


function [LL,f_m,f_mm,f_w]=order(L,f)
    [pts,~]=size(L);

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

function [Vn,diam]=von(LL)
    [n,sp]=size(LL);
    mat=ones(n);
    
    mat(:,1:2)=LL;
    DT=abs(det(mat)/factorial(sp));
    
    dist(n,n)=0;
    for i=1:n
        for k=1:n
            for j=1:sp
                dist(i,k)=dist(i,k)+(LL(i,j)-LL(k,j))^2;
            end
            dist(i,k)=sqrt(dist(i,k));
        end
    end
    diam=max(max(dist));
    
    Vn=DT/diam^sp;
end

function [LL2,f_m,f_mm,f_w]=rota2D(LL,near,x_a,x,beta,diam)

    [~,sp]=size(LL);

    DIS=max(diam/10,eps*2);
    
    
    LL1(1,:)=LL(1,:);
    LL2(1,:)=LL(1,:);
    
    LL1(2,1)=LL(2,1)+DIS;
    LL1(2,2)=LL(2,2);
    LL1(3,1)=LL(3,1);
    LL1(3,2)=LL(3,2)-DIS;
    
    LL2(2,1)=LL(2,1)-DIS;
    LL2(2,2)=LL(2,2);
    LL2(3,1)=LL(3,1);
    LL2(3,2)=LL(3,2)+DIS;
    
    f_1(1)=1e10;
    f_2(1)=1e10;
    [~,val,~,~]=Gamma_(sp,x_a,x,beta,LL1(2,:),near);
    f_1(2)=norm(val);
    [~,val,~,~]=Gamma_(sp,x_a,x,beta,LL1(3,:),near);
    f_1(3)=norm(val);
    [~,val,~,~]=Gamma_(sp,x_a,x,beta,LL2(2,:),near);
    f_2(2)=norm(val);
    [~,val,~,~]=Gamma_(sp,x_a,x,beta,LL2(3,:),near);
    f_2(3)=norm(val);
    
    f1=min(f_1);
    f2=min(f_2);
    
    if f1<f2
        LL2=LL1;
        f_2=f_1;
    end
    [~,val,~,~]=Gamma_(sp,x_a,x,beta,LL1(1,:),near);
    f_2(1)=norm(val);
    
    [LL2,f_m,f_mm,f_w]=order(LL2,f_2);
end

