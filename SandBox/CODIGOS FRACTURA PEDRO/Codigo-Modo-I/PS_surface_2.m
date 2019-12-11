
function[a,t,c_list]=PS_surface_2(X,lmass,v,a,h,SUP,c_list)

    global mu_2 K PS_LIST X_PS
    
    mu=mu_2;
    tol=min(h)*0.01;
    
    % LIST
    LIST=PS_LIST{SUP};
    x_ps=X_PS{SUP};
    nPS=length(x_ps)-1;

    %In or Out
    [nP,sp]=size(X);
    t=0;
    CC_m(1)=0;
    k=1;
    for i=1:nP
        [I,Cont]=IoO(X(i,:),nPS,LIST,i,x_ps,tol);

        if I
            
            ll=length(c_list);
            if ll==1 && c_list==0
                c_list(1)=i;
            else
                g=0;
                for j=1:ll
                    if i==c_list(j)
                        g=1;
                    end
                end
                if g==0
                    c_list(ll+1)=i;
                end
            end
            
            t=t+1;
            CONTACT(t,1)=i;        % Node
            CONTACT(t,2)=Cont(1);  % Contour
            CONTACT(t,3)=Cont(2);  % Delta

            
            list=LIST{CONTACT(t,2)};
            v2=x_ps(list(length(list)),2)-x_ps(list(1),2);
            v1=x_ps(list(length(list)),1)-x_ps(list(1),1);
            n_vec=[-v2 v1]/sqrt(v1^2+v2^2);
            VEC(t,:)=n_vec;
        
            l=0;
            for j=1:k-1
                if CONTACT(t,2)==CONTACT(CC_m(j),2);
                    l=1;
                end
            end
            if l==0
                CC_m(k)=t;
                k=k+1;
            end 
            
        end
    end
    
    if t>0
        mass=zeros(nP,1);
        for nod=1:nP
            mass(nod,1)=lmass(nod*sp,nod*sp);
        end
        [as_n,as_t]=acce_contact(t,VEC,CONTACT,X,mass,mu,v,K);
        a(:,2)=a(:,2)+as_n+as_t;
    end
     
    t;

end

function [I,Cont]=IoO(X,nCC,LIST,J,x_ps,tol)
    I=1;
    for i=1:nCC
        list=LIST{i};
        x0=x_ps(list(1),:);
        xf=x_ps(list(length(list)),:);
        M=(xf(2)-x0(2))/(xf(1)-x0(1));
        if isinf(M) || abs(M)>1e2
            if (xf(2)-x0(2))>0
                M=1;
            else
                M=-1;
            end
            d(i)=M*(X(1,1)-x0(1,1));
        else
            if (xf(1)-x0(1))>0
                M0=1;
            else
                M0=-1;
            end
            rM=sqrt(M^2+1);
            d(i)=M0*(M*X(1,1)-X(1,2)+x0(1,2)-M*x0(1,1))/rM;
        end
        
        if i~=1 && d(1)*d(i)<=0
            I=0;
        end
    end
    
    Cont=zeros(1,2);
    if I==1 
        [a,b]=min(abs(d));
        if a<tol
            I=0;
        else
            Cont(1,1)=b;
            Cont(1,2)=a;
        end
    end

end

function [as_n,as_t]=acce_contact(t,VEC,CONTACT,X,mass,mu,V,K)

    [nP,sp]=size(X);
    as_n=zeros(sp*nP,1);
    as_t=zeros(sp*nP,1);
    
    %%%% NORMAL ACCE
    for i=1:t
        nod=CONTACT(i,1);
        %%% Vector
        n_vec=VEC(i,:);
        %%% Acceleration
        A=K*CONTACT(i,3)/mass(nod);
        as_n(sp*nod-1:sp*nod)=A*n_vec;
    end
    
    %%% TANGENTIAL ACCELERATIONS
    for i=1:t
        nod=CONTACT(i,1);
        
        %%% Vector
        n_vec=VEC(i,:);
        t_vec=[n_vec(2) -n_vec(1)];
        
        %%% Normal forces - friction
        N = abs(mu * mass(nod) * n_vec * as_n(sp*nod-1:sp*nod));
        
        %%% Stick forces
        f1=t_vec*V(sp*nod-1:sp*nod,1);
        
        %%% Force
        if f1==0
            Fs=0;
        else
            Fs=-f1/abs(f1)*min(N,abs(f1));
        end
        
        % Slave Accelerations
        as_t(sp*nod-1:sp*nod)=t_vec*Fs/mass(nod);      
        
        clear n_vec t_vec;
    end
    t;

end