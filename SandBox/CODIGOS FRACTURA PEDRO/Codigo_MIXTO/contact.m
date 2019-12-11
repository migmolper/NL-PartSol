
function[a,t]=contact(X,lmass,v,a)

    global Mat_nds MASTER nCC M_LIST time_step mu
    delta_t=time_step;

    %In or Out
    [nP,sp]=size(X);
    t=0;
    CC_m(1)=0;
    k=1;
    for i=1:nP    
        
        if i==122
            i;
        end
        
        [I,Cont]=IoO(X(i,:),X,nCC,M_LIST,i);

        if I && Mat_nds(i)~=MASTER
            t=t+1;
            CONTACT(t,1)=i;        % Node
            CONTACT(t,2)=Cont(1);  % Contour
            CONTACT(t,3)=Cont(2);  % Delta

            
            list=M_LIST{CONTACT(t,2)};
            v2=X(list(length(list)),2)-X(list(1),2);
            v1=X(list(length(list)),1)-X(list(1),1);
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
        if k>2
            k;
        end
        [am_n,as_n,am_t,as_t]=acce_contact(...
            t,CC_m,VEC,CONTACT,M_LIST,X,delta_t,mass,mu,v);
        a(:,2)=a(:,2)+am_n+as_n+am_t+as_t;
    end
     
    t;

end

function [I,Cont]=IoO(X,x_a,nCC,LIST,J)
    I=1;
    for i=1:nCC
        list=LIST{i};
        x0=x_a(list(1),:);
        xf=x_a(list(length(list)),:);
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
        Cont(1,1)=b;
        Cont(1,2)=a;
        if isnan(a)
            a;
        end
    end

end

function [am_n,as_n,am_t,as_t]=acce_contact(...
                t,CC_m,VEC,CONTACT,LIST,X,delta_t,mass,mu,V)

    [nP,sp]=size(X);
    am_n=zeros(sp*nP,1);
    as_n=zeros(sp*nP,1);
    am_t=zeros(sp*nP,1);
    as_t=zeros(sp*nP,1);
    %%%% P FORCES
    P=zeros(t,1);
    for i=1:t
        list=LIST{CONTACT(i,2)};
        nod=CONTACT(i,1);
        %%% W slave
        dt=0;
        d1=zeros(length(list),1);
        for j=1:length(list)
            d1(j)=1/sqrt((X(list(j),1)-X(nod,1))^2+(X(list(j),2)-X(nod,2))^2);
            dt=dt+d1(j);
        end
        d1=d1/dt;
        Ws(i)={d1};
        clear list d1;
        %%% Forces
        P(i)=2*mass(nod)*CONTACT(i,3)/delta_t^2;
    end
    
    %%% NORMAL ACCELERATIONS
    for k=1:length(CC_m)
        list=LIST{CONTACT(CC_m(k),2)};
        n_vec=VEC(CC_m(k),:);
        for j=1:length(list)        
            %%% W master
            dt=0;
            d1(1,1)=0;
            nds(1,1)=0;
            l=0;
            for i=1:t
                if CONTACT(i,2)==CONTACT(CC_m(k),2)
                    l=l+1;
                    nod=CONTACT(i,1);
                    d1(l)=1/sqrt((X(list(j),1)-X(nod,1))^2+(X(list(j),2)-X(nod,2))^2);
                    dt=dt+d1(l);
                    nds(l)=i;
                end
            end
            d1=d1/dt;
            
            %%% Master accelerations
            a1=0;
            a2=mass(list(j));
            for i=1:l
                a1=a1+d1(i)*P(nds(i));
                a2=a2+d1(i)*mass(CONTACT(nds(i),1));
            end
            am_n(sp*list(j)-1:sp*list(j))=a1/a2*n_vec;
            
            Wm(j)={d1};
            Nm(j)={nds};
            
            clear d1 nds 
        end
        WC(k)={Wm};
        NC(k)={Nm};
        clear list Wm Nm
    end  
    for i=1:t
        list=LIST{CONTACT(i,2)};
        nod=CONTACT(i,1);
        %%% W slave
        ws=Ws{i};
        %%% Vector
        n_vec=VEC(i,:);
        %%% Slave accelerations
        a2=-P(i)/mass(nod);
        for j=1:length(list)
            a1=n_vec*am_n(sp*list(j)-1:sp*list(j));
            a2=a2+ws(j)*a1;
        end
        as_n(sp*nod-1:sp*nod)=a2*n_vec;
        clear listn ws n_vec;
    end
    
    %%% TANGENTIAL ACCELERATIONS
    for i=1:t
        list=LIST{CONTACT(i,2)};
        nod=CONTACT(i,1);
        
        %%% W slave
        ws=Ws{i};
        
        %%% Vector
        n_vec=VEC(i,:);
        t_vec=[n_vec(2) -n_vec(1)];
        
        %%% Normal forces - friction
        N = abs(mu * mass(nod) * n_vec * as_n(sp*nod-1:sp*nod));
        
        %%% Stick forces
        f1=t_vec*V(sp*nod-1:sp*nod,1);
        for j=1:length(list)
            f1=f1-t_vec*V(sp*list(j)-1:sp*list(j),1);
        end
        
        %%% Force
        Fs=-f1/abs(f1)*min(N,abs(f1));
        
        % Slave Accelerations
        as_t(sp*nod-1:sp*nod)=t_vec*Fs/mass(nod);      
        
        clear list ws n_vec t_vec;
    end
    % Master accelerations
    for k=1:length(CC_m)
        list=LIST{CONTACT(CC_m(k),2)};
        %%% W master
        Wm=WC{k};
        Nm=NC{k};   
        for j=1:length(list)
            wm=Wm{j};
            nds=Nm{j};
            for i=1:length(nds)
                nod=CONTACT(nds(i),1);
                am_t(sp*list(j)-1:sp*list(j))=am_t(sp*list(j)-1:sp*list(j))-...
                    as_t(sp*nod-1:sp*nod)*wm(i)*mass(nod)/mass(list(j));
            end
            clear wm nds
        end
        clear list Wm Nm
    end
    t;

end