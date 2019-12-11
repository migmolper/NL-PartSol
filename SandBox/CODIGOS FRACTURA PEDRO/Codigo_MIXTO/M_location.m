
function [LIST,k]=M_location(x_a,x0,xc,CC,tol)

    global Mat_nds MASTER nodes

    [nCC,~]=size(CC);

    k=0;
    for i=1:nCC
        l=0;
        clear list
        if CC(i,1)==4
            R=CC(i,2);
            A=CC(i,3);
            [list2]=circ(xc,x0,x_a,R,A,tol);
            for j=1:length(list2)-1
                ll=[list2(j); list2(j+1)];
                k=k+1;
                LIST{k}=ll;
            end 
            x0=x_a(list2(length(list2)),:);
        else
            if CC(i,1)==1
                for j=1:nodes
                    d=x0(1)-x_a(j,1);
                    if abs(d)<tol && Mat_nds(j)==MASTER
                        l=l+1;
                        list(l)=j;
                    end               
                end 
            else

                M=CC(i,2);
                rM=sqrt(M^2+1);
                for j=1:nodes
                    d=(M*x_a(j,1)-x_a(j,2)+x0(1,2)-M*x0(1,1))/rM;
                    if abs(d)<tol && Mat_nds(j)==MASTER
                        l=l+1;
                        list(l)=j;
                    end               
                end       
            end 
            [list]=order(list,x0,x_a);
            x0=x_a(list(length(list)),:);
            k=k+1;
            LIST{k}=list;
        end
    end

end


function [list2]=circ(xc,x0,x_a,R,A,tol)

    [nodes,~]=size(x_a);
    
    %%%Select
    l=0;
    for i=1:nodes
        d=sqrt((x_a(i,1)-xc(1))^2+(x_a(i,2)-xc(2))^2);
        if abs(d)<R+tol  && abs(d)>R-tol
            l=l+1;
            list(l)=i;
        end 
    end
    
    %%%Order
    
    %%% 1st and 2nd
    for i=1:l
        if x0(1)==x_a(list(i),1) && x0(2)==x_a(list(i),2)
            list2(1)=list(i);
            list3(1)=i;
            break;
        end
    end
    [j,k]=dist2(x_a,list,i);
    
    V=[x_a(list2(1),1)-xc(1) x_a(list2(1),2)-xc(2)];
    
    P_vec=(x_a(list2(1),1)-xc(1))*(x_a(list(j),2)-xc(2))-...
        (x_a(list(j),1)-xc(1))*(x_a(list2(1),2)-xc(2));
    
    if P_vec>0
        list2(2)=list(j);
        list3(2)=j;
        V1=[x_a(list(j),1)-xc(1) x_a(list(j),2)-xc(2)];
    else
        list2(2)=list(k);
        list3(2)=k;
        V1=[x_a(list(k),1)-xc(1) x_a(list(k),2)-xc(2)];
    end
    
    i=2;
    Angle_old=acos((V(1)*V1(1)+V(2)*V1(2))/...
            sqrt(V(1)^2+V(2)^2)/sqrt(V1(1)^2+V1(2)^2));
    while i
        [j]=dist(x_a,list,list3(i),list3(i-1));
        V1=[x_a(list(j),1)-xc(1) x_a(list(j),2)-xc(2)];
        Angle=acos((V(1)*V1(1)+V(2)*V1(2))/...
            sqrt(V(1)^2+V(2)^2)/sqrt(V1(1)^2+V1(2)^2));
        %Exact angle
        if Angle<Angle_old
            if Angle_old>=pi
                Angle=2*pi-Angle;
            end
        end
        
        %Final point
        if Angle<A+tol && Angle>A-tol
            i=i+1;
            list2(i)=list(j);
            list3(i)=j;
            i=0;
        elseif Angle<A
            i=i+1;
            list2(i)=list(j);
            list3(i)=j;
            Angle_old=Angle;
        else
            i=0;
        end
    end    
end

function [i]=dist(x_a,nd,j,k)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    d=zeros(length(nd),1);
    for i=1:length(nd)
        if(i~=j) && (i~=k)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
        else
            d(i)=1.0e32;
        end
    end
    [~,i]=min(d);  
end

function [i,k]=dist2(x_a,nd,j)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    d=zeros(length(nd),1);
    for i=1:length(nd)
        if(i~=j)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
        else
            d(i)=1.0e32;
        end
    end
    [~,i]=min(d);
    d(i)=1.0e32;
    [~,k]=min(d);  
end

function [l2]=order(l1,x0,x_a)

    N=length(l1);

    d=zeros(N,1);
    for i=1:N
        d(i)=sqrt((x0(1,1)-x_a(l1(i),1))^2+(x0(1,2)-x_a(l1(i),2))^2);
    end

    l2=zeros(N,1);
    for i=1:N
        [~,b]=min(d);
        l2(i)=l1(b);
        d(b)=1e32;
    end
    N;
end

