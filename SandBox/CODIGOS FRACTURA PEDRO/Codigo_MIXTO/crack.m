
function [NO]=crack(x_g,x_a,NO,Cr,A)

    global nodes elements

    d(nodes,1)=0;
    d2(nodes,1)=0;
    clas(nodes,1)=0;
    dg(elements,1)=0;
    dg2(elements,1)=0;
    clasg(elements,1)=0;

    h=Cr(A,3)*4;
    M=Cr(A,4);
    x0=Cr(A,1:2);

    for i=1:nodes
        % Identify side of the nodes
        if isinf(M) || abs(M)>1e2 
            d(i)=x_a(i,1)-x0(1,1);
        else
            rM=sqrt(M^2+1);
            d(i)=(M*x_a(i,1)-x_a(i,2)+x0(1,2)-M*x0(1,1))/rM;
        end
        d2(i)=sqrt((x_a(i,1)-x0(1,1))^2+(x_a(i,2)-x0(1,2))^2);
        if abs(d(i))>h/2 || d2(i)>h/2
            clas(i)=0;
        elseif d(i)>=0
            clas(i)=1;
        else
            clas(i)=2;
        end
    end

    for i=1:elements
        
        if i==438
            i;
        end
        % Identify side of the material point
        if isinf(M) || abs(M)>1e2 
            dg(i)=x_g(i,1)-x0(1,1);
        else
            rM=sqrt(M^2+1);
            dg(i)=(M*x_g(i,1)-x_g(i,2)+x0(1,2)-M*x0(1,1))/rM;
        end
        dg2(i)=sqrt((x_g(i,1)-x0(1,1))^2+(x_g(i,2)-x0(1,2))^2);
        if abs(dg(i))>h/2 || dg2(i)>h/2 
            clasg(i)=0;
        elseif dg(i)>=0
            clasg(i)=1;
        else
            clasg(i)=2;
        end

        % Initialize list of incompatible nodes
        no=NO{i};
        t=length(no);
        % Exclude nodes in the other side
        if clasg(i)~=0
            for j=1:nodes
                if clas(j)~=0 && clas(j)~=clasg(i)
                    l=0;
                    for k=1:t
                        if j==no(k)
                            l=1;
                        end
                    end
                    if l==0
                        t=t+1;
                        no(t)=j;
                    end
                end
            end
            NO(i)={no};
            %plot_nb(i,NO,x_a,x_g,elem,1,0);
            %A;
        end

        clear no
    end
    
end

