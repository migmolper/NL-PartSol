function [ext_forces_s,ext_forces_w]=ext_forces(nod,x_a,fd,Dir,UW,AGUA,Hf)

    if AGUA
        if Hf
            [ext_forces_s,ext_forces_w]=dist_f_x_w(nod,x_a,fd,Dir,UW);
        else
            [ext_forces_s,ext_forces_w]=dist_f_w(nod,x_a,fd,Dir,UW);
        end
    else
        ext_forces_w=0;
        if Hf
            [ext_forces_s]=dist_f_x(nod,x_a,fd,Dir);
        else
            [ext_forces_s]=dist_f(nod,x_a,fd,Dir);
        end
    end


end
function [ext_forces_s,ext_forces_w]=dist_f_x_w(nod,x_a,fd,Dir,UW)

    [nodes,sp]=size(x_a);
    [i,~]=size(nod);
    ext_forces_s=zeros(sp*nodes,1);
    ext_forces_w=zeros(sp*nodes,1);
    ft=0;

    maxy=1.0e-32;
    miny=1.0e32;
    for j=1:i
         if (x_a(nod(j),1)>maxy)
             maxy=x_a(nod(j),1);
         end
         if (x_a(nod(j),1)<miny)
             miny=x_a(nod(j),1);
         end
    end
    for j=1:i
        [d1,d2]=dist(x_a,nod,j);
        if (maxy==1.0e-32)
            d=0;
        elseif (x_a(nod(j),1)==maxy) || (x_a(nod(j),1)==miny)
            d=d1/2;
        else
            d=d1/2+d2/2;
        end
        f=fd*d;
        if Dir==1
            if UW==0 || UW==2
                ext_forces_s(nod(j)*sp-1)=f;
            elseif UW==1 || UW==2
                ext_forces_w(nod(j)*sp-1)=f;
            end
        elseif Dir==2
            if UW==0 || UW==2
                ext_forces_s(nod(j)*sp)=f;
            elseif UW==1 || UW==2
                ext_forces_w(nod(j)*sp)=f;
            end
        end
        ft=ft+f;
    end
    ft
end

function [ext_forces]=dist_f_x(nod,x_a,fd,Dir)

    [nodes,sp]=size(x_a);
    [i,~]=size(nod);
    ext_forces=zeros(sp*nodes,1);
    ft=0;

    maxy=1.0e-32;
    miny=1.0e32;
    for j=1:i
         if (x_a(nod(j),1)>maxy)
             maxy=x_a(nod(j),1);
         end
         if (x_a(nod(j),1)<miny)
             miny=x_a(nod(j),1);
         end
    end
    for j=1:i
        [d1,d2]=dist(x_a,nod,j);
        if (maxy==1.0e-32)
            d=0;
        elseif (x_a(nod(j),1)==maxy) || (x_a(nod(j),1)==miny)
            d=d1/2;
        else
            d=d1/2+d2/2;
        end
        f=fd*d;
        if Dir==1
            ext_forces(nod(j)*sp-1)=f;
        elseif Dir==2
            ext_forces(nod(j)*sp)=f;
        end
        ft=ft+f;
    end
    ft
end

function [ext_forces]=dist_f(nod,x_a,fd,Dir)

    [nodes,sp]=size(x_a);
    [i,~]=size(nod);
    ext_forces=zeros(sp*nodes,1);
    ft=0;

    maxy=1.0e-32;
    miny=1.0e32;
    for j=1:i
         if (x_a(nod(j),2)>maxy)
             maxy=x_a(nod(j),2);
         end
         if (x_a(nod(j),2)<miny)
             miny=x_a(nod(j),2);
         end
    end
    for j=1:i
        [d1,d2]=dist(x_a,nod,j);
        if (maxy==1.0e-32)
            d=0;
        elseif (x_a(nod(j),2)==maxy) || (x_a(nod(j),2)==miny)
            d=d1/2;
        else
            d=d1/2+d2/2;
        end
        f=fd*d;
        if Dir==1
            ext_forces(nod(j)*sp-1)=f;
        elseif Dir==2
            ext_forces(nod(j)*sp)=f;
        end
        ft=ft+f;
    end
    ft
end

function [ext_forces_s,ext_forces_w]=dist_f_w(nod,x_a,fd,Dir,UW)

    [nodes,sp]=size(x_a);
    [i,~]=size(nod);
    ext_forces_s=zeros(sp*nodes,1);
    ext_forces_w=zeros(sp*nodes,1);
    ft=0;

    maxy=1.0e-32;
    miny=1.0e32;
    for j=1:i
         if (x_a(nod(j),2)>maxy)
             maxy=x_a(nod(j),2);
         end
         if (x_a(nod(j),2)<miny)
             miny=x_a(nod(j),2);
         end
    end
    for j=1:i
        [d1,d2]=dist(x_a,nod,j);
        if (maxy==1.0e-32)
            d=0;
        elseif (x_a(nod(j),2)==maxy) || (x_a(nod(j),2)==miny)
            d=d1/2;
        else
            d=d1/2+d2/2;
        end
        f=fd*d;
        if Dir==1
            if UW==0 || UW==2
                ext_forces_s(nod(j)*sp-1)=f;
            elseif UW==1 || UW==2
                ext_forces_w(nod(j)*sp-1)=f;
            end
        elseif Dir==2
            if UW==0 || UW==2
                ext_forces_s(nod(j)*sp)=f;
            elseif UW==1 || UW==2
                ext_forces_w(nod(j)*sp)=f;
            end
        end
        ft=ft+f;
    end
    ft
end

function [d1,d2]=dist(x_a,nd,j)
    x1=x_a(nd(j),1);
    y1=x_a(nd(j),2);
    for i=1:length(nd)
        if(i~=j)
            d(i)=sqrt((x_a(nd(i),1)-x1)^2+(x_a(nd(i),2)-y1)^2);
        else
            d(i)=1.0e32;
        end
    end
    d1=min(d);
    t=0;
    i=0;
    while t==0&&i<length(nd)
        i=i+1;
        if d(i)==d1;
            t=1;
            d(i)=1.0e32;
        end  
    end
    d2=min(d);   
end