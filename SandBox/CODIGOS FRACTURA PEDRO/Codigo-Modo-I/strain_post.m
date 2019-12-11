function strain_post(mult)

    load FILE
    
    [elements,~]=size(elem);
   
%------    STRAIN RATE of every M.P.    -------  
    
    E1_dot(elements,ste_p)=0;
    E2_dot(elements,ste_p)=0;
    
    for i=2:ste_p
       for j=1:elements
           E1_dot(j,i)=(Es(j*4-3,i)-Es(j*4-3,i-1))/(tp(i)-tp(i-1));
           E2_dot(j,i)=(Es(j*4-2,i)-Es(j*4-2,i-1))/(tp(i)-tp(i-1));           
       end       
    end
    

%------     STRAIN at strain gauges     -------

SG=[200 50;
    200 60;
    200 70;
    200 80] * mult;

D1=6*mult;
D2=2*mult;
tol=mult*0.5;

E_SG(4,ste_p)=0;
l(4,1)=0;
for e=1:elements
    for i=1:4
        if xg(e,1)<SG(i,1)+tol && xg(e,1)>SG(i,1)-D1-tol
            if xg(e,2)<SG(i,2)+tol && xg(e,2)>SG(i,2)-D2-tol
                l(i)=l(i)+1;
                for j=1:ste_p
                    E_SG(i,j)=E_SG(i,j)+Es(e*4-3,j);
                end
            end
        end
    end  
end
for i=1:4
     E_SG(i,:)=E_SG(i,:)/l(i);
     figure;plot(tp,E_SG(i,:));
end

%------    CRACK POSITION - VELOCITY    -------

    pos=50*ones(ste_p,1);
    vcrack(ste_p,1)=0;
    
    t=0;
    for i=2:ste_p
        p=50*mult;
        for j=1:elements
            if Status(j,i)==1 && xg(j,2)<96*mult
                p=max(p,xg(j,2));
            end
        end
        pos(i)=p;
        if t==0 && pos(i)~=50*mult
            t=i;
        end
        vcrack(i)=(pos(i)-pos(i-1))/(tp(i)-tp(i-1));
    end
    
    [xplot,yplot1,yplot2]=deriv(7,tp,t,ste_p,pos);
    
    vcrack2(:,1)=xplot;
    vcrack2(:,2)=yplot2;
    
    figure;plot(xplot,yplot1,tp(t-1:ste_p),pos(t-1:ste_p));
    figure;plot(xplot,yplot2,tp,vcrack);
    
    save STRAIN E1_dot E2_dot pos vcrack vcrack2 E_SG
    
end

function [xplot,yplot1,yplot2]=deriv(N,tp,t,ste_p,pos)

    fitResults1 = polyfit(tp(t-1:ste_p),pos(t-1:ste_p), N);
    
    ff=length(fitResults1);
    
    fitResults2=zeros(ff-1,1);
    for i=1:ff-1
        fitResults2(i)=fitResults1(i)*(ff-i);
    end
    
    xplot = linspace(tp(t-1),tp(ste_p));
    yplot1 = polyval(fitResults1, xplot);
    yplot2 = polyval(fitResults2, xplot);

end