

function plot_mp(x_0,d,ste_p,Mat_nds,Status,xg)

    rr=1;

    [Elements,~]=size(Status);
    [DOF,contador]=size(d);
    contador=min(contador,ste_p);
    Nodos=DOF/2;
    sp=2;
    
    DDD=max(x_0(:,1));
    HHH=max(x_0(:,2));

    %
    for cont=1:rr:contador
        s1=0;
        s=0;
        for nodo=1:Nodos
            if Mat_nds(nodo)==1
                s=s+1;
                x2(s,1)=x_0(nodo,1)+d(sp*nodo-1,cont);
                y2(s,1)=x_0(nodo,2)+d(sp*nodo,cont);
            else
                s1=s1+1;
                x1(s1,1)=x_0(nodo,1)+d(sp*nodo-1,cont);
                y1(s1,1)=x_0(nodo,2)+d(sp*nodo,cont);
            end
        end
        t=0;
        for e=1:Elements
            if Status(e,cont)==1
                t=t+1;
                x3(t,1)=xg(e,1);
                y3(t,1)=xg(e,2);
            end
        end
        
        scatter(x2,y2,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]), ...
        hold on,...
        scatter(x1,y1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]), ...
        if t
            scatter(x3,y3,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'Marker','.'),...
        end
        hold off
        axis([-20,DDD+20,0,HHH/2])
        drawnow
        
        clear x2 y2 y1 x1 x3 y3
    end
    

end
