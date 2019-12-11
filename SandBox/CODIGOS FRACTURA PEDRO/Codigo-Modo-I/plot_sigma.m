    
function plot_sigma(Ss,x_0,elem_0,ste_p,Area,d)

    rt=1;

    [Nodos,sp]=size(x_0);    
    [El,NNE]=size(elem_0);

    Sigma=zeros(El,ste_p);
    
    for i=1:El
        Sigma(i,:)=Ss(i*4-3,1:ste_p);
    end
    
    for nd=1:Nodos
        x(nd,1)=x_0(nd,1);
        y(nd,1)=x_0(nd,2);
    end

    DDD=max(x(:));
    HHH=max(y(:));

    PWnodo=zeros(Nodos,ste_p);
    x1=zeros(Nodos,1);
    y1=zeros(Nodos,1);
    % Nodal Pw2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    for nodo=1:Nodos
        for cont=1:ste_p
            suma1=0;
            suma2=0;
            for i=1:El
                for j=1:NNE
                    if nodo==elem_0(i,j)
                        suma1=suma1+Sigma(i,cont)*Area(i);
                        suma2=suma2+Area(i);
                    end
                end
            end;
            PWnodo(nodo,cont)=suma1/suma2;
        end;
    end;
    
    
    figure
    for cont=1:rt:ste_p
        for i=1:Nodos
            x1(i)=x(i)+d(i*sp-1,cont);
            y1(i)=y(i)+d(i*sp,cont);
        end
        pr=TriRep(elem_0,x1(:),y1(:),PWnodo(:,cont));
        trisurf(pr,'LineStyle','none')
        colorbar
        caxis([-30 30])
        axis([-2,DDD+2,0,HHH+2])
        drawnow
        pause;
    end;


end