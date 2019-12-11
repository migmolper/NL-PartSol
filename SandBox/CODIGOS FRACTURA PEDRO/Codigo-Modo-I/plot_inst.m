function plot_inst(elem,x_0,elem_0,cont,d,Area,Ss,mem)

hh=1;     %Amplification in pressure plot

%%    
[larg2,~]=size(Ss);
[El,NNE]=size(elem_0);
[Nodos,sp]=size(x_0);

Sigma=zeros(El,1);

x=zeros(Nodos,1);
y=zeros(Nodos,1);

for nd=1:Nodos
    x(nd,1)=x_0(nd,1);
    y(nd,1)=x_0(nd,2);
end

DDD=max(x(:));
HHH=max(y(:));

%%
large=larg2/El;
plo=large-mem;
for i=1:El
    Sigma(i,1)=Ss(i*large-plo,cont);
end


PWnodo=zeros(Nodos,1);
x1=zeros(Nodos,1);
y1=zeros(Nodos,1);
% Nodal Pw2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for nodo=1:Nodos
    %for cont=1:contador
        suma1=0;
        suma2=0;
        for i=1:El
            for j=1:NNE
                if nodo==elem_0(i,j)
                    suma1=suma1+Sigma(i,1)*Area(i);
                    suma2=suma2+Area(i);
                end
            end
        end;
        PWnodo(nodo,1)=suma1/suma2;
    %end;
end;


figure
%for cont=1:rr:contador
    for i=1:Nodos
        x1(i)=x(i)+hh*d(i*sp-1,cont);
        y1(i)=y(i)+hh*d(i*sp,cont);
    end
    pr=TriRep(elem,x1(:),y1(:),PWnodo(:,1));
    trisurf(pr)
    colorbar
    axis([-2,DDD+2,0,HHH+2])
    drawnow
%end;

end
