function plot_inst_nds(elem,x_0,cont,d,Ss,mem)

hh=1;     %Amplification in pressure plot

%%    
[larg2,~]=size(Ss);
[Nodos,sp]=size(x_0);

x=zeros(Nodos,1);
y=zeros(Nodos,1);

for nd=1:Nodos
    x(nd,1)=x_0(nd,1);
    y(nd,1)=x_0(nd,2);
end

DDD=max(x(:));
HHH=max(y(:));

%%
large=larg2/Nodos;
plo=large-mem;
for i=1:Nodos
    SSnodo(i,1)=Ss(i*large-plo,cont);
end

x1=zeros(Nodos,1);
y1=zeros(Nodos,1);
figure
%for cont=1:rr:contador
    for i=1:Nodos
        x1(i)=x(i)+hh*d(i*sp-1,cont);
        y1(i)=y(i)+hh*d(i*sp,cont);
    end
    pr=TriRep(elem,x1(:),y1(:),SSnodo(:,1));
    trisurf(pr)
    colorbar
    axis([-2,DDD+2,0,HHH+2])
    drawnow
%end;

end
