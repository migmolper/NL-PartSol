function plot_f(elem,x_0,ste_p,d,DOF,hh,rr)

film=0; %film or not?
flag2=1; %Displacements of the whole
%hh=1;     %Amplification in pressure plot
%rr=1;   %Frequency in pressure plot

%%            
contador=ste_p;

[Elements,NNE]=size(elem);
[Nodos,sp]=size(x_0);
df=DOF*sp;
x=zeros(Nodos,1);
y=zeros(Nodos,1);

for nd=1:Nodos
    x(nd,1)=x_0(nd,1);
    y(nd,1)=x_0(nd,2);
end

DDD=max(x(:));
HHH=max(y(:));


if film
    movie = VideoWriter('video.avi');
    open(movie);
end


%%
if flag2==1
    x2=zeros(Nodos,1);
    y2=zeros(Nodos,1);
    figure
    for cont=1:rr:contador
        for nodo=1:Nodos
            if DOF==2
                x2(nodo)=x(nodo)+hh*d(df*nodo-3,cont);
                y2(nodo)=y(nodo)+hh*d(df*nodo-2,cont);
            else
                x2(nodo)=x(nodo)+hh*d(df*nodo-1,cont);
                y2(nodo)=y(nodo)+hh*d(df*nodo,cont);
            end
        end
        if NNE==3
            triplot(elem,x2,y2)
        elseif NNE==4
            quadplot(elem,x2,y2)
        end
        axis([-0.05,DDD*1.2,-0.05,HHH*1.1])
        drawnow
        
        if film
               frame = getframe;
               writeVideo(movie,frame);  
        end

        pause;
        if cont==42
            cont
        end
    end
    
    if film
        close(movie);
    end
end

end
