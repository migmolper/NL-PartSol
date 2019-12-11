

function plot_nb(NODO,near,x_0,xg,elem,ds,h)

[~,NNE]=size(elem);
[nodes,sp]=size(x_0);
x_a=x_0;
if h~=0
    for i=1:nodes
        x_a(i,1)=x_a(i,1)+h*ds(i*sp-1,1);
        x_a(i,2)=x_a(i,2)+h*ds(i*sp,1);
    end
end


if NODO
    nds=near{NODO};
    if nds(1)==0
        XX=zeros(length(nds)-1,2);
        for i=2:length(nds)
            XX(i,1)=x_a(nds(i),1);
            XX(i,2)=x_a(nds(i),2);
        end
    else
        XX=zeros(length(nds),2);
        for i=1:length(nds)
            XX(i,1)=x_a(nds(i),1);
            XX(i,2)=x_a(nds(i),2);
        end
    end
    if NNE==3
    hold on, triplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        scatter(xg(NODO,1),xg(NODO,2),'MarkerFaceColor',[1 1 0]), ...
        scatter(XX(:,1),XX(:,2),'MarkerFaceColor',[1 0 0]), hold off
    elseif NNE==4
    hold on, quadplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        scatter(xg(NODO,1),xg(NODO,2),'MarkerFaceColor',[1 1 0]), ...
        scatter(XX(:,1),XX(:,2),'MarkerFaceColor',[1 0 0]), hold off
    end
else
    if NNE==3
    hold on, triplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        hold off
    elseif NNE==4
    hold on, quadplot(elem,x_a(:,1),x_a(:,2)), scatter(xg(:,1),xg(:,2)), ...
        hold off
    end
    xlabel('X'),ylabel('Y')
end
    

end
