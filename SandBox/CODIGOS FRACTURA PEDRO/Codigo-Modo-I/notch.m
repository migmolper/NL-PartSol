
function [n_nds,n_mps]=notch(x,xg)


    [nodes,~]=size(x);
    [elements,~]=size(xg);

    M=210;
    H=48;

    n_nds(nodes,1)=0;
    n_mps(elements,1)=0;

    for e=1:elements
        if xg(e,2)<=H
            if xg(e,1)<M
                n_mps(e)=1;
            else
                n_mps(e)=2;
            end
        end   
    end

    for i=1:nodes
        if x(i,2)<=H
            if x(i,1)<M
                n_nds(i)=1;
            else
                n_nds(i)=2;
            end
        end   
    end



end