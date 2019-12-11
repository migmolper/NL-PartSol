function [R,I]=reaction(residual)

    global M1 H1 s1 s2 ws wi x_0

    R=0;
    I=0;

    [nodes,sp]=size(x_0);

    for i=1:nodes
        if x_0(i,2)==min(x_0(:,2))
            if (x_0(i,1) <= (s1+ws/2) && x_0(i,1) >= (s1-ws/2)) || ...
                    (x_0(i,1) <= (s2+ws/2) && x_0(i,1) >= (s2-ws/2))
                
                R=R+residual(sp*i);
            end
        elseif x_0(i,2)==H1
            if (x_0(i,1) <= (M1+wi/2) && x_0(i,1) >= (M1-wi/2))                
                I=I+residual(sp*i);
            end
            
        end
    end


end



