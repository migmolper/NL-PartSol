function [R,I]=reaction_2(residual,c_list,p_list)

    global s1 s2 ws x_0 % M1 H1 wi 

    R=0;
    I=0;

    [nodes,sp]=size(x_0);

    for i=1:nodes
%         if x_0(i,2)==min(x_0(:,2))
%             if (x_0(i,1) <= (s1+ws/2) && x_0(i,1) >= (s1-ws/2)) || ...
%                     (x_0(i,1) <= (s2+ws/2) && x_0(i,1) >= (s2-ws/2))
%                 
%                 R=R+residual(sp*i);
%             end
% %         elseif x_0(i,2)==H1
% %             if (x_0(i,1) <= (M1+wi/2) && x_0(i,1) >= (M1-wi/2))                
% %                 I=I+residual(sp*i);
% %             end
%         else
            l1=length(p_list);
            if c_list~=0
                for j=1:l1
                    if i==p_list(j)
                        R=R+residual(sp*i);
                    end
                end
            end

            ll=length(c_list);
            if c_list~=0
                for j=1:ll
                    if i==c_list(j)
                        I=I+residual(sp*i);
                    end
                end
            end
    end


end



