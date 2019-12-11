function reaction_post

    load FILE_OK
    
    load EXP2640
    
    M1 = 210;       % Middle point  (mm)
    H1 = 100;       % Region Top (mm)
    s1 = 60;        % Middle support 1 (mm)
    s2 = 360;       % Middle support 2 (mm)
    ws = 20;        % Wide of the support (mm)
    wi = 2;        % Wide of the support (mm)


    [nodes,sp]=size(x_0);
    [~,cont]=size(FINT);
    
    velo=2640;
    del=0.1/velo;
    
    React_p(cont,1)=0;
    Impact_p(cont,1)=0;

    for k=1:cont
        for i=1:nodes
            if x_0(i,2)==min(x_0(:,2))
                if (x_0(i,1) <= (s1+ws/2) && x_0(i,1) >= (s1-ws/2)) || ...
                        (x_0(i,1) <= (s2+ws/2) && x_0(i,1) >= (s2-ws/2))

                    React_p(k)=React_p(k)+FINT(sp*i,k);
                end
            elseif x_0(i,2)==H1
                 if (x_0(i,1) <= (M1+wi/2) && x_0(i,1) >= (M1-wi/2))                
                     Impact_p(k)=Impact_p(k)+FINT(sp*i,k);
                 end
%             else
%                 ll=length(c_list);
%                 if ll>1 || c_list~=0
%                     g=0;
%                     for j=1:ll
%                         if i==c_list(j)
%                             g=1;
%                         end
%                     end
%                     if g==0
%                         Impact(k)=Impact(k)+residual(sp*i);
%                     end
%                 end
            end
        end
    end

    figure1 = figure;
    axes1 = axes('Parent',figure1);
    plot(tp(1:ste_p),React_p(1:ste_p)/1000,tp(1:ste_p),-Impact_p(1:ste_p)/1000,...
        EXP(:,1)+del,EXP(:,3),EXP(:,1)+del,EXP(:,2));
    xlabel('Time (s)');
    ylabel('Forces [kN]');
    legend('Reaction','Impact','EXP Reaction','EXP Impact');
    xlim(axes1,[0 0.0008]);

end



