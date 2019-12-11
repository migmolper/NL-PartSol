
function [epsilon]=eps_nb(xg)

    global elements h_ini MAT MAT_el W1 CR1 H1 H2 Material
    
    M1=CR1;
    
    h=h_ini;

    for e=1:elements
        
        Ceps=MAT(8,MAT_el(e));
        
        if xg(e,1)>(M1-W1/2) && xg(e,1)<(M1+W1/2) &&...
                xg(e,2)<H1 && xg(e,2)>H2
            
            C=Ceps*h(e);
            
            t=1;
            for i=1:elements
                d=sqrt((xg(e,1)-xg(i,1))^2+(xg(e,2)-xg(i,2))^2);
                if d<=C && i~=e && Material(e)==Material(i)
                    list(t)=i;
                    t=t+1;
                end
            end
            epsilon(e)={list};
            clear list;  
            
        else
            epsilon(e)={0};
        end
        
    end


end