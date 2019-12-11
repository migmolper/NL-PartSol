

function [K_beam,K_hammer,L_beam,L_hammer]=initial_K(near,p,volume,dens,v1)
    
    global nodes sp elements Material
    
    v2=zeros(nodes,1);
    v3=zeros(nodes,1);
    K_beam=0;
    K_hammer=0;
    L_beam=0;
    L_hammer=0;

    %% Velocity magnitude **********************
    for i=1:nodes
        for j=1:sp
            v2(i)=v2(i)+v1(i*sp+1-j,1)^2;
        end
        v2(i)=sqrt(v2(i));
        
        v3(i)=v1(i*sp,1);
    end

    %% Kinetic energy **********************
    for e=1:elements
        m=length(near{e});
        nd = near{e};
        sh= p{e};
        
        K=0;
        L=0;
        for t1=1:m
        	K = K +dens(e)*volume(e)*sh(t1)*v2(nd(t1))*v2(nd(t1))/2;
            L = L +dens(e)*volume(e)*sh(t1)*v3(nd(t1));
        end
        
        if Material(e)==1
            K_beam = K_beam + K;
            L_beam = L_beam + L;
        else
            K_hammer = K_hammer + K;
            L_hammer = L_hammer + L;
        end
    end
end