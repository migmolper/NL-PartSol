function [L_beam,L_hammer,I,R]=...
    forces(v1,dens,volume,near,p,ste_p,L_beam,L_hammer)

    global elements nodes SAVE_I time_step sp Material
    v2=zeros(nodes,1);
    
    %% Velocity magnitude **********************
    for i=1:nodes
        v2(i)=v1(i*sp,1);
    end
    
%     for i=1:nodes
%         for j=1:sp
%             v2(i)=v2(i)+v1(i*sp+1-j,1)^2;
%         end
%         v2(i)=sqrt(v2(i));
%     end
    
    for e=1:elements
        %Kinetic energy
        m=length(near{e});
        nd = near{e};
        sh= p{e};
        
        L=0;
        for t1=1:m
        	L = L +dens(e)*volume(e)*sh(t1)*v2(nd(t1));
        end
        
        if Material(e)==1
            L_beam(ste_p) = L_beam(ste_p) + L;
        else
            L_hammer(ste_p) = L_hammer(ste_p) + L;
        end
    end
    
    I=(L_hammer(ste_p)-L_hammer(ste_p-1))/time_step/SAVE_I;
    R=I+(L_beam(ste_p)-L_beam(ste_p-1))/time_step/SAVE_I;
    
end