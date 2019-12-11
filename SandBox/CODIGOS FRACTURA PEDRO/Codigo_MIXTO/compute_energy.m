
function [W_beam,W_hammer,K_beam,K_hammer,Sup_energy]=compute_energy...
    (def_G,jacobians,volume,status,near,p,dens,v1,d,Sup_energy,React,ste_p)

    global sp nodes MAT elements Material MAT_el x_0 AMP
    
    f_v=zeros(5,1);
    v2=zeros(nodes,1);
    W_beam   = 0;
    W_hammer = 0;
    K_beam=0;
    K_hammer=0;

    %% Velocity magnitude **********************
    for i=1:nodes
        for j=1:sp
            v2(i)=v2(i)+v1(i*sp+1-j,1)^2;
        end
        v2(i)=sqrt(v2(i));
    end

    %% Calculation of energy **********************
    for e=1:elements
        
        G  = MAT(4,MAT_el(e));
        Lam= MAT(5,MAT_el(e)); 

        %Build matrix
        for i=1:5
            f_v(i,1)=def_G((e-1)*5 + i,1);
        end           
        [F]=v2m(f_v,sp);
        C=F'*F;
        
        %Strain energy
        if Material(e)==1
            if status(e,1)==0            
                I1=trace(C);
                W=Lam/2*log(jacobians(e))^2 + G/2*(I1-3)-...
                  G*log(jacobians(e));
                W_beam = W_beam + W * volume(e);
            end
        else
            I1=trace(C);
            W=Lam/2*log(jacobians(e))^2 + G/2*(I1-3)-G*log(jacobians(e));
            W_hammer = W_hammer + W * volume(e);
        end
        
        %Kinetic energy
        m=length(near{e});
        nd = near{e};
        sh= p{e};
        
        K=0;
        for t1=1:m
        	K = K +dens(e)*volume(e)*sh(t1)*v2(nd(t1))*v2(nd(t1))/2;
        end
        
        if Material(e)==1
            K_beam = K_beam + K;
        else
            K_hammer = K_hammer + K;
        end
    end

    %% Calculation of support energy **********************
    D=0;
    D1=0;
    for i=1:nodes
        if x_0(i,2)==min(x_0(:,2))
            if (x_0(i,1) <= 72.0*AMP && x_0(i,1) >= 48.0*AMP) || ...
                (x_0(i,1) <= 372.0*AMP && x_0(i,1) >= 348.0*AMP)

                D=min(D,d(sp*i,ste_p));
                D1=min(D1,d(sp*i,ste_p-1));
            end
        end
    end
    Sup_energy(ste_p)=Sup_energy(ste_p-1)+(React(ste_p)+React(ste_p-1))/2*(D-D1);
end