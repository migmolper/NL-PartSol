function [A,T,W,Ee]=Saint_Venant(Kt,e,F)
    
    % St. Venant Material

    global MAT MAT_el sp
    
    G  = MAT(4,MAT_el(e));
    Lam= MAT(5,MAT_el(e));   
    I=eye(sp+1);

    A=0;
    W=0;
    
    C=F'*F;
    Ee=logm(C)/2;
    %Ee(2,1)=Ee(2,1)/2;
    %Ee(1,2)=Ee(1,2)/2;

    T=Lam*trace(Ee)*I+2*G*Ee;
    
    if Kt==1 || Kt==2
        
           A=[Lam+2*G    Lam     0   Lam ;
                Lam    Lam+2*G   0   Lam ;
                0          0     G    0  ;
                Lam      Lam     0   Lam+2*G]; 
            
    end
    
    if w
        I1=trace(Ee);
        sum=0;
        for i=1:3
            for j=1:3
                sum=sum+Ee(i,j)*Ee(i,j);
            end
        end
        W=Lam/2*I1^2 + G*sum;
        
    end
    
%     E_vec(1,1)=Ee(1,1);
%     E_vec(2,1)=Ee(2,2);
%     E_vec(4,1)=Ee(3,3);
%     E_vec(3,1)=Ee(2,1);
%     T2=A*E_vec;
%     e;



end