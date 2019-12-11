function [A,T,W,Ee]=Neo_Hookean(Kt,e,F,jacobians,w)
    
%     % Neo-Hookean Wriggers
%     A=0;
%     global MAT Material sp
%     
%     G  = MAT(4,Material(e));
%     Lam= MAT(5,Material(e)); 
%     J2=jacobians(e)^2;
%     
%     I=eye(sp+1);
%     b=F*F';
%     C=F'*F;
%     Ee=logm(C)/2;
% 
%     T=Lam/2*(J2-1)*I+G*(b-I);
%     T=T/jacobians(e);
%     
%     if Kt==1 || Kt==2      
%            A=[Lam+2*G    J2*Lam     0       J2*Lam ;
%                 J2*Lam    Lam+2*G   0       J2*Lam ;
%                 0          0   (1-J2)*Lam/2+G    0  ;
%                 J2*Lam   J2*Lam     0       Lam+2*G]; 
%     end
% 
% 
% 
% end

    % Neo-Hookean Bonet
    A=0;
    W=0;
    global MAT MAT_el sp
    
    G  = MAT(4,MAT_el(e));
    Lam= MAT(5,MAT_el(e)); 
    
    I=eye(sp+1);
    C=F'*F;
    b=F*F';
    Ee=logm(C)/2;
    
    
    T=Lam*log(jacobians(e))*I+G*(b-I);
    T=T/jacobians(e);
    %C_1=I/C;
    %T=Lam*log(jacobians(e))*C_1+G*(I-C_1);
    %T=F*T*F'/jacobians(e);
    
    if Kt                   
            lam_p=Lam/jacobians(e);
            mu_p=(G-Lam*log(jacobians(e)))/jacobians(e);
            
            A=[lam_p+2*mu_p    lam_p       0    lam_p;
                lam_p       lam_p+2*mu_p   0    lam_p;
                  0               0       mu_p     0;
                lam_p          lam_p       0     lam_p+2*mu_p];
    end

    if w
        I1=trace(C);
        W=Lam/2*log(jacobians(e))^2 + G/2*(I1-3)-G*log(jacobians(e));

    end


end

% 
%     % Neo-Hookean Ehlers
%     A=0;
%     global MAT Material sp
%     
%     G  = MAT( 4,Material(e));
%     Lam= MAT( 5,Material(e)); 
%     n0 = MAT(16,Material(e));
%     
%     I=eye(sp+1);
%     C=F'*F;
%     b=F*F';
%     Ee=logm(C)/2;
%     
%     T=Lam*n0^2*(jacobians(e)/n0-jacobians(e)/(jacobians(e)-1+n0))*I+G*(b-I);
%     T=T/jacobians(e);
%     
%     if Kt                   
%         lam_p=Lam/jacobians(e);
%         mu_p=(G-Lam*n0^2*(jacobians(e)/n0-jacobians(e)/(jacobians(e)-1+n0)))/jacobians(e);
% 
%         A=[lam_p+2*mu_p    lam_p       0    lam_p;
%             lam_p       lam_p+2*mu_p   0    lam_p;
%               0               0       mu_p     0;
%             lam_p          lam_p       0     lam_p+2*mu_p];
%     end
% 
% 
% 
% end