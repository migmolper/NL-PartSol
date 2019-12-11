function [FT]=G_forces(a1,v1,d1,K,Mm,Cm,load,load1)
   
    
    global TI_param time_step

    af=TI_param(1);
    am=TI_param(2);
	delta=TI_param(3);
    alpha=TI_param(4);
    theta=TI_param(5);
    
    %A=1/alpha/time_step/time_step/theta/theta;
    B=1/alpha/time_step/theta;
    C=1/2/alpha;
    %D=delta/alpha/time_step;
    E=delta/alpha;
    F=(1-delta/2/alpha)*time_step*theta;
    
%     if TIS==1
%         delta_R=zeros(dof,1);
%     else
        delta_R=load1-Mm*a1(:,2)-Cm*v1(:,2)-K*d1(:,2);
    %end

    FT= theta*(1-af)*(load-load1)+ delta_R...
            +(1-af)*Cm*(E*v1(:,2)-F*a1(:,2))...
            +(1-am)*Mm*(B*v1(:,2)+C*a1(:,2));

end

