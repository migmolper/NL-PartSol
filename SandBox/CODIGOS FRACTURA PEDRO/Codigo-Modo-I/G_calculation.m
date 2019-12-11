function [GT]=G_calculation(d1,a1,v1,Fint,Mm,Cm,load,load1)
   
    
    global TI_param time_step

    af=TI_param(1);
    am=TI_param(2);
	delta=TI_param(3);
    alpha=TI_param(4);
    theta=TI_param(5);
    
    A=(1-am)/alpha/time_step/time_step/theta/theta;
    B=(1-am)/alpha/time_step/theta;
    C=(1-am)/2/alpha-1;
    D=(1-af)*delta/alpha/time_step/theta;
    E=1-(1-af)*delta/alpha;
    F=(1-af)*(1-delta/2/alpha)*time_step*theta;
    
    G=theta*(1-af);
    H=(1-af);
   
    du=d1(:,1)-d1(:,2);

    GT= G*(load-load1)+load1 - H*Fint...
            -Mm*(A*du-B*v1(:,2)-C*a1(:,2))...
            -Cm*(D*du+E*v1(:,2)+F*a1(:,2));

end

