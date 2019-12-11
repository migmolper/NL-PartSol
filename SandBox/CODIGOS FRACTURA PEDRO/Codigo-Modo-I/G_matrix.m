
function [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx)

    global TI_param time_step

    af=TI_param(1);
    am=TI_param(2);
	delta=TI_param(3);
    alpha=TI_param(4);  
    theta=TI_param(5);
    
    A=(1-am)/alpha/time_step/time_step/theta/theta;
    %B=(1-am)/alpha/time_step/theta;
    %C=(1-am)/2/alpha;
    D=(1-af)*delta/alpha/time_step/theta;
    %E=1-(1-af)*delta/alpha;
    %F=(1-af)*(1-delta/2/alpha)*time_step*theta;

    
    matrix = (1-af)*stiff_mtx + A*mass_mtx + D*damp_mtx;
    
end


% delta=beta1=gamma
% alpha=beta2=beta