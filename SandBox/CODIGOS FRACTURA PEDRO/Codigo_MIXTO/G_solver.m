function [d1,a1,v1]=G_solver(InvK,FT,d1,a1,v1)

    global TI_param time_step

    %af=TI_param(1);
    %am=TI_param(2);
	delta=TI_param(3);
    alpha=TI_param(4);
    theta=TI_param(5);
    
    A=1/alpha/time_step/time_step/theta/theta;
    B=1/alpha/time_step/theta;
    C=1/2/alpha;
    
    incr_d_th=InvK*FT;
    incr_a_th=A*incr_d_th-B*v1(:,2)-C*a1(:,2);
    
%    incr_d=incr_d_th-(time_step*(theta-1)*v1+time_step^2/2*(theta^2-1)*a1...
%                      +alpha*time_step^2*(theta^2-1)*incr_a_th);
    
     a1(:,1)=a1(:,2)+incr_a_th/theta;
     v1(:,1)=v1(:,2)+time_step*a1(:,2)+delta*time_step/theta*incr_a_th;
     d1(:,1)=d1(:,2)+time_step*v1(:,2)+time_step^2/2*a1(:,2)+...
         alpha*time_step^2/theta*incr_a_th;
    
end