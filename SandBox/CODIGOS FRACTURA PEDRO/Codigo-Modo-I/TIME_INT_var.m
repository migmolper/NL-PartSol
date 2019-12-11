function TIME_INT_var(IMPLICIT,TIS,DAMP)

    global step_final Time_final time_step NNE elements h_ini dim TI_param
    global SAVE_I MAT MAT_el Area SAVE_F
    
    af=0;
    am=0;
    delta=0;
    alpha=0;
    theta=0;

    tt(elements,1)=0;
    h(elements,1)=0;
    for e=1:elements
        if NNE==3
            h(e)=sqrt(2*Area(e));
        elseif NNE==4
            h(e)=sqrt(Area(e));
        end
        tt(e,1)=h(e)/MAT(6,MAT_el(e));
    end
    TT=min(tt);
    CFL=time_step/TT;
    fprintf('%f of CFL\n',CFL);
    h_ini=h;


    step_final= round(Time_final/time_step);


    if SAVE_I==1
        dim=floor(step_final/SAVE_I);
    else
        dim=floor(step_final/SAVE_I)+1;
    end
    fprintf('%i plot steps\n',dim);
    fprintf('Save %i times before the final\n',round(dim/SAVE_F));
        
    %% EXPLICIT
    if IMPLICIT==0
        
        delta=0.5;     %gamma
        
    else
    %% IMPLICIT
        if TIS==1                           %NEWMARK 1
            beta1=0.6;

            delta=beta1;
            
            theta=1;
            
        elseif TIS==2                       %NEWMARK
            if DAMP
                beta1=0.6;
                beta2=0.605;
            else
                beta1=0.5;
                beta2=0.5;
            end
            af=0;
            am=0;

            delta=beta1;
            alpha=beta2/2;

            theta=1;

        elseif TIS==3                   %GENERALIZED ALPHA
            rho=0.300;

            af=rho/(1+rho);
            am=(2*rho-1)/(1+rho);

            delta=0.5+af-am;
            alpha=0.25*(1-am+af)^2;

            theta=1;

        elseif TIS==4                  %HHT

            rho=0.800;
            am=0;
            af=(1-rho)/(1+rho);

            delta=(1+2*af)/2;
            alpha=0.25*(1+af)^2;

            theta=1;        

        elseif TIS==6                   %WBZ
            rho=0.800;

            af=0;
            am=(rho-1)/(1+rho);

            delta=0.5-am;
            alpha=0.25*(1-am)^2;

            theta=1;

        elseif TIS==5                  %WILSON-THETA

            theta=1.5;
            alpha=1/6;
            delta=1/2;

            af=0;
            am=0;

        elseif TIS==7                   %COLLOCATION METHOD

            theta=1.5;

            %a1=theta/2/(theta+1);
            %a2=(2*theta^2-1)/4/(2*theta^3-1);

            alpha=0.273;
            delta=1/2; 

            %alpha=0.305;
            %delta=0.6;

            af=0;
            am=0;
        end
    end
    TI_param=[af,am,delta,alpha,theta];



end