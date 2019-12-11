function [load_s,ds,as,vs,x_a,x_a1,ss,es,es_p,def_G,def_G_p,Eps_p,Sy,volume,...
    dens,jacobians,a,v,d,Ps,Es,Es_p,Ss,Def_G,Def_G_p,Gamma_tot,Sy_tot, ...
    Gamma_nds,React,tp,EP,JACO,...
    dw,aw,vw,pw,def_G_w,Pw,Def_G_w,n,perm,load_w,dt]=init
 

    global nodes sp elements df Area x_0
    global INIT_file DOF IMPLICIT
    global rho_w g MAT vad dad Material thickness tf
    global time_step dim

    %--------------------------------------------------------------------------
    % VECTORS OF PARAMETERS
    %--------------------------------------------------------------------------

    Sy=zeros(elements,2);
    dens(elements,1)=0;
    perm(elements,1)=0;
    n(elements,1)=0;
    for i=1:elements
         Sy(i,2)     = MAT(7,Material(i));
         if DOF==2
            n(i)      = MAT(16,Material(i));
            perm(i)   = MAT(15,Material(i))/rho_w/g;
            dens(i)   = n(i)*rho_w+(1-n(i))*MAT(3,Material(i));
         else
            dens(i)   = MAT(3,Material(i));
         end
    end
    
    
    %--------------------------------------------------------------------------
    % VECTORS OF ZEROS and ONES   or taken from FILE
    %--------------------------------------------------------------------------
    x_a1=zeros(nodes,2);                    % Predictive position
    
    if INIT_file==1
        load FILE d v a Ps Es Es_p Ss Def_G Def_G_p Gamma_tot Sy_tot ...
              Gamma_nds React tp EP JACO dens x_a ste_p
          
        jacobians=JACO(:,ste_p);
        volume=Area.*jacobians;
        
        load FILE Def_G_w Pw
        for e=1:elements
            n(e)=1-(1-MAT(16,Material(e)))/jacobians(e);
        end
        pw(:,2)=Pw(:,ste_p);
        def_G_w(:,2)=Def_G_w(:,ste_p);
 

        
        if DOF==2 && IMPLICIT==0
            for i=1:nodes
                dw(i*sp-1:i*sp,2)=d(i*df-1:i*df,ste_p);
                aw(i*sp-1:i*sp,2)=a(i*df-1:i*df,ste_p);
                vw(i*sp-1:i*sp,2)=v(i*df-1:i*df,ste_p);
                ds(i*sp-1:i*sp,2)=d(i*df-3:i*df-2,ste_p);
                as(i*sp-1:i*sp,2)=a(i*df-3:i*df-2,ste_p);
                vs(i*sp-1:i*sp,2)=v(i*df-3:i*df-2,ste_p);
            end
            dt=0;
        else
            aw=0;
            vw=0;
            as(:,2)=a(:,ste_p);
            vs(:,2)=v(:,ste_p);
            if DOF==2
                dt(:,2)=d(:,ste_p);
                for i=1:nodes
                    dw(i*sp-1:i*sp,2)=d(i*df-1:i*df,ste_p);
                    ds(i*sp-1:i*sp,2)=d(i*df-3:i*df-2,ste_p);
                end
            else
                ds(:,2)=d(:,ste_p);
                dw=0;
                dt=0;
            end
        end

        Eps_p(:,2)=Gamma_tot(:,ste_p);
        Sy(:,2)=Sy_tot(:,ste_p);

        es(:,2)=Es(:,ste_p);                %STRAIN
        es_p(:,2)=Es_p(:,ste_p);            %STRAIN pl
        ss(:,2)=Ss(:,ste_p);                %STRESS
        def_G(:,2)=Def_G(:,ste_p);          %DEFORMATION GRADIENT
        def_G_p(:,2)=Def_G_p(:,ste_p);      %DEFORMATION GRADIENT pl        
        
    else
        x_a=x_0;
        
        [def_G,def_G_w,def_G_p]=ini_F;
        jacobians=ones(elements,1);
        JACO=zeros(elements,dim);
        JACO(:,1)=jacobians;
        volume = Area * thickness;
        
        d(nodes*df,dim)=0;       %DISPLACEMENTS
        a(nodes*df,dim)=0;       %ACCELERATIONS
        v(nodes*df,dim)=0;       %VELOCITIES
        
        load_w = zeros(nodes*sp,2);
        pw     = zeros(elements,2);

        
        tp(dim,1)=0;
        React(dim,1)=0;
        
        Eps_p=zeros(elements,2); 
        load_s=zeros(nodes*sp,2);

        ss=zeros(4*elements,2);                  %STRESS
        es=zeros(4*elements,1);                  %STRAIN
        es_p=zeros(4*elements,1);                %STRAIN_P
        
        Ps(elements,dim)=0;                   %PRESSURE
        Es(4*elements,dim)=0;                 %STRAIN
        Es_p(4*elements,dim)=0;               %STRAIN_P
        Ss(4*elements,dim)=0;                 %STRESS
        Def_G(elements*5,dim)=0;              %DEFORMATION GRADIENT
        Def_G_p(elements*5,dim)=0;            %PLASTIC DEFORMATION GRADIENT

        Sy_tot(elements,dim)=0;
        Gamma_tot(elements,dim)=0;

        Gamma_nds(nodes,dim)=0;

        EP(elements*3,2)=0;                   %Principal strain
        
        
        Pw(elements,dim)=0;                   %PORE PRESSURE
        Def_G_w(elements*5,dim)=0;            %water DEFORMATION GRADIENT
        

        if tf>0
            d(:,1)=dad(:)*(time_step/tf);
            v(:,1)=vad(:)*(time_step/tf);
        else
            d(:,1)=dad(:);
            v(:,1)=vad(:);
        end
        
        if DOF==2 && IMPLICIT==0
            for i=1:nodes
                ds(i*sp-1:i*sp,2)=d(i*df-3:i*df-2,1);
                dw(i*sp-1:i*sp,2)=d(i*df-1:i*df,1);
                vs(i*sp-1:i*sp,2)=v(i*df-3:i*df-2,1);
                vw(i*sp-1:i*sp,2)=v(i*df-1:i*df,1);
            end
            aw     = zeros(nodes*sp,2);
            as     = zeros(nodes*sp,2);
            dt     = 0;
        else
            vw=0;
            aw=0;
            vs(:,2)=v(:,1);
            as(:,2)=a(:,1);
            if DOF==2
                dt(:,2)=d(:,1);
                for i=1:nodes
                    dw(i*sp-1:i*sp,2)=d(i*df-1:i*df,1);
                    ds(i*sp-1:i*sp,2)=d(i*df-3:i*df-2,1);
                end
            else
                dw=0;
                dt=0;
                ds(:,2)=d(:,1);
            end

        end
       
    end

end

function [def_G,def_G_w,def_G_p]=ini_F

    global elements

    def_G=zeros(elements*5,2);
    def_G_w=zeros(elements*5,2);
    def_G_p=zeros(elements*5,2);
    for e=1:elements       
        def_G((e-1)*5+1,1)=1; 
        def_G((e-1)*5+4,1)=1;
        def_G((e-1)*5+5,1)=1;
        
        def_G_p((e-1)*5+1,1)=1; 
        def_G_p((e-1)*5+4,1)=1;
        def_G_p((e-1)*5+5,1)=1;  
        
        def_G_w((e-1)*5+1,1)=1; 
        def_G_w((e-1)*5+4,1)=1;
        def_G_w((e-1)*5+5,1)=1;
    end
    def_G(:,2)=def_G(:,1);
    def_G_p(:,2)=def_G_p(:,1);
    def_G_w(:,2)=def_G_w(:,1);
end