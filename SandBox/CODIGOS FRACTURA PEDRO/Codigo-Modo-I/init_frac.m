function [load_s,ds,as,vs,x_a,x_a1,ss,es,def_G,volume,status,...
    dens,jacobians,a,v,d,Ps,Es,Ss,Def_G,Status,W_beam,W_hammer,D_energy,...
    K_beam,K_hammer,React,Impact,tp,EP,JACO,NO,...
    React2,Impact2,L_beam,L_hammer,Sup_energy,W_elem]=init_frac
 

    global nodes sp elements df Area x_0
    global INIT_file %DOF IMPLICIT g
    global MAT vad dad MAT_el thickness tf
    global time_step dim

    %--------------------------------------------------------------------------
    % VECTORS OF PARAMETERS
    %--------------------------------------------------------------------------

    dens(elements,1)=0;
    for i=1:elements
    	dens(i) = MAT(3,MAT_el(i));
    end
    
    
    %--------------------------------------------------------------------------
    % VECTORS OF ZEROS and ONES   or taken from FILE
    %--------------------------------------------------------------------------
    x_a1=zeros(nodes,2);                    % Predictive position
    
    if INIT_file==1
        load FILE d v a Ps Es Ss Def_G W_beam W_hammer D_energy ...
              K_beam K_hammer Status React Impact tp EP JACO dens x_a ste_p ...
              React2 Impact2 L_beam L_hammer Sup_energy W_elem NO
          
        jacobians=JACO(:,ste_p);
        volume=Area.*jacobians*thickness;
        
        as(:,2)=a(:,ste_p);
        vs(:,2)=v(:,ste_p);
        ds(:,2)=d(:,ste_p);

        status(:,2)=Status(:,ste_p);

        es(:,2)=Es(:,ste_p);                %STRAIN
        ss(:,2)=Ss(:,ste_p);                %STRESS
        def_G(:,2)=Def_G(:,ste_p);          %DEFORMATION GRADIENT     
        
    else
        x_a=x_0;
        
        [def_G,~,~]=ini_F;
        jacobians=ones(elements,1);
        JACO=zeros(elements,dim);
        JACO(:,1)=jacobians;
        volume = Area * thickness;
        
        d(nodes*df,dim)=0;       %DISPLACEMENTS
        a(nodes*df,dim)=0;       %ACCELERATIONS
        v(nodes*df,dim)=0;       %VELOCITIES
               
        tp(dim,1)=0;
        React(dim,1)=0;
        Impact(dim,1)=0;
        
        status=zeros(elements,2); 
        load_s=zeros(nodes*sp,2);

        ss=zeros(4*elements,2);                  %STRESS
        es=zeros(4*elements,1);                  %STRAIN
        
        Ps(elements,dim)=0;                   %PRESSURE
        Es(4*elements,dim)=0;                 %STRAIN
        Ss(4*elements,dim)=0;                 %STRESS
        Def_G(elements*5,dim)=0;              %DEFORMATION GRADIENT

        Status(elements,dim)=0;
        W_elem(elements,dim)=0;
        
        W_beam(dim,1)=0;
        W_hammer(dim,1)=0;
        K_beam(dim,1)=0;
        K_hammer(dim,1)=0;  
        D_energy(dim,1)=0;
        Sup_energy(dim,1)=0;

        React2(dim,1)=0;
        Impact2(dim,1)=0;
        L_beam(dim,1)=0;
        L_hammer(dim,1)=0;


        EP(elements*3,2)=0;                   %Principal strain     
               

        if tf>0
            d(:,1)=dad(:)*(time_step/tf);
            v(:,1)=vad(:)*(time_step/tf);
        else
            d(:,1)=dad(:);
            v(:,1)=vad(:);
        end
        
        vs(:,2)=v(:,1);
        as(:,2)=a(:,1);
        ds(:,2)=d(:,1);
        
        NO(elements,1)={0};
        %for i=1:elements
        %    NO(i)={0};
        %end
       
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