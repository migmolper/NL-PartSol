
function EXPL_frac_ener(xg)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global INIT_file SAVE_I SAVE_F FB_BAR DOF
    global nodes sp elements NNE df x_0 elem elem_0 Area patch_con patch_el
    global g MAT boundary vad dad Material Mat_nds thickness tf ...
    ext_forces_s h_ini
    global time_step step_final
    global MASTER X_PS
    
    [~,SUP]=size(X_PS);
   
    [load_s,ds,as,vs,x_a,x_a1,ss,es,def_G,volume,status,...
    dens,jacobians,a,v,d,Ps,Es,Ss,Def_G,Status,W_beam,W_hammer,D_energy,...
    K_beam,K_hammer,React,Impact,tp,EP,JACO,NO,...
    React2,Impact2,L_beam,L_hammer,Sup_energy,W_elem]=init_frac;
    Inv=zeros(nodes*sp,nodes*sp);   % ???
    int_Fs_0=zeros(nodes*sp,1);
    h=h_ini.*sqrt(jacobians);
    
    E_ini(elements,1)=0;
    p_list=0;
    c_list=0;
    
    
    %--------------------------------------------------------------------------
    % Epsilon neighborhood & notch
    %--------------------------------------------------------------------------  
    global M1 H1
    h1=H1*0.05;
    Cr(1,:)=[M1 0 h1 inf];  % Initial notch
    [NO]=crack(xg,x_a,NO,Cr,1);
    for i=1:8
        Cr(i+1,:)=[M1 i*h1 h1 inf];  % Initial notch
        [NO]=crack(xg,x_a,NO,Cr,i+1);
    end
    
    [epsilon]=eps_nb(xg);
    

    %--------------------------------------------------------------------------
    % Initial state & initial shape functions
    %--------------------------------------------------------------------------

    if INIT_file==0 
        ste=1;
        t(1)=time_step;
        tp(1)=t(1);
    else
        load FILE ste ste_p xg t
    end
    t(step_final,1)=0;

    global TI_param
    gamma=TI_param(3);
    
    [B,near,p,gamma_,lam_LME,EP,FAIL]=LME_EP_ini(jacobians,volume,x_a,xg,NO);

    [lmass]=lumped_mass(near,p,volume,dens); %Initial mass
    [K_beam(1),K_hammer(1),L_beam(1),L_hammer(1)]=...
        initial_K(near,p,volume,dens,v(:,1)); %Initial kinetic energy
    d_energy=0;   %Initial dissipated energy
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    acont(nodes*sp,step_final)=0;
    asup(nodes*sp,step_final)=0;
    for ste=ste+1:step_final

        if ste>1
            t(ste)=t(ste-1)+time_step;
        end
        
%         if rem(ste,3000)==0
%             fprintf('ste %i \n',ste);
%         end
%         if rem(ste,100)==0
%             fprintf('ste %i \n',ste);
%         end

        load_s(:,2)=load_s(:,1);
        if t(ste)<tf
            load_s(:,1)=ext_forces_s*t(ste)/tf;
        else
            load_s(:,1)=ext_forces_s;
        end;

        % 1. Predictor
        C1=1;
        C2=1;
        iter=1;
        while C1 || C2
            for i=1:nodes
                for j=1:sp
                    % solid
                    if (boundary(i*sp+1-j)==0)
                        ds(i*sp+1-j,1)=ds(i*sp+1-j,2)+time_step*vs(i*sp+1-j,2)+0.5*time_step^2*as(i*sp+1-j,2);
                        vs(i*sp+1-j,1)=vs(i*sp+1-j,2)+(1-gamma)*time_step*as(i*sp+1-j,2);
                    else
                        if t(ste)<tf
                            ds(i*sp+1-j,1)=dad(i*sp+1-j)*(t(ste)/tf);
                        else
                            ds(i*sp+1-j,1)=dad(i*sp+1-j);
                        end
                        vs(i*sp+1-j,1)=(ds(i*sp+1-j,1)-ds(i*sp+1-j,2))/time_step;
                    end
                end
            end

            for j=1:nodes
                for i=1:sp
                    x_a1(j,i)=x_a(j,i)+(ds((j-1)*sp+i,1)-ds((j-1)*sp+i,2));
                end
            end

            if MASTER
                [as1,C1,c_list]=contact_2(x_a1,lmass,vs,as,c_list);
            else
                as1=as*0;
                C1=0;
            end
            
            acont(:,ste)=acont(:,ste)+as(:,2)-as1(:,2);
            as=as1;
             
            C2=0;
            if iter==1
                as1=as;
                for j=1:SUP
                    [as1,Cc,p_list]=PS_surface_2(x_a1,lmass,vs,as1,h,j,p_list);
                    if Cc
                        C2=1;
                    end
                end
            end
            asup(:,ste)=asup(:,ste)+as(:,2)-as1(:,2);
            as=as1;
            iter=iter+1;

        end

        x_a=x_a1;

        % 2. REMAPPING
%         REMAP=1;  %Flag
%         iter=1;
%         error_tol=min(h_ini.*sqrt(jacobians))*1e-6;
%         wrap=zeros(elements,1);
%         while REMAP==1
            [xg2]=update_mp(ds,near,p,xg);

            [def_G,jacobians,volume,dens,EP]=update_F(ds,near,B,def_G,volume,EP);
            h=h_ini.*sqrt(jacobians);
%             if iter==1
%                 xg1=xg2;
%             else
%                 error=norm(abs(xg1-xg2));
%                 if error<error_tol
%                     REMAP=0;
%                 else
%                     xg1=xg2;
%                     if iter>=3
%                         iter;
%                     end
%                 end
%             end
% 
%             if REMAP==1
%                 [B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
%                     volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);
%             end
%             iter=iter+1;
%         end
        xg=xg2;
        def_G(:,2)=def_G(:,1);


        % 3. Residual    
        [~,int_Fs,status,d_energy,ss,es,FAIL,Wlist,E_ini]=Constitutive_ft(3,def_G,...
            near,volume,B,jacobians,h,ss,es,status,epsilon,d_energy,E_ini,FAIL);

        [R,I]=reaction_2(int_Fs,c_list,p_list);
        
        % 4. Recompute lumped mass
        [lmass]=lumped_mass(near,p,volume,dens);

        %%%%  5.  SOLVER   %%%%
        %residual=-(int_Fs-int_Fs_0)+(load_s(:,1)-load_s(:,2));
        residual=-int_Fs+load_s(:,1);

        for i=1:nodes
            for j=1:sp
                if (boundary(i*sp+1-j)~=0)
                    lmass(i*sp+1-j,i*sp+1-j)=1;
                    Inv(i*sp+1-j,i*sp+1-j)=1;
                    residual(i*sp+1-j,1)=(vs(i*sp+1-j,1)-vs(i*sp+1-j,2))/time_step;
                else
                    Inv(i*sp+1-j,i*sp+1-j)=1/lmass(i*sp+1-j,i*sp+1-j);
                end
            end
        end

        %as(:,1)=Inv*residual(:,1)+as(:,2);
        as(:,1)=Inv*residual(:,1);

        % 6. Corrector
        for i=1:nodes
            for j=1:sp
                if (boundary(i*sp+1-j)==0)
                    vs(i*sp+1-j,1)=vs(i*sp+1-j,1)+gamma*time_step*as(i*sp+1-j,1);
                end
            end
        end

        % 7. Storage
        if rem(ste,SAVE_I)==0;

            if SAVE_I==1
                ste_p=ste
            else
                ste_p=ste/SAVE_I+1
            end
    
            d(:,ste_p)=ds(:,1);
            %a(:,ste_p)=as(:,1);
            %v(:,ste_p)=vs(:,1);
            
            % Impact and Reaction
            React(ste_p,1)=R;
            Impact(ste_p,1)=I; 
            
            %[L_beam,L_hammer,Impact2(ste_p,1),React2(ste_p,1)]=forces...
            %    (vs,dens,volume,near,p,ste_p,L_beam,L_hammer);
            
            % Fracture 
            %Status(:,ste_p)=status(:,1);
            D_energy(ste_p)=d_energy;
            [~,~,K_beam(ste_p),K_hammer(ste_p),...
                Sup_energy]= compute_energy(def_G,jacobians,volume,status,...
                near,p,dens,vs,d,Sup_energy,React,ste_p);
            %W_elem(:,ste_p)=Wlist;
            
            W_beam(ste_p)   = 0;
            W_hammer(ste_p) = 0;
            for e=1:elements
            	if Material(e)==1
                    if status(e,1)==0            
                        W_beam(ste_p) = W_beam(ste_p) + Wlist(e) * volume(e);
                    end
                else
                    W_hammer(ste_p) = W_hammer(ste_p) + Wlist(e) * volume(e);
                end
            end

            %JACO(:,ste_p)=jacobians(:);
            %FINT(:,ste_p)=int_Fs;
        
            % STRESS and STRAIN
%             for e=1:elements
%                 Ps(e,ste_p)=...                                  %PRESSURE
%                     (ss((e-1)*4+1,1)+ss((e-1)*4+2,1)+ss((e-1)*4+3,1))/3;                    
%             end

            %Es(:,ste_p)=es(:,1);                %STRAIN
            %Ss(:,ste_p)=ss(:,1);                %STRESS
            %Def_G(:,ste_p)=def_G(:,1);          %DEFORMATION GRADIENT

            tp(ste_p,1)=t(ste);
            
        end
        
        % 8. Save info
        if ((rem(ste/SAVE_I,SAVE_F)==0) || (ste==step_final) || (FAIL==1))
            %save FILE elem x_0 x_a xg Area ste d v a Ps t tp Es Status ...
            %    Ss Def_G near p ste_p Impact React EP dens JACO DOF ...
            %    K_beam K_hammer W_beam W_hammer D_energy Mat_nds ...
            %    React2 Impact2 L_beam L_hammer Sup_energy W_elem NO FINT
            save FILE K_beam K_hammer W_beam W_hammer D_energy ...
                ste_p Impact React d tp Sup_energy Status xg x_0
        end

        % 9. Update
        ds(:,2)=ds(:,1);
        vs(:,2)=vs(:,1);
        as(:,2)=as(:,1);
        
        ss(:,2)=ss(:,1);
        es(:,2)=es(:,1);
        
        %int_Fs_0=int_Fs;

    end
    
end

