
function EXPL_solver(xg)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global INIT_file SAVE_I SAVE_F FB_BAR DOF UW MASTER
    global nodes sp elements NNE df x_0 elem elem_0 Area patch_con patch_el
    global rho_w g K_w K_s MAT boundary vad dad Material Mat_nds thickness tf ...
    ext_forces_s ext_forces_w h_ini
    global time_step step_final
    
    [load_s,ds,as,vs,x_a,x_a1,ss,es,es_p,def_G,def_G_p,Eps_p,Sy,volume,...
    dens,jacobians,a,v,d,Ps,Es,Es_p,Ss,Def_G,Def_G_p,Gamma_tot,Sy_tot, ...
    Gamma_nds,React,tp,EP,JACO,...
    dw,aw,vw,pw,def_G_w,Pw,Def_G_w,n,perm,load_w]=init;

    Inv=zeros(nodes*sp,nodes*sp);   % ???
    int_Fs_0=zeros(nodes*sp,1);

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
    
    [B,near,p,gamma_,lam_LME,EP,FAIL]=LME_EP_ini(jacobians,volume,x_a,xg);

    %[lmass]=lumped_mass(near,p,volume,dens); %Initial mass

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ste=ste+1:step_final

        if ste>1
            t(ste)=t(ste-1)+time_step;
        end
        
        if rem(ste,500)==0
            fprintf('ste %i \n',ste);
        end

        load_s(:,2)=load_s(:,1);
        if t(ste)<tf
            load_s(:,1)=ext_forces_s*t(ste)/tf;
        else
            load_s(:,1)=ext_forces_s;
        end;

        % 1. Predictor
        C1=1;
        C2=1;
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

    %         if MASTER
    %             [a1,C1]=contact(...
    %                 nCC,MASTER,M_LIST,Mat_nds,x_a1,time_step,lmass,vs,as);
    %         else
    %             C1=0;
    %         end
    %         
    %         [a1,C2]=PS_surface(nPS,PS_LIST,x_a1,lmass,vs,as,x_ps,K,h_min);
              C1=0;
              C2=0;

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

        % 3. Recompute lumped mass
        [lmass]=lumped_mass(near,p,volume,dens);

        % 4. Residual
        [~,int_Fs,Eps_p,Sy,def_G_p,ss,es,es_p,FAIL]=Constitutive...
         (3,def_G,def_G_p,Eps_p,Sy,near,volume,B,jacobians,ss,es,es_p,FAIL);


        %%%%  5.  SOLVER   %%%%
        residual=-(int_Fs-int_Fs_0)+(load_s(:,1)-load_s(:,2));

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

        as(:,1)=Inv*residual(:,1)+as(:,2);


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

    %        React(ste_p,1)=R;
    
            d(:,ste_p)=ds(:,1);
            a(:,ste_p)=as(:,1);
            v(:,ste_p)=vs(:,1);
            
            Gamma_tot(:,ste_p)=Eps_p(:,1);
            Sy_tot(:,ste_p)=Sy(:,1);

            [Gamma_nds(:,ste_p)]=Ep2Ep_n(Gamma_tot,p,near,ste_p);

            for e=1:elements
                Ps(e,ste_p)=...                                  %PRESSURE
                    (ss((e-1)*4+1,1)+ss((e-1)*4+2,1)+ss((e-1)*4+3,1))/3;                    
            end

            JACO(:,ste_p)=jacobians(:);

            Es(:,ste_p)=es(:,1);                %STRAIN
            Es_p(:,ste_p)=es_p(:,1);            %STRAIN pl
            Ss(:,ste_p)=ss(:,1);                %STRESS
            Def_G(:,ste_p)=def_G(:,1);          %DEFORMATION GRADIENT
            Def_G_p(:,ste_p)=def_G_p(:,1);      %DEFORMATION GRADIENT pl

            tp(ste_p,1)=t(ste);
            

        end
        
        % 8. Save info
        if ((rem(ste/SAVE_I,SAVE_F)==0) || (ste==step_final) || (FAIL==1))
            save FILE elem x_0 x_a xg Area ste d v a Ps t tp Es Es_p ...
                Ss Def_G Def_G_p Gamma_tot Sy_tot near p ste_p ...
                React EP dens Gamma_nds JACO DOF
        end


        % 9. Update
        ds(:,2)=ds(:,1);
        vs(:,2)=vs(:,1);
        as(:,2)=as(:,1);
        
        ss(:,2)=ss(:,1);
        es(:,2)=es(:,1);
        es_p(:,2)=es_p(:,1);
        
        Eps_p(:,2)=Eps_p(:,1);
        Sy(:,2)=Sy(:,1);
        
        int_Fs_0=int_Fs;

    end
    
end

