
function Implicit_def(xg,NR)

    %--------------------------------------------------------------------------
    % Initialize variables
    %--------------------------------------------------------------------------
    global INIT_file SAVE_I SAVE_F DOF 
    global nodes sp elements df x_0 elem Area 
    global boundary dad tf ext_forces_s
    global time_step step_final
    
    [load_s,ds,as,vs,x_a,~,ss,es,es_p,def_G,def_G_p,Eps_p,Sy,volume,...
    dens,jacobians,a,v,d,Ps,Es,Es_p,Ss,Def_G,Def_G_p,Gamma_tot,Sy_tot, ...
    Gamma_nds,React,tp,EP,JACO,...
    ~,~,~,~,~,~,~,~,~,~]=init;

    %--------------------------------------------------------------------------
    % Initial state & initial shape functions and matrixes
    %--------------------------------------------------------------------------

    if INIT_file==0 
        t=zeros(step_final,1);
        ste=1;
        t(1)=time_step;
        tp(1)=t(1);
    else
        load FILE ste ste_p xg t
    end

    I=eye(df*nodes);
    
    [B,near,p,~,~,~,FAIL]=LME_EP_ini(jacobians,volume,x_a,xg);
    
    [mass_mtx]=mass_mat(near,p,dens,volume);
    damp_mtx=zeros(df*nodes);
    [stiff_mtx,int_Fs,~,~,~,~,~,~,FAIL]=Constitutive...
        (1,def_G,def_G_p,Eps_p,Sy,near,volume,B,jacobians,...
        ss,es,es_p,FAIL);

    [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx);
    matrix_1=matrix;
    
    for i=1:nodes*df
        if (boundary(i)~=0)
            for j=1:nodes*df
                matrix(i,j)=0;
                matrix(j,i)=0;
            end
            matrix(i,i)=1;
        end
    end
    
    InvK=matrix\I;
    %tic; InvK=matrix\I; t1=toc
    %tic; InvK_2=lusolve(matrix,0); t2=toc
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ste=ste+1:step_final

        if ste>1
            t(ste)=t(ste-1)+time_step;
        end

        if rem(ste,1000)==0
            fprintf('ste %i \n',ste);
        end

        % 1. Forces and imposed displacements
        load_s(:,2)=load_s(:,1);
        if t(ste)<tf
            load_s(:,1)=ext_forces_s*t(ste)/tf;
        else
            load_s(:,1)=ext_forces_s;
        end
        
        if t(ste)<tf
            i_disp=dad*time_step/tf;
        else
            i_disp=zeros(df*nodes);
        end
              
        [GT]=G_calculation(ds,as,vs,int_Fs,mass_mtx,...
            damp_mtx,load_s(:,1),load_s(:,2));     

        for i=1:nodes*df
            if (boundary(i)~=0)               
                GT(i)=i_disp(i);
            else
                for j=1:nodes*df
                    if (boundary(j)~=0)
                        GT(i)=GT(i)-matrix_1(i,j)*i_disp(j);
                    end
                end                     
            end
        end

        % --------------------------------------------------------
        % 2. Newton-Raphson
        clear iter error du
        tolerance=5e-5;
        TOL=tolerance*norm(GT);
        error(1)=1e32;
        iter=1;
        x_a1=x_a;

        NR1=NR;
        while error(iter) > TOL

            % 2.1 Solver        
            [du(:,iter)]=G_solver_1(InvK,GT,as,vs);
            ds(:,1)=ds(:,1)+du(:,iter);
             for j=1:nodes
                 for i=1:sp
                     x_a1(j,i)=x_a(j,i)+du((j-1)*sp+i,iter);
                 end
             end

            % 2.2 NO REMAPPING
            [def_G,jacobians,volume,dens,EP]=update_F(ds,near,B,def_G,volume,EP);

    %         [B,near,p,gamma_,lam_LME,REMAP,wrap,EP]=LME_EP(jacobians,...
    %             volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP,ste);

    
            % 2.3. Internal forces &/o stiffness matrix
            if rem(iter,NR1)==0
                
                [stiff_mtx,int_Fs,~,~,~,~,~,~,FAIL]=Constitutive...
                (1,def_G,def_G_p,Eps_p,Sy,near,volume,B,jacobians,...
                ss,es,es_p,FAIL);

                [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx);
                %matrix_1=matrix;

                for i=1:nodes*df
                    if (boundary(i)~=0)
                        for j=1:nodes*df
                            matrix(i,j)=0;
                            matrix(j,i)=0;
                        end
                        matrix(i,i)=1;
                    end
                end
                InvK=matrix\I;
            else
                [~,int_Fs,~,~,~,~,~,~,FAIL]=Constitutive...
                (0,def_G,def_G_p,Eps_p,Sy,near,volume,B,jacobians,...
                ss,es,es_p,FAIL);
            end

            % 2.4 G calculation
            [GT]=G_calculation(ds,as,vs,int_Fs,mass_mtx,...
                damp_mtx,load_s(:,1),load_s(:,2));

            for i=1:nodes*df
                if (boundary(i)~=0)               
                    GT(i)=0;%i_disp(i);
                else
                    for j=1:nodes*df
                        if (boundary(j)~=0)
                            GT(i)=GT(i);%-matrix_1(i,j)*i_disp(j);
                        end
                    end                     
                end
            end
            
            %
            iter=iter+1;
            if iter>100
                fprintf('iter %i \n',iter);
                if NR1>1
                    NR1=floor(NR1/2);
                end
            end
            error(iter)=norm(GT);
            %nu(iter)=norm(du(:,iter-1));
        
        end
        % --------------------------------------------------------

        % 3. Update variables
        def_G(:,2)=def_G(:,1);
        def_G_p(:,2)=def_G_p(:,1);
        x_a=x_a1;
        [as,vs]=G_solver_2(ds,as,vs);
        [xg]=update_mp(ds,near,p,xg);
        

        % 4. Recompute mass
        [mass_mtx]=mass_mat(near,p,dens,volume);

        % 5. Constitutive & Stiffness_mat
        [stiff_mtx,int_Fs,Eps_p,Sy,def_G_p,ss,es,es_p,FAIL]=Constitutive...
        (2,def_G,def_G_p,Eps_p,Sy,near,volume,B,jacobians,...
        ss,es,es_p,FAIL);

        [R]=reaction(int_Fs,x_0,DOF);
        %  6. Assemble time integration matrix

        [matrix]=G_matrix(mass_mtx,stiff_mtx,damp_mtx);
        matrix_1=matrix;

        for i=1:nodes*df
            if (boundary(i)~=0)
                for j=1:nodes*df
                    matrix(i,j)=0;
                    matrix(j,i)=0;
                end
                matrix(i,i)=1;
            end
        end
        InvK=matrix\I;


        % 7. Storage
        if rem(ste,SAVE_I)==0;

            if SAVE_I==1
                ste_p=ste;
            else
                ste_p=ste/SAVE_I+1;
            end
            fprintf('ste_p %i \n',ste_p);

            React(ste_p,1)=R;
    
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

    end
    


end

