function [stiff_mtx,int_Fs,status,D_energy,ss,es,FAIL]=Constitutive_frac...
    (Kt,def_G,near,volume,B,jacobians,h,ss,es,status,epsilon,D_energy,FAIL)

    global nodes sp elements MODEL MAT MAT_el
    
    f_v       = zeros(sp*sp+1,1);
    int_Fs    = zeros(nodes*sp,1);
    stiff_mtx = zeros(sp*nodes,sp*nodes);
    
    Wlist=zeros(elements,1);

    if FAIL==0
        for e=1:elements
            
            %Vector F, Fp to Matrixes
            for i=1:5
                f_v(i,1)=def_G((e-1)*5 + i,1);
            end           
            [F]=v2m(f_v,sp);
            
            if status(e,2)==0
                w=1;
                if MODEL==0
                    [A,T,W,Ee]=Saint_Venant(Kt,e,F,w);
                elseif MODEL==1
                    [A,T,W,Ee]=Neo_Hookean(Kt,e,F,jacobians,w);
                end
            else
                T=zeros(3);
                W=0;
            end
            
            Cauchy(e)={T};
            Wlist(e)=W;
        end

        for e=1:elements
            
            Ceps=MAT(8,MAT_el(e));
 
            %Derivatives           
            nn=length(near{e});
            nd = near{e};

            B_=B{e};
            sh=zeros(2,nn);
            for i=1:nn
                sh(1,i)=B_(1,i*2-1);
                sh(2,i)=B_(2,i*2);
            end
            
            
            T=Cauchy{e};
            
            AA=isreal(T);
            if  AA(1)==0
                FAIL=1;
                break;
            else    
                if status(e,2)==0 
                    eps_list=epsilon{e};
                    
                    if eps_list(1)>0
                        [S_prin]=Principal(T);
                        crit=S_prin(1)+S_prin(2)+S_prin(3);
                        
                        if crit>0                      
                            % ------ G calculation -----
                            G=Wlist(e)*volume(e);
                            M=volume(e);
                            for i=1:length(eps_list)
                                G=G+Wlist(eps_list(i))*volume(eps_list(i));
                                M=M+volume(eps_list(i));
                            end
                            G=G*Ceps*h(e)/M;

                            % ------- G > Gc ?? --------
                            if G > MAT(7,MAT_el(e))
                                status(e,1)=1;
                                T=zeros(3);
                                D_energy = D_energy + Wlist(e)*volume(e);
                            end
                        end
                    end
                else
                    T=zeros(3);
                end  
                
                % ----------------------------
                % Internal forces
                % ----------------------------
                Tt=T(1:2,1:2);
                int_forces_1=Tt*sh*volume(e);

                for i=1:nn
                   nod=nd(i);
                   for j=1:sp
                        int_Fs(nod*sp+1-j,1)=int_Fs(nod*sp+1-j,1)+int_forces_1(3-j,i);
                   end
                end
                
                % ----------------------------
                % Stiffness matrix
                % ----------------------------
                if Kt==1 || Kt==2
                    [stiff_mtx]=stiff_mat(B,near,volume,e,stiff_mtx,T,A);
                end
                
                % ----------------------------
                % Update
                % ----------------------------
                if Kt==3 || Kt==2
                    [sig]=E2e(T);
                    [ee]=E2e(Ee);
                    for i=1:4
                        ss((e-1)*4+i,1)=sig(i,1);
                        es((e-1)*4+i,1)=ee(i,1);
                    end
                   
                end 
            end
            clear sh
        end
        status(:,2)=status(:,1);
    end   
end

function [e]=E2e(E)
   
    e=zeros(4,1);

    %Build vector
    e(1)=E(1,1);
    e(2)=E(2,2);
    e(3)=E(3,3);
    e(4)=E(1,2);%+E(2,1);
end

function [EP]=Principal(C)

    [eps]=eig(C);

    [EP(1),i]=max(eps);
    eps(i)=-1e32;
    [EP(2),i]=max(eps);
    eps(i)=-1e32;
    EP(3)=max(eps);
end