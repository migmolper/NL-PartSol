function [stiff_mtx,int_Fs,Gamma,Sy,def_G_p,ss,es,es_p,FAIL]=Constitutive...
    (Kt,def_G,def_G_p,Gamma,Sy,near,volume,B,jacobians,ss,es,es_p,FAIL)

    global nodes sp elements MODEL
    
    f_v       = zeros(sp*sp+1,1);
    f_p       = zeros(sp*sp+1,1);
    int_Fs    = zeros(nodes*sp,1);
    stiff_mtx = zeros(sp*nodes,sp*nodes);
    
    Ep=zeros(3);
    Fp=zeros(3);

    if FAIL==0
        for e=1:elements
            
            %Vector F, Fp to Matrixes
            for i=1:5
                f_v(i,1)=def_G((e-1)*5 + i,1);
                f_p(i,1)=def_G_p((e-1)*5 + i,2);
            end           
            [F]=v2m(f_v,sp);
            [Fp]=v2m(f_p,sp);
            
            nn=length(near{e});
            nd = near{e};

            %Derivatives
            B_=B{e};
            sh=zeros(2,nn);
            for i=1:nn
                sh(1,i)=B_(1,i*2-1);
                sh(2,i)=B_(2,i*2);
            end

            if MODEL==0
                [A,T,Ee]=Saint_Venant(Kt,e,F);
            elseif MODEL==1
                [A,T,Ee]=Neo_Hookean(Kt,e,F,jacobians);
            elseif MODEL>1 
                [A,T,Gamma,Sy,Ee,Ep,Fp]=Drucker_prager(Kt,e,Gamma,Sy,F,Fp);
            end

            AA=isreal(T);
            if  AA(1)==0
                FAIL=1;
                break;
            else
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
                    [fp]=m2v(Fp,sp);
                    [ee]=E2e(Ee);
                    [ep]=E2e(Ep);
                    for i=1:5
                        def_G_p((e-1)*5+i,1)=fp(i,1);
                    end 
                    for i=1:4
                        ss((e-1)*4+i,1)=sig(i,1);
                        es((e-1)*4+i,1)=ee(i,1);
                        es_p((e-1)*4+i,1)=ep(i,1);
                    end
                   
                end 
            end
            clear sh
        end
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