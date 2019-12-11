
function [def_G,jacobians,volume,dens,EP]=update_F(d,near,B,def_G,volume,EP)
    
    global elements sp Area Area_p FB_BAR patch_con MAT MAT_el thickness
  
    [aa,bb]=size(patch_con);
    area_p_n=zeros(aa,1);
        
    m_I=[1 0 0 1 1];
    f_v=zeros(sp*sp+1,1);
    dens=zeros(elements,1);
    jacobians=zeros(elements,1);
        
    % Update **********************
    for j=1:aa
        for k=1:bb
            e=patch_con(j,k);
            nn=length(near{e});
            nd = near{e};

            b=B{e};
            b2=zeros;
            for i=2:2:nn*2
                b2(1,i-1)=b(1,i-1);
                b2(2,i-1)=b(2,i);
                b2(3,i)=b(3,i);
                b2(4,i)=b(3,i-1);
            end 

            u2=zeros(nn*sp,1);     %Incremental displacements

            for i=1:nn
                u2(i*sp-1,1)=d(nd(i)*sp-1,1)-d(nd(i)*sp-1,2);
                u2(i*sp,1)=d(nd(i)*sp,1)-d(nd(i)*sp,2);
            end

            dF_=b2*u2;                % Incremental F
            dF_(sp*sp+1)=0;
            dF_=dF_+m_I';

            %Vector F to Matrix F
            for i=1:5
                f_v(i,1)=def_G((e-1)*5 + i,2);
            end           
            [F]=v2m(f_v,sp);
            [dF]=v2m(dF_,sp);
            F=dF*F;
            jacobians(e)=det(F);

            if FB_BAR==2
                %%%%%%   F_BAR   %%%%%%
                F_store(k)={F};
                % New vol, dens, jaco
                area_p_n(j)=area_p_n(j)+Area(e)*jacobians(e)*thickness;
            else
                % New vol, dens, jaco
                volume(e)=Area(e)*jacobians(e)*thickness;
                dens(e)=MAT(3,MAT_el(e))/jacobians(e);

                %Storage of vector F
                [f]=m2v(F,sp);

                for i=1:5
                    def_G(e*5+1-i,1)=f(6-i);  
                end
                
                C=F'*F;
                %Principal strecht
                [EP_]=Principal(C);
                for i=1:3
                    EP((e-1)*3+i,1)=EP_(i);
                end
            end
        end
        
        %%%%%%   F_BAR   %%%%%%
        if FB_BAR==2
            J_bar=area_p_n(j)/Area_p(j);
            for k=1:bb
                F=F_store{k};
                e=patch_con(j,k);

                F_=(J_bar/jacobians(e))^(1/2)*F(1:2,1:2);
                F_(3,3)=1;

                % New vol, dens, jaco
                jacobians(e)=det(F_);
                volume(e)=Area(e)*jacobians(e)*thickness;
                dens(e)=MAT(3,MAT_el(e))/jacobians(e);
                
                %Storage of vector F
                [f]=m2v(F_,sp);

                for i=1:5
                    def_G(e*5+1-i,1)=f(6-i);  
                end
                
                C=F_'*F_;
                %Principal strecht
                [EP_]=Principal(C);
                for i=1:3
                    EP((e-1)*3+i,1)=EP_(i);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end



function [EP]=Principal(C)

    [eps]=eig(C);

    [EP(1),i]=max(eps);
    eps(i)=-1e32;
    [EP(2),i]=max(eps);
    eps(i)=-1e32;
    EP(3)=max(eps);


end




