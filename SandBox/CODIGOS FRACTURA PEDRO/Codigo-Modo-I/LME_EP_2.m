

function [B,near,p,gamma_,lam_LME,REMAP,wrap,EP,FAIL]=LME_EP_2(jacobians,...
    volume,x_a,xg,B,near,p,gamma_,lam_LME,wrap,EP)


    global elements sp h_ini Material Mat_nds FB_BAR LME_param ...
    patch_con patch_el
    
    h=h_ini.*sqrt(jacobians);

    % This sets the numerical threshold for the support of the shape functions
    target_zero= LME_param(2); 
    % This is the tolerance for the Newton iteration to compute the shape
    % functions
    TolLag= LME_param(3);
    % Nelder Mead method to obtain minimum lambda in shape functions
    Nelder=LME_param(4);
    %Tolerance to make new re-mapping or not and proportion to do it
    tol_search = LME_param(5);
    prop = LME_param(6);
    %Minimum and initial gamma 
    gamma_top=LME_param(7);
    gamma_lme=LME_param(1);

    %% WRAP or not?
    beta_=zeros(elements,1);
    REMAP=0;
    dis(3)=0;
    for i=1:elements
        if wrap(i)==1
            REMAP=1;
        else
            Ep=EP((i-1)*3+1:(i-1)*3+3,:);
            for j=1:3
                dis(j)=abs(sqrt(Ep(j,1))-sqrt(Ep(j,2)));
            end
            Dis=max(dis);
            if Dis>tol_search
                m_dis=max(Ep(:,1));
                for j=1:3
                    EP((i-1)*3+j,2)=EP((i-1)*3+j,1);
                end
                wrap(i)=1;
                %plot_nb(i, near, x_a, xg, elem_0)
                REMAP=1;
                if m_dis>0
                    gamma_(i)=max(gamma_lme-sqrt(m_dis)*prop,gamma_top);
                end
            end 
        end
        beta_(i)=gamma_(i)/h(i)^2;
    end
    
    %% MAKE NEAR
    range1=zeros(elements,1);
    range_lme=zeros(elements,1);
    for i=1:elements
        range1(i)=h(i)*sqrt(-1/gamma_(i)*log(target_zero));
        range_lme(i)=max(range1(i)*1, h(i));
    end
    [near]=make_near(x_a,xg,range_lme,wrap,near,Material,Mat_nds);


    %% MAKE SHAPE FUNCTION
    if sp==2
        [p,dp,lam_LME,FAIL]=shapef...
            (x_a,beta_,xg,near,TolLag,Nelder,wrap,p,lam_LME);
    end

    %% MAKE Bbar - Patch
    if FAIL==0
        if FB_BAR==1
        [B,near,p]=bbar_matrix...
            (x_a,patch_con,patch_el,elements,near,dp,p,volume,B,wrap,FB_BAR);
        else
        [B]=b_m(elements,near,dp,sp,wrap,B);
        end
    end

end