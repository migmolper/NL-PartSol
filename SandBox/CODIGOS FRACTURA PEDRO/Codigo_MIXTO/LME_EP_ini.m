
function [B,near,p,gamma_,lam_LME,EP,FAIL]=...
    LME_EP_ini(jacobians,volume,x_a,xg,NO)
    
    global elements sp h_ini elem Material Mat_nds FB_BAR LME_param ...
        patch_con patch_el
    
    h=h_ini.*sqrt(jacobians);
    
    B={0};
    near={0};
    [near]=first_near(elem,near);
    p={0};

    % As described in [1], gamma measures the degree of locality, the larger,
    % the closer to Delaunay finite elements

    gamma_lme=1.8;
    gamma_top=0.9;
    gamma_=gamma_lme*ones(elements,1);
        
    % This sets the numerical threshold for the support of the shape functions
    target_zero=1.e-6; 
    % This is the tolerance for the Newton iteration to compute the shape
    % functions
    TolLag=max(2*eps,1.e-18);
    % Nelder Mead method to obtain minimum lambda in shape functions
    Nelder=1;
    
    %Tolerance to make new re-mapping or not and proportion to do it
    tol_search = 0.5;       % Optimum [0.4-0.7]
    prop = 0.1;             % Proportion of gamma decreasing
    
    LME_param(1)=gamma_lme;
    LME_param(2)=target_zero;
    LME_param(3)=TolLag;
    LME_param(4)=Nelder;
    LME_param(5)=tol_search;
    LME_param(6)=prop;
    LME_param(7)=gamma_top;
    
    %% WRAP or not?
    beta_=zeros(elements,1);
    wrap=zeros(elements,1);
    EP=zeros(3*elements,2);
    for i=1:elements
        wrap(i)=1;
        for j=1:3
            EP((i-1)*3+j,1)=1;
            EP((i-1)*3+j,2)=1;
        end  
        beta_(i)=gamma_(i)/h(i)^2;
    end
    
    [lam_LME]=first_lambda(elem,xg,x_a,beta_);   %First lambda
    
    %% MAKE NEAR
    range1=zeros(elements,1);
    range_lme=zeros(elements,1);
    for i=1:elements
        range1(i)=h(i)*sqrt(-1/gamma_(i)*log(target_zero));
        range_lme(i)=max(range1(i)*1, 2*h(i));
    end
    [near]=make_near(x_a,xg,range_lme,wrap,near,Material,Mat_nds,NO);


    %% MAKE SHAPE FUNCTION
    if sp==2
        [p,dp,lam_LME,FAIL]=shapef...
            (x_a,beta_,xg,near,TolLag,Nelder,wrap,p,lam_LME);
    end

    %% MAKE Bbar - Patch
    if FB_BAR==1
        [B,near,p]=bbar_matrix...
            (x_a,patch_con,patch_el,elements,near,dp,p,volume,B,wrap,FB_BAR);
    else
        [B]=b_m(elements,near,dp,sp,wrap,B);
    end

end

function [near]=first_near(elem,near)
    [n_sample,~]=size(elem);
    for i=1:n_sample
        near(i)={elem(i,:)};
    end
end

function [lambda]=first_lambda(elem,x_sample,x_a,beta_)

[n_sample,bg]=size(elem);
[~,sp]=size(x_a);

lambda=zeros(sp*n_sample,1);

    for i=1:n_sample
        
          x=x_sample(i,:);
          beta=beta_(i);

        % Initialize lambda
          y=zeros(bg,sp);
          for j=1:bg
            for k=1:sp
                y(j,k)=x(k)-x_a(elem(i,j),k);
            end
          end
          A(1,1)=y(2,1)-y(1,1);
          A(1,2)=y(2,2)-y(1,2);
          A(2,1)=y(3,1)-y(1,1);
          A(2,2)=y(3,2)-y(1,2);
          B(1,1)=-beta*(y(1,1)^2+y(1,2)^2-(y(2,1)^2+y(2,2)^2));
          B(2,1)=-beta*(y(1,1)^2+y(1,2)^2-(y(3,1)^2+y(3,2)^2));
          lam=A\B;
          
          lambda(i*sp-1:i*sp)=lam;
    end
end