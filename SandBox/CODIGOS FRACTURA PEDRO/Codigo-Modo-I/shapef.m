function [p,dp,lambda,FAIL]=shapef...
    (x_a,beta_,x_sample,near_,TolLag,Nelder,wrap,p,lambda)
    % Calculation of the local maximum-entropy shape functions with Newton's
    % method as described in section 4.2 of [1]
    %
    % INPUT
    % =====
    % ndim:     spacial dimension
    % x_a:      nodal set (nnode, ndim)
    % beta_:    value of the thermalization parameter at each sample point
    % x_sample: sample point where the shape functions are evaluated
    % near_:    near_{i} contains the list of nodes with non-zero shape funcion
    %           at the i-th sample point
    % TolLag:   tolerance for Newton's iterations
    %
    % OUTPUT
    % ======
    % p:        p{i} contains the values of the shape functions corresponding
    %           to the nodes in near_{i} at the i-th sample point
    % dp:       dp{i} contains the values of the ndim spacial derivatives of 
    %           the shape functions corresponding to the nodes in near_{i} at 
    %           the i-th sample point
    %
    % Reference:
    % [1] Marino Arroyo and Michael Ortiz, "Local maximum-entropy approximation
    %     schemes: a seamless bridge between finite elements and meshfree methods", 
    %     International Journal for Numerical Methods in Engineering, 65:2167-2202 (2006). 


    [n_sample,ndim]=size(x_sample);
    FAIL=0;

    for i=1:n_sample        
        if wrap(i)==1
              if i==24
                  i;
              end
            
            clear p(i) p_a
        
              if (ndim==1)
                x=x_sample(i);
              else
                x=x_sample(i,:);
              end
              beta=beta_(i);
              near=near_{i};

              lama=lambda(i*ndim-1:i*ndim);
              lam=lama';

              if Nelder
                  [p_a,J,cnt,lam,tol]=Nelder_LME_mod(i,x,x_a,beta,near,TolLag,lam,ndim);
              else
                  [p_a,J,~]=NR_LME(x,x_a,beta,near,TolLag,lam,ndim);
              end
              
              if cnt>=250
                  fprintf('Convergence problem in element %i with tol %i'...
                      ,i,tol);
              end
              
              lambda(i*ndim-1:i*ndim)=lam';
              p(i)={p_a}; 

              if isnan(p_a)
                  fprintf('Not a number in element %i \n',i);
                  FAIL=1;
                  break;
              elseif isnan(J)
                  fprintf('Not a number in element %i \n',i);
                  FAIL=1;
                  break;
              else
                    %Spacial Gradients
                  dp(i)={ zeros(length(near),ndim) };
                  for ia=1:length(near)
                    if ndim==1
                      dp{i}(ia,1)= -(p_a(ia) * (x-x_a(near(ia))') )/J(1,1);      
                    else
                      dp{i}(ia,:)= -J\(p_a(ia) * (x(:)-x_a(near(ia),:)') );
                    end
                  end
              end
        else
        	dp(i)={0};
        end
    end   %nsample points
    
end

