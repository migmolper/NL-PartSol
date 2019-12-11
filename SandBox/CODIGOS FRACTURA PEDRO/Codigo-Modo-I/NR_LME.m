function [p_a,J,niter]=NR_LME(x,x_a,beta,near,TolLag,lam,ndim)


  R=10*ones(1,ndim);
  niter=0;
  I=eye(ndim);
  
  %Newton iteration
  iflag=0;
  dlam=10*ones(1,ndim);
  while (norm(R)>TolLag)
    [gam,R,J,p_a]=Gamma_(ndim,x_a,x,beta,lam,near);
    if (abs(rcond(J)))<1e-8
      iflag=1;
      disp('Newton Failed, near to singular matrix')
    end
    inv=-J\I;
    dlam=-J\R';
    lam=lam+dlam';
    niter=niter+1;
    if (niter>100) 
      i, niter
      iflag=1;
      disp('Newton Failed 2, no convergence in 100 iterations')
    end
  end
  
end
  

        


