#include "nl-partsol.h"

/*********************************************************************/

Fields allocate_U_vars__Fields__(int NumParticles)
{
  int Ndim = NumberDimensions;
  Fields Phi;

  /*!
    Global coordinates 
  */
  Phi.x_GC = allocZ__MatrixLib__(NumParticles,Ndim);
  strcpy(Phi.x_GC.Info,"Global Coordinates");
  
  /*!
    Natural coordinates (Vectorial) 
  */
  Phi.x_EC = allocZ__MatrixLib__(NumParticles,Ndim);
  strcpy(Phi.x_EC.Info,"Element Coordinates GP");
  
  /*!
    Displacement field (Vectorial) 
  */
  Phi.dis = allocZ__MatrixLib__(NumParticles,Ndim);
  strcpy(Phi.dis.Info,"Displacement field GP");
  Phi.D_dis = allocZ__MatrixLib__(NumParticles,Ndim);
  strcpy(Phi.D_dis.Info,"Increment of displacement field GP");

  /*!
    Velocity field (Vectorial) 
  */
  Phi.vel = allocZ__MatrixLib__(NumParticles,Ndim);
  strcpy(Phi.vel.Info,"Velocity field GP");
  
  /*!
    Acceleration field (Vectorial) 
  */
  Phi.acc = allocZ__MatrixLib__(NumParticles,Ndim);
  strcpy(Phi.acc.Info,"Acceleration field GP");
  
  /*!
    Strain field (Tensor)
  */
  Phi.Strain = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.Strain.Info,"Strain field GP");

  /*!
    Deformation gradient field (Tensor) + Initialise it with the indentity
  */  
  Phi.F_n = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.F_n.Info,"Deformation gradient at t = n");
  Phi.F_n1 = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.F_n1.Info,"Deformation gradient at t = n + 1");
  Phi.DF = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.DF.Info,"Increment deformation gradient");

  for(int p = 0 ; p<NumParticles ; p++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      Phi.F_n.nM[p][i + i*Ndim] = 1.0;	  
      Phi.F_n1.nM[p][i + i*Ndim] = 1.0;
      Phi.DF.nM[p][i + i*Ndim] = 1.0;
    }  
  }

  Phi.dt_F_n = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.dt_F_n.Info,"Rate of deformation gradient at t = n");
  Phi.dt_F_n1 = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.dt_F_n1.Info,"Rate of deformation gradient at t = n + 1");
  Phi.dt_DF = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.dt_DF.Info,"Rate of increment deformation gradient");

  /*!
    Inverse plastic deformation gradient field (Tensor) + Initialise it with the indentity
  */
  Phi.F_m1_plastic = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.F_m1_plastic.Info,"Inverse plastic deformation gradient");
 
  for(int p = 0 ; p<NumParticles ; p++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      Phi.F_m1_plastic.nM[p][i + i*Ndim] = 1.0;	  
    }
  }

  /*!
    Jacobian field (Scalar) 
  */
  Phi.J_n = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.J_n.Info,"Jacobian of the particle");  
  Phi.J_n1 = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.J_n1.Info,"Jacobian of the particle");  

  for(int p = 0 ; p<NumParticles ; p++)
  {
    Phi.J_n.nV[p] = 1.0;
    Phi.J_n1.nV[p] = 1.0;  
  }
  
  /*!
   * F-bar variables
   * */
  Phi.Fbar = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.Fbar.Info,"Fbar deformation gradient");

  for(int p = 0 ; p<NumParticles ; p++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      Phi.Fbar.nM[p][i + i*Ndim] = 1.0;   
    }
  }

  Phi.Jbar = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.Jbar.Info,"Jacobian of the particle");  

  for(int p = 0 ; p<NumParticles ; p++)
  {
    Phi.Jbar.nV[p] = 1.0;
  }

  /*!
    Strain_If field (Scalar) 
  */
  Phi.Strain_If = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.Strain_If.Info,"Strain in fracture GP");

  /*!
    Stress field (Tensor)
  */
  Phi.Stress = allocZ__MatrixLib__(NumParticles,Ndim*Ndim + (Ndim == 2? 1 : 0));
  strcpy(Phi.Stress.Info,"Stress field GP");

  /*!
    Deformation Energy (Scalar) 
  */
  Phi.W = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.W.Info,"Deformation Energy GP");

  /*!
    Mass 
  */
  Phi.mass = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.mass.Info,"Mass GP");

  /*!
    Density 
  */
  Phi.rho = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.rho.Info,"Density GP");

  /*!
    Inital volume
  */
  Phi.Vol_0 = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.Vol_0.Info,"Inital volume GP");

  /*!
    Damage parameter (Fracture) 
  */
  Phi.chi = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.chi.Info,"Damage parameter GP");

  /*!
    Cohesion (Plasticity)
  */
  Phi.cohesion = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.cohesion.Info,"Cohesion GP");
 
  /*!
    Equivalent plastic strain (Plasticity)
  */
  Phi.EPS = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.EPS.Info,"EPS GP");

  /*!
   Isotropic hardeing parameter
   */
  Phi.Kappa_hardening = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.Kappa_hardening.Info,"Kappa hardening GP");

  /*! 
  * Back stress for kinematic hardening (plasticity)
  */
  Phi.Back_stress = allocZ__MatrixLib__(NumParticles,3);
  strcpy(Phi.EPS.Info,"Back stress GP");

  return Phi;
}

/*********************************************************************/

void free_U_vars__Fields__(Fields Phi)
{
  free__MatrixLib__(Phi.rho);
  free__MatrixLib__(Phi.mass);
  free__MatrixLib__(Phi.x_GC);  
  free__MatrixLib__(Phi.x_EC);
  free__MatrixLib__(Phi.dis);
  free__MatrixLib__(Phi.D_dis);
  free__MatrixLib__(Phi.vel);
  free__MatrixLib__(Phi.acc);
  free__MatrixLib__(Phi.Stress);
  free__MatrixLib__(Phi.Strain);
  free__MatrixLib__(Phi.Strain_If);
  free__MatrixLib__(Phi.F_m1_plastic);
  free__MatrixLib__(Phi.F_n);
  free__MatrixLib__(Phi.F_n1);
  free__MatrixLib__(Phi.J_n);
  free__MatrixLib__(Phi.J_n1);
  free__MatrixLib__(Phi.dt_F_n);
  free__MatrixLib__(Phi.dt_F_n1);
  free__MatrixLib__(Phi.dt_DF);
  free__MatrixLib__(Phi.Fbar);
  free__MatrixLib__(Phi.Jbar);
  free__MatrixLib__(Phi.DF);
  free__MatrixLib__(Phi.W);
  free__MatrixLib__(Phi.Vol_0);
  free__MatrixLib__(Phi.chi);
  free__MatrixLib__(Phi.cohesion);
  free__MatrixLib__(Phi.EPS);
  free__MatrixLib__(Phi.Back_stress);
}

/*********************************************************************/