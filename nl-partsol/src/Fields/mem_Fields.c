#include "nl-partsol.h"


Fields allocate_Fields(int NumParticles)
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
    Deformation gradient field (Tensor)
  */  
  Phi.F_n = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.F_n.Info,"Deformation gradient at t = n");
  Phi.F_n1 = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.F_n1.Info,"Deformation gradient at t = n + 1");

  /* Initialize the deformation gradient with the identity */
  for(int p = 0 ; p<NumParticles ; p++)
    {
      for(int i = 0 ; i<Ndim ; i++)
	{
	  Phi.F_n.nM[p][i + i*Ndim] = 1.0;	  
	}
      
    }
  
  
  /*!
    Strain_If field (Scalar) 
  */
  Phi.Strain_If = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.Strain_If.Info,"Strain in fracture GP");

  /*!
    Stress field (Tensor)
  */
  Phi.Stress = allocZ__MatrixLib__(NumParticles,Ndim*Ndim);
  strcpy(Phi.Stress.Info,"Stress field GP");

  /*!
    Deformation Energy (Scalar) 
  */
  Phi.W = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.W.Info,"Deformation Energy GP");

  /*!
    Damage parameter (fracture) 
  */
  Phi.chi = allocZ__MatrixLib__(NumParticles,1);
  strcpy(Phi.chi.Info,"Damage parameter GP");

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

  return Phi;
}


/*********************************************************************/

void free_Fields(Fields Phi)
{

  free__MatrixLib__(Phi.rho);
  free__MatrixLib__(Phi.mass);
  free__MatrixLib__(Phi.x_GC);  
  free__MatrixLib__(Phi.x_EC);
  free__MatrixLib__(Phi.dis);
  free__MatrixLib__(Phi.vel);
  free__MatrixLib__(Phi.acc);
  free__MatrixLib__(Phi.Stress);
  free__MatrixLib__(Phi.Strain);
  free__MatrixLib__(Phi.Strain_If);
  free__MatrixLib__(Phi.W);
  free__MatrixLib__(Phi.chi);
  
}

/*********************************************************************/
