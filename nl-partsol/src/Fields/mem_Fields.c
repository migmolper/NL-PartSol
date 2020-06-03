#include "nl-partsol.h"


Fields allocate_Fields(int NumParticles)
{
  int Ndim = NumberDimensions;
  Fields Phi;

  /*!
    Global coordinates 
  */
  Phi.x_GC = MatAllocZ(NumParticles,Ndim);
  strcpy(Phi.x_GC.Info,"Global Coordinates");
  
  /*!
    Natural coordinates (Vectorial) 
  */
  Phi.x_EC = MatAllocZ(NumParticles,Ndim);
  strcpy(Phi.x_EC.Info,"Element Coordinates GP");
  
  /*!
    Displacement field (Vectorial) 
  */
  Phi.dis = MatAllocZ(NumParticles,Ndim);
  strcpy(Phi.dis.Info,"Displacement field GP");
  
  /*!
    Velocity field (Vectorial) 
  */
  Phi.vel = MatAllocZ(NumParticles,Ndim);
  strcpy(Phi.vel.Info,"Velocity field GP");
  
  /*!
    Acceleration field (Vectorial) 
  */
  Phi.acc = MatAllocZ(NumParticles,Ndim);
  strcpy(Phi.acc.Info,"Acceleration field GP");
  
  /*!
    Strain field (Tensor)
  */
  Phi.Strain = MatAllocZ(NumParticles,Ndim*Ndim);
  strcpy(Phi.Strain.Info,"Strain field GP");

  /*!
    Strain_If field (Scalar) 
  */
  Phi.Strain_If = MatAllocZ(NumParticles,1);
  strcpy(Phi.Strain_If.Info,"Strain in fracture GP");

  /*!
    Stress field (Tensor)
  */
  Phi.Stress = MatAllocZ(NumParticles,Ndim*Ndim);
  strcpy(Phi.Stress.Info,"Stress field GP");

  /*!
    Deformation Energy (Scalar) 
  */
  Phi.W = MatAllocZ(NumParticles,1);
  strcpy(Phi.W.Info,"Deformation Energy GP");

  /*!
    Damage parameter (fracture) 
  */
  Phi.ji = MatAllocZ(NumParticles,1);
  strcpy(Phi.ji.Info,"Damage parameter GP");

  /*!
    Mass 
  */
  Phi.mass = MatAllocZ(NumParticles,1);
  strcpy(Phi.mass.Info,"Mass GP");

  /*!
    Density 
  */
  Phi.rho = MatAllocZ(NumParticles,1);
  strcpy(Phi.rho.Info,"Density GP"); 

  return Phi;
}


/*********************************************************************/

void free_Fields(Fields Phi)
{

  FreeMat(Phi.rho);
  FreeMat(Phi.mass);
  FreeMat(Phi.x_GC);  
  FreeMat(Phi.x_EC);
  FreeMat(Phi.dis);
  FreeMat(Phi.vel);
  FreeMat(Phi.acc);
  FreeMat(Phi.Stress);
  FreeMat(Phi.Strain);
  FreeMat(Phi.Strain_If);
  FreeMat(Phi.W);
  FreeMat(Phi.ji);
  
}

/*********************************************************************/
