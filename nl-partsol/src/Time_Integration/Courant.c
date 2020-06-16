#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

double DeltaT_CFL(GaussPoint MPM_Mesh, double h)
{

  double DeltaT;
  double CEL_MAX = 0;
  double C[3] = {0,0,0};
  int Ndim = NumberDimensions;
  int Nmat = MPM_Mesh.NumberMaterials;

  for(int i = 0 ; i<Nmat ; i++)
    {
      CEL_MAX = MAXVAL(MPM_Mesh.Mat[i].Cel,CEL_MAX);
    }
  
  /*!
    Get the maximum wave speed
  */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
  	{
  	  /* C[j] = MAXVAL(C[j],CEL_MAX+MPM_Mesh.Phi.vel.nM[i][j]); */
	  C[j] = MAXVAL(C[j],CEL_MAX);
  	}
    }

  /*! 
    Get the minimum value of the time step 
  */
  DeltaT = h/C[0];
  for(int j = 1 ; j<Ndim ; j++)
    {
      DeltaT = MINVAL(DeltaT,h/C[j]);
    }

  /*!
    Step multiplier ranging 0-1 
  */
  DeltaT *= CFL;

  /*! 
    Return new time step 
  */
  return DeltaT;
}

/*********************************************************************/
