#include "nl-partsol.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

double DeltaT_CFL(GaussPoint MPM_Mesh, double h)
{

  double DeltaT;
  double CEL_MAT = 0;
  double CEL_MAX = 0;
  double C[3] = {0,0,0};
  int Ndim = NumberDimensions;
  int Nmat = MPM_Mesh.NumberMaterials;
  bool DynamicTimeStep = true;

  /*
    Get the maximum material celerity
   */
  for(int i = 0 ; i<Nmat ; i++)
    {
      CEL_MAT = MAXVAL(MPM_Mesh.Mat[i].Cel,CEL_MAT);
    }

  /*
    Consider the velocity of the particles for the courant. In
    some cases, for instance Fr>1 is important.
   */
  if(DynamicTimeStep)
    {
    /*
      Get the maximum wave speed in any direction
    */
    for(int i = 0 ; i<MPM_Mesh.NumGP ; i++)
      {
	for(int j = 0 ; j<Ndim ; j++)
	  {
	    C[j] = MAXVAL(C[j],CEL_MAT+fabs(MPM_Mesh.Phi.vel.nM[i][j]));
	  }
      }

    /*
      Get the absolute maximum value of the celerity
    */
    for(int j = 0 ; j<Ndim ; j++)
      {
	CEL_MAX = MAXVAL(CEL_MAX,C[j]);
      }
  }
  else
    {
      CEL_MAX = CEL_MAT;
    }

  
  /*
    Get the minimum value of the time step 
  */
  DeltaT = CFL*h/CEL_MAX;


  /*
    Return new time step 
  */
  return DeltaT;
}

/*********************************************************************/
