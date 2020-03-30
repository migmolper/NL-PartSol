#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

#define MAXVAL(A,B) ((A)>(B) ? (A) : (B))
#define MINVAL(A,B) ((A)<(B) ? (A) : (B))

/*********************************************************************/

double DeltaT_CFL(GaussPoint MPM_Mesh, double h)
/*
  Get the time step using :
  C.E. Anderson Jr. An overview of the theory of hydrocodes (1987) 
*/
{

  double DeltaT;
  double C[3] = {0,0,0};
  int N_dim = 3;

  /* Get the maximum wave speed */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    for(int j = 0 ; j<N_dim ; j++){
      C[j] = MAXVAL(C[j],CEL+MPM_Mesh.Phi.vel.nM[i][j]);
    }
  }

  /* Get the minimum value of the time step */
  DeltaT = h/C[0];
  for(int j = 1 ; j<N_dim ; j++){
    DeltaT = MINVAL(DeltaT,h/C[j]);
  }

  /* Step multiplier ranging 0-1 */
  DeltaT *= CFL;

  /* Return new time step */
  return DeltaT;
}

/*********************************************************************/
