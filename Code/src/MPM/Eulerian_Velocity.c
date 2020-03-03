#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*******************************************************/

void UpdateGridNodalMomentum(Mesh FEM_Mesh,
			     Matrix Nodal_MOMENTUM,
			     Matrix Nodal_TOT_FORCES)
/*!
 * \brief Brief description of UpdateGridNodalMomentum.
 *        Compute the nodal contribution of each GP to the total forces.
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param Nodal_MOMENTUM : Nodal value of the momentum.
 *  @param Nodal_TOT_FORCES : Nodal value of the total forces.
 *
 */
{
  /* Update the grid nodal momentum */
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<3 ; j++){
      if(FEM_Mesh.ActiveNode[i] > 0){
	Nodal_MOMENTUM.nM[j][i] +=
	  DeltaTimeStep*Nodal_TOT_FORCES.nM[j][i];
      }
    }
  }  
}

/*******************************************************/

void GA_UpdateNodalKinetics(Mesh FEM_Mesh,
			    Matrix Nodal_Kinetics,
			    Matrix Nodal_Forces,
			    Time_Int_Params Params)
/*!
 * \brief Brief description of Update_Nodal_Acceleration_Velocity.
 *
 *  The parameters for this functions are  :
 *  @param FEM_Mesh
 *  @param Nodal_Kinetics = {m, a0, a1, v}
 *  @param Nodal_Forces
 *  @param Params
 */
{
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = 3;
  double SizeTable = N_dim*sizeof(double *);
  double Mass_I;

  /* Time integration parameters */
  double alpha = Params.GA_alpha;
  double gamma = Params.GA_gamma;

  /* Asign forces to an auxiliar variable */
  Matrix F = Nodal_Forces;
  
  /* Nodal values the fields */
  Matrix Nodal_Acceleration_t0 =
    MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));
  Matrix Nodal_Acceleration_t1 =
    MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));
  Matrix Nodal_Velocity =
    MatAssign(N_dim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));

  /* 1º Asign Kinetics values to matricial tables */
  for(int i = 0 ; i<N_dim ; i++){
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1+i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1+N_dim+i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+2*N_dim+i];
  }
  /* 2º Update the grid nodal variales */
  for(int i = 0 ; i<N_Nodes ; i++){
    Mass_I = Nodal_Kinetics.nM[0][i];
    if(Mass_I > 0){
      for(int j = 0 ; j<N_dim ; j++){
	/* Get the nodal acceleration t1 */
	Nodal_Acceleration_t1.nM[j][i] =
	  (F.nM[j][i]/Mass_I - alpha*Nodal_Acceleration_t0.nM[j][i])/(1-alpha);
	/* Update nodal velocity */
	Nodal_Velocity.nM[j][i] +=
	  ((1-gamma)*Nodal_Acceleration_t0.nM[j][i] +
	   gamma*Nodal_Acceleration_t1.nM[j][i])*DeltaTimeStep;
      }
    }
  }

  /* 3º Free tables */
  free(Nodal_Acceleration_t0.nM);
  free(Nodal_Acceleration_t1.nM);
  free(Nodal_Velocity.nM);
  
}
