#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

/*******************************************************/
void update_NodalMomentum(Mesh FEM_Mesh, Matrix Phi_I, Matrix F_I)
/*!
 * \brief Brief description of UpdateGridNodalMomentum.
 *        Compute the nodal contribution of each GP to the total forces.
 *
 *  The parameters for this functions are  :
 *  @param MPM_Mesh : Mesh with the material points.
 *  @param Phi_I : Nodal value of the momentum.
 *  @param F_I : Nodal value of the total forces.
 *
 */
{
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  
  /* Update the grid nodal momentum */
  for(int i = 0 ; i<Nnodes ; i++){
    for(int j = 0 ; j<Ndim ; j++){
      if(FEM_Mesh.ActiveNode[i] > 0){
	Phi_I.nM[i][j] += DeltaTimeStep*F_I.nM[i][j];
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
  int Ndim = NumberDimensions;
  double SizeTable = Ndim*sizeof(double *);
  double Mass_I;

  /* Time integration parameters */
  double alpha = Params.GA_alpha;
  double gamma = Params.GA_gamma;

  /* Asign forces to an auxiliar variable */
  Matrix F = Nodal_Forces;
  
  /* Nodal values the fields */
  Matrix Nodal_Acceleration_t0 =
    MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));
  Matrix Nodal_Acceleration_t1 =
    MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));
  Matrix Nodal_Velocity =
    MatAssign(Ndim,N_Nodes,NAN,NULL,(double**)malloc(SizeTable));

  /* 1º Asign Kinetics values to matricial tables */
  for(int i = 0 ; i<Ndim ; i++){
    Nodal_Acceleration_t0.nM[i] = Nodal_Kinetics.nM[1+i];
    Nodal_Acceleration_t1.nM[i] = Nodal_Kinetics.nM[1+Ndim+i];
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+2*Ndim+i];
  }
  /* 2º Update the grid nodal variales */
  for(int i = 0 ; i<N_Nodes ; i++){
    Mass_I = Nodal_Kinetics.nM[0][i];
    if(Mass_I > 0){
      for(int j = 0 ; j<Ndim ; j++){
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
