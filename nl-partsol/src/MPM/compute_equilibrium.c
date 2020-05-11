#include "nl-partsol.h"

/*************************************************************/

Matrix compute_equilibrium_U(Matrix V_I, GaussPoint MPM_Mesh,
			     Mesh FEM_Mesh, double TimeStep)
/*

*/
{
  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;
  
  update_LocalState(V_I, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
       
  Matrix F_I = MatAllocZ(Nnodes,Ndim);    

  F_I = compute_InternalForces(F_I, MPM_Mesh, FEM_Mesh);    

  F_I = compute_BodyForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);

  F_I = compute_ContacForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);

  return F_I;
}

/*************************************************************/
