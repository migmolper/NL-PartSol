
#ifndef TypeDefinitions
#define TypeDefinitions
#endif

/*
  Displacement formulation with a Forward-Euler 
  integration scheme.
*/
void u_ForwardEuler(Mesh FEM_Mesh, GaussPoint GP_Mesh);

/* 
   Stress-velocity formulation with a two-step Taylor-Galerkin
   integration sheme.
*/
/* void SigmaV_TSTG(Mesh FEM_Mesh, GaussPoint GP_Mesh); */
