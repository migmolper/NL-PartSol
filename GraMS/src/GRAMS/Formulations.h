
#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* Dourant condition */
double DeltaT_CFL(GaussPoint MPM_Mesh, double h);

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
