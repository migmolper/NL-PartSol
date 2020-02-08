
#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* Dourant condition */
double DeltaT_CFL(GaussPoint MPM_Mesh, double h);

/*!
 * Displacement formulation with a Forward-Euler 
 * integration scheme.
 */
void u_ForwardEuler(Mesh FEM_Mesh, GaussPoint GP_Mesh);

/*! 
 * Displacement formulation with a Generalized-alpha
 * integration scheme.
 */
void u_GeneralizedAlpha(Mesh FEM_Mesh, GaussPoint MPM_Mesh);
