
#ifndef TypeDefinitions
#define TypeDefinitions
#endif


/* Dourant condition */
double DeltaT_CFL(GaussPoint, double);

/*!
 * Displacement formulation with a Forward-Euler 
 * integration scheme.
 */
void u_ForwardEuler(Mesh, GaussPoint);

/*! 
 * Displacement formulation with a Generalized-alpha
 * integration scheme.
 */
void U_GA(Mesh, GaussPoint);


/*!
 * Explicit predictor corrector 
 */ 
void U_PCE(Mesh, GaussPoint);
