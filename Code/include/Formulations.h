#ifndef _FORMULATIONS_H_
#define _FORMULATIONS_H_

/* Dourant condition */
double DeltaT_CFL(GaussPoint, double);

/*!
 * Displacement formulation with a Forward-Euler 
 * integration scheme.
 */
void U_FE(Mesh, GaussPoint);

/*! 
 * Displacement formulation with a Generalized-alpha
 * integration scheme.
 */
void U_GA(Mesh, GaussPoint);


/*!
 * Explicit predictor corrector 
 */ 
void U_PCE(Mesh, GaussPoint);



#endif
