/*! \file Formulations.h
    \brief File with the prototype of time integration scheme
*/

#ifndef _FORMULATIONS_H_
#define _FORMULATIONS_H_


/*!
  \fn double DeltaT_CFL(GaussPoint MPM_Mesh, double DeltaX)

  \brief  Get the time step using \cite Anderson_1987

  Inputs :
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param double DeltaX : Minimum mesh size
*/
double DeltaT_CFL(GaussPoint, double);

/*!
  Displacement formulation with a Forward-Euler 
  integration scheme.
 */
void U_FE(Mesh, GaussPoint);

/*! 
  Displacement formulation with a Generalized-alpha
  integration scheme.
 */
void U_GA(Mesh, GaussPoint);


/*!
  Explicit predictor corrector 
 */ 
void U_PCE(Mesh, GaussPoint);

#endif
