/*! \file Formulations.h
    \brief File with the prototype of time integration scheme
*/

#ifndef _FORMULATIONS_H_
#define _FORMULATIONS_H_


/*!
  \fn double DeltaT_CFL(GaussPoint MPM_Mesh, double DeltaX)

  \brief  Get the time step using \cite Anderson_1987

  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param double DeltaX : Minimum mesh size
*/
double DeltaT_CFL(GaussPoint, double);

/*!
  \fn void U_FE(Mesh FEM_Mesh, GaussPoint MPM_Mesh)

  \brief Displacement formulation of the MPM with a Forward Euler as 
  time integrator scheme. The algorithm was taken from \cite Zhang_book_2016

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
 */
void U_FE(Mesh, GaussPoint);


/*!
  \fn void U_GA(Mesh FEM_Mesh, GaussPoint MPM_Mesh)

  \brief The generalized-alpha algorithm here implemented is analogous
  to the one described in \cite Tran2019e

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
 */
void U_GA(Mesh, GaussPoint);

/*!
  \fn void U_PCE(Mesh FEM_Mesh, GaussPoint MPM_Mesh)

  \brief Explicit predictor-corrector gamma = 0.5 and beta = 0.25
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
 */
void U_PCE(Mesh, GaussPoint);

#endif
