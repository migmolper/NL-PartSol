/*! \file Solvers.h
    \brief File with the prototype of time integration scheme
*/

#ifndef _SOLVERS_H_
#define _SOLVERS_H_


/*!
  \fn double DeltaT_CFL(GaussPoint MPM_Mesh, double DeltaX)

  \brief  Get the time step using \cite Anderson_1987

  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param double DeltaX : Minimum mesh size
*/
double DeltaT_CFL(GaussPoint, double);

/*!
  \fn void U_Forward_Euler(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Displacement formulation of the MPM with a Forward Euler as 
  time integrator scheme. The algorithm was taken from \cite Zhang_book_2016

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Forward_Euler(Mesh, GaussPoint, int);


/*!
  \fn void SV_Two_Steps_Taylor_Galerkin(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief The stress-velocity formulation with a two steps Taylor Galerkin time intergration

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void SV_Two_Steps_Taylor_Galerkin(Mesh, GaussPoint, int);

/*!
  \fn void U_GA(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief The generalized-alpha algorithm here implemented is analogous
  to the one described in \cite Tran2019e

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_GA(Mesh, GaussPoint, int);

/*!
  \fn void U_Newmark_Predictor_Corrector(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Explicit predictor-corrector gamma = 0.5 and beta = 0.25
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Newmark_Predictor_Corrector(Mesh, GaussPoint, int);

/*
  \fn void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Discrete energy-momentum method. Implicit iterative solver proposed 
  by \cite Simo_and_Tarnow_1992 , preserves linear and angular momentum.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void U_Discrete_Energy_Momentum(Mesh, GaussPoint, int);

#endif
