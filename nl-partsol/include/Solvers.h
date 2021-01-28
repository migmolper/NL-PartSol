/*! \file Solvers.h
    \brief File with the prototype of time integration scheme
*/

#ifndef _SOLVERS_H_
#define _SOLVERS_H_


/*!
  \fn void NonLinear_Gauss_Point_Analysis(GaussPoint PointAnalysis)

  \brief Function to test advanced constitutive models

  \param GaussPoint PointAnalysis :  

*/
void NonLinear_Gauss_Point_Analysis(GaussPoint PointAnalysis);


/*!
  \fn double DeltaT_CFL(GaussPoint MPM_Mesh, double DeltaX)

  \brief  Get the time step using \cite Anderson_1987

  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param double DeltaX : Minimum mesh size
*/
double DeltaT_CFL(GaussPoint, double);

/*!
  \fn double DeltaT_Coussy__SolversLib__(GaussPoint MPM_Mesh, double h, double xi)

  \brief  Get the time step using \cite Coussy_1995

  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param double DeltaX : Minimum mesh size
  \param double xi : Parameter to include the tortuosity of flow
*/
double DeltaT_Coussy__SolversLib__(GaussPoint, double, double);

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

  \brief Explicit predictor-corrector gamma = 0.5
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Newmark_Predictor_Corrector(Mesh, GaussPoint, int);

/*!
  \fn void U_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Finite strains explicit predictor-corrector gamma = 0.5
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Newmark_Predictor_Corrector_Finite_Strains(Mesh, GaussPoint, int);

/*
  \fn void U_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Finite strains Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
void U_Newmark_beta_Finite_Strains(Mesh, GaussPoint, int);

/*
  \fn void U_Newmark_beta_Finite_Strains_BDB(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Finite strains Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
void U_Newmark_beta_Finite_Strains_BDB(Mesh, GaussPoint, int);

/*
  \fn void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief Discrete energy-momentum method. Implicit iterative solver proposed 
  by \cite Simo_and_Tarnow_1992, preserves linear and angular momentum.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void U_Discrete_Energy_Momentum(Mesh, GaussPoint, int);


/*
  \fn void upw_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, GaussPoint MPM_Mesh, int InitialStep)

  \brief u-pw formulation with one single set of material points. Explicit solver for finite strains based in the 
    Newmark Predictor-Corrector, see \cite Molinos_et_al_2021_CMAME.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param GaussPoint MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void upw_Newmark_Predictor_Corrector_Finite_Strains(Mesh, GaussPoint, int);


#endif
