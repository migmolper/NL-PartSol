/*! \file Solvers.h
    \brief File with the prototype of time integration scheme
*/

#ifndef _SOLVERS_H_
#define _SOLVERS_H_


/*!
  \fn void NonLinear_Gauss_Point_Analysis(Particle PointAnalysis)

  \brief Function to test advanced constitutive models

  \param Particle PointAnalysis :  

*/
void NonLinear_Gauss_Point_Analysis(Particle PointAnalysis);


/*!
  \fn double U_DeltaT__SolversLib__(Particle MPM_Mesh, double DeltaX, Time_Int_Params Parameters_Solver)

  \brief  Get the time step using \cite Anderson_1987

  \param MPM_Mesh : Variable with the particle information
  \param DeltaX : Minimum mesh size
  \param Parameters_Solver : Parameters
*/
double U_DeltaT__SolversLib__(Particle, double, Time_Int_Params);

/*!
  \fn double DeltaT_Coussy__SolversLib__(Particle MPM_Mesh, double h, double xi)

  \brief  Get the time step using \cite Coussy_1995

  \param Particle MPM_Mesh : Variable with the particle information
  \param double DeltaX : Minimum mesh size
  \param double xi : Parameter to include the tortuosity of flow
*/
double DeltaT_Coussy__SolversLib__(Particle, double, double, double);

/*!
  \fn void U_Forward_Euler(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Displacement formulation of the MPM with a Forward Euler as 
  time integrator scheme. The algorithm was taken from \cite Zhang_book_2016

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Forward_Euler(Mesh, Particle, Time_Int_Params);

/*!
  \fn void U_Generalized_alpha(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief The generalized-alpha algorithm here implemented is analogous
  to the one described in \cite Tran2019e

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Generalized_alpha(Mesh, Particle, Time_Int_Params);

/*!
  \fn void U_Newmark_Predictor_Corrector(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Explicit predictor-corrector gamma = 0.5
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Newmark_Predictor_Corrector(Mesh, Particle, Time_Int_Params);

/*!
  \fn void U_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Finite strains explicit predictor-corrector gamma = 0.5
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
 */
void U_Newmark_Predictor_Corrector_Finite_Strains(Mesh, Particle, Time_Int_Params);

/*
  \fn void U_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Finite strains Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
void U_Newmark_beta_Finite_Strains(Mesh, Particle, Time_Int_Params);


/*
  \fn void U_Newmark_beta_Finite_Strains_HPC(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Finite strains Newmark-beta with HPC techniques
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
void U_Newmark_beta_Finite_Strains_HPC(Mesh,Particle,Time_Int_Params);


/*
  \fn void Up_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Incompresssible formulation of the finite strain Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
void Up_Newmark_beta_Finite_Strains(Mesh,Particle,Time_Int_Params);

/*
  \fn void U_Newmark_beta_Finite_Strains_BDB(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Finite strains Newmark-beta
 
  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep
*/
void U_Newmark_beta_Finite_Strains_BDB(Mesh, Particle, Time_Int_Params);

/*
  \fn void U_Discrete_Energy_Momentum(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief Discrete energy-momentum method. Implicit iterative solver proposed 
  by \cite Simo_and_Tarnow_1992, preserves linear and angular momentum.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void U_Discrete_Energy_Momentum(Mesh, Particle, Time_Int_Params);

/*
  \fn void upw_Newmark_Predictor_Corrector_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief u-pw formulation with one single set of material points. Explicit solver for finite strains based in the 
    Newmark Predictor-Corrector, see \cite Molinos_et_al_2021_CMAME.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void upw_Newmark_Predictor_Corrector_Finite_Strains(Mesh, Particle, Time_Int_Params);

/*
  \fn void upw_Newmark_beta_Finite_Strains(Mesh FEM_Mesh, Particle MPM_Mesh, Time_Int_Params Parameters_Solver)

  \brief u-pw formulation with one single set of material points. Implicit solver for finite strains based in the 
    Newmark-beta, see \cite Molinos_et_al_2021_CMAME.

  \param Mesh FEM_Mesh : Variable with the nodal information
  \param Particle MPM_Mesh : Variable with the particle information
  \param InitialStep  
*/
void upw_Newmark_beta_Finite_Strains(Mesh, Particle, Time_Int_Params);

#endif
