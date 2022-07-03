
#ifndef _COURANT_H_
#define _COURANT_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include "Macros.h"
#include "Types.h"
#include "Matlib.h"

/*
  Call global variables
*/
//Mixture *Soil_Water_Mixtures;
//int Number_Soil_Water_Mixtures;


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

#endif