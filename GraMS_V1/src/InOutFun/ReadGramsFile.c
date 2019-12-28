#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h> 
#include "../GRAMS/grams.h"


/***************************************/
/*********** GraMS Interface ***********/
/***************************************/

GaussPoint DefineParticleMesh(char * GramsFile)
/*
  GramsParticles (GID=MPM_Mesh.msh) {
  idx=0
  analysis=U
  shapefun=uGIMP
  material=0
  }
*/
{
  GaussPoint ParticleMesh;
  return ParticleMesh;  
}

Fields InitializeParticles(char * GramsFile)
/*
  GramsInit () {
  rho = 1540
  }
*/
{
  Fields FieldsParticles;
  return FieldsParticles;
}

/**********************************************************************/

/* SimParam DefineSimulationParameters(char * GramsFile) */
/* /\* */
/* GramsTime (Euler) { */
/* 	  cfl=0.1 */
/* 	  dtmax=2e-8 */
/* 	  Steps=6000 */
/* } */
/* *\/ */
/* { */
/*   SimParam Simulation; */
/*   return Simulation; */
/* } */

/* Outputs * DefineOutputs(char * GramsFile) */
/* /\* */
/* GramsOutput (istep=25, Particles=0) { */
/* 	v = Stress , Velocity */
/* } */
/* *\/ */
/* {} */


