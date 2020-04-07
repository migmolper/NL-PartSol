#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

void U_FE(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*
  Displacement formulation of the MPM with a Forward Euler as 
  time integrator scheme
*/
{

  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;
  /* Time step */
  int TimeStep;

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix Phi_I;
  strcpy(Phi_I.Info,"MOMENTUM;MASS");
  Matrix V_I;
  Matrix F_I;
  Matrix R_I;
  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);

    puts("*************************************************");
    puts(" First step : Get the nodal fields ... WORKING");
    puts(" \t Nodal mass and momentum");    
    Phi_I = compute_NodalMomentumMass(MPM_Mesh,FEM_Mesh);
    puts(" \t DONE !!!");
    puts(" \t Essential boundary conditions over P");
    imposse_NodalMomentum(FEM_Mesh,Phi_I,TimeStep);    
    puts(" \t DONE !!!");
    puts(" \t Compute nodal velocity");
    V_I = compute_NodalVelocity(FEM_Mesh, Phi_I);
    puts(" \t DONE !!!");

    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Phi_I,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }

    puts("*************************************************");
    puts(" Second step : Update local state ... WORKING");
    update_LocalState(V_I, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Third step : Compute internal forces ... WORKING");
    F_I = MatAllocZ(Nnodes,Ndim);    
    puts(" \t Compute internal forces");
    F_I = compute_InternalForces(F_I, MPM_Mesh, FEM_Mesh);    
    puts(" \t DONE !!!");
    puts(" \t Compute body forces");
    F_I = compute_BodyForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute contact forces");
    F_I = compute_ContacForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute reactions");
    R_I = compute_Reactions(FEM_Mesh, F_I);
    puts(" \t DONE !!!");

    puts("*************************************************");
    puts(" Four step : Update nodal momentum ... WORKING");
    update_NodalMomentum(FEM_Mesh,Phi_I,F_I);
    puts(" DONE !!!");    
    puts("*************************************************");
    puts(" Five step : Update lagrangian ... WORKING");
    update_Particles_FE(MPM_Mesh, FEM_Mesh, Phi_I, F_I, DeltaTimeStep);
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Six step : Reset nodal values ... WORKING");
    FreeMat(Phi_I);
    FreeMat(V_I);
    FreeMat(F_I);
    FreeMat(R_I);
    puts(" DONE !!!");

  } /* End of temporal integration */

}
