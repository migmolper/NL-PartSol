#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

void u_ForwardEuler(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*
  Displacement formulation of the MPM with a Forward Euler as 
  time integrator scheme
*/
{

  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;
  /* Time step */
  int TimeStep;
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = NumberDimensions;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix Nodal_MASS_MOMENTUM;

  Matrix Nodal_MASS = MatAssign(1,N_Nodes,NAN,NULL,NULL);
  strcpy(Nodal_MASS.Info,"MASS");

  Matrix Nodal_MOMENTUM = MatAssign(N_dim,N_Nodes,NAN,NULL,
				    (double **)malloc(N_dim*sizeof(double*)));
  strcpy(Nodal_MOMENTUM.Info,"MOMENTUM");

  Matrix Nodal_Velocity;
  strcpy(Nodal_Velocity.Info,"VELOCITY");
  
  Matrix Nodal_Forces;
  Matrix Nodal_Reactions;
  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);
    
    puts("*************************************************");
    puts(" First step : Get the nodal mass and the momentum");
    puts(" \t WORKING ...");
    Nodal_MASS_MOMENTUM = GetNodalMassMomentum(MPM_Mesh,FEM_Mesh);
    Nodal_MASS.nV = Nodal_MASS_MOMENTUM.nM[0];
    Nodal_MOMENTUM.nM[0] = Nodal_MASS_MOMENTUM.nM[1];
    Nodal_MOMENTUM.nM[1] = Nodal_MASS_MOMENTUM.nM[2];
    puts(" \t DONE !!! \n");
    
    puts("*************************************************");
    puts(" Second step : Set the essential BCC (over P)");
    puts(" \t WORKING ...");
    BCC_Nod_VALUE(FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    puts(" \t DONE !!!");
    
    puts("*************************************************");
    puts(" Third step : Update the particle stress state");
    puts(" \t a) Get the grid nodal velocity ... WORKING");
    Nodal_Velocity = GetNodalVelocity(FEM_Mesh,Nodal_MOMENTUM,Nodal_MASS);
    puts(" \t DONE !!!");
    
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_Velocity,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    } 
    
    puts(" \t b) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(MPM_Mesh,FEM_Mesh,Nodal_Velocity);
    puts(" \t DONE !!!");
    
    puts(" \t c) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(MPM_Mesh);
    ComputeDamage(MPM_Mesh,FEM_Mesh);
    puts(" \t DONE !!!");
    
    puts("*************************************************");
    puts(" Four step : Calculate total forces forces");
    puts(" \t WORKING ...");
    Nodal_Forces = GetNodalForces(MPM_Mesh,FEM_Mesh,TimeStep);
    Nodal_Reactions = GetNodalReactions(FEM_Mesh,Nodal_Forces);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Five step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    FE_Update_Momentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_Forces);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Six step : Update the particle velocity and position");
    puts(" \t WORKING ...");
    FE_Update_Lagrangian(MPM_Mesh,FEM_Mesh,Nodal_MASS,
			  Nodal_MOMENTUM,Nodal_Forces);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Seven step : Search the GP in the mesh");
    puts(" \t WORKING ...");
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Eight step : Reset nodal values of the mesh");
    puts(" \t WORKING ...");
    FreeMat(Nodal_MASS_MOMENTUM);
    FreeMat(Nodal_Velocity);
    FreeMat(Nodal_Forces);
    FreeMat(Nodal_Reactions);
    puts(" DONE !!!");

  } /* End of temporal integration */

}
