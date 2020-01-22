#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

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

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix Nodal_MASS_MOMENTUM;

  Matrix Nodal_MASS; 
  Nodal_MASS.N_rows = 1;
  Nodal_MASS.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MASS.nM = NULL;
  strcpy(Nodal_MASS.Info,"MASS");

  Matrix Nodal_MOMENTUM; 
  Nodal_MOMENTUM.N_rows = NumberDimensions;
  Nodal_MOMENTUM.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MOMENTUM.nV = NULL;
  Nodal_MOMENTUM.nM =
    (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  strcpy(Nodal_MOMENTUM.Info,"MOMENTUM");

  Matrix Nodal_VELOCITY;
  Matrix Nodal_TOT_FORCES;
  
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
    
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_MOMENTUM,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }
    puts("*************************************************");
    puts(" Third step : Update the particle stress state");
    puts(" \t a) Get the grid nodal velocity ... WORKING");
    Nodal_VELOCITY = GetNodalVelocity(FEM_Mesh,
				      Nodal_MOMENTUM,
				      Nodal_MASS);
    puts(" \t DONE !!!");
    puts(" \t b) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(MPM_Mesh,
			   FEM_Mesh,
			   Nodal_VELOCITY);
    puts(" \t DONE !!!");
    puts(" \t c) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(MPM_Mesh);
    puts(" \t DONE !!!");
    /* if(MPM_Mesh.Mat[0].Fracture){ */
    /*   printf(" \t d) %s %i %s \n", */
    /* 	     "Update particle status for material",0, */
    /* 	     "... WORKING"); */
    /*   UpdateBeps(MPM_Mesh,FEM_Mesh); */
    /*   MPM_Mesh.Phi.ji = ComputeDamage(MPM_Mesh.Phi.ji, MPM_Mesh.Phi.W, */
    /* 				      MPM_Mesh.Phi.mass,MPM_Mesh.MatIdx, */
    /* 				      MPM_Mesh.Mat,MPM_Mesh.Beps, */
    /* 				      FEM_Mesh.DeltaX); */
    /*   puts(" \t DONE !!!"); */
    /* }   */
    puts("*************************************************");
    puts(" Four step : Calculate total forces forces");
    puts(" \t WORKING ...");
    Nodal_TOT_FORCES = GetNodalForces(MPM_Mesh,FEM_Mesh,TimeStep);
    puts(" DONE !!!");    
    puts("*************************************************");
    puts(" Five step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    UpdateGridNodalMomentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_TOT_FORCES);
    BCC_Nod_VALUE(FEM_Mesh,Nodal_TOT_FORCES,TimeStep);
    puts(" DONE !!!");
    puts("*************************************************");
    puts(" Six step : Update the particle velocity and position");
    puts(" \t WORKING ...");
    UpdateVelocityAndPositionGP(MPM_Mesh,FEM_Mesh,Nodal_MASS,
				Nodal_MOMENTUM,Nodal_TOT_FORCES);
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
    FreeMat(Nodal_VELOCITY);
    FreeMat(Nodal_TOT_FORCES);
    puts(" DONE !!!");

  } /* End of temporal integration */

}
