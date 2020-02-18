#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

void U_PCE(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*!
 * Explicit predictor-corrector gamma = 0.5 and beta = 0.25
 */
{

  /* Some auxiliar variables for the outputs */
  Matrix List_Fields;

  int N_dim = NumberDimensions;

  /* Control parameters of the generalized-alpha algorithm 
     all the parameters are controled by a simple parameter :
     SpectralRadius */
  Time_Int_Params Params;
  Params.GA_gamma = 0.5;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix Nodal_MASS_MOMENTUM;

  Matrix Nodal_MASS = MatAssign(1, FEM_Mesh.NumNodesMesh, NAN, NULL, NULL);
  strcpy(Nodal_MASS.Info,"MASS");

  Matrix Nodal_MOMENTUM =
    MatAssign(N_dim,FEM_Mesh.NumNodesMesh,
	      NAN, NULL, (double **)malloc((unsigned)N_dim*sizeof(double*)));
  strcpy(Nodal_MOMENTUM.Info,"MOMENTUM");

  Matrix Nodal_VELOCITY;
  Matrix Nodal_TOT_FORCES;
  
  /*********************************************************************/
  /*********************************************************************/

  for(int TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);
    
    if(TimeStep > 0){
      puts("*************************************************");
      puts(" First step : Predictor stage");
      puts(" \t WORKING ...");
      PCE_Predictor(MPM_Mesh, FEM_Mesh, Nodal_MASS,
		    Nodal_MOMENTUM, Nodal_TOT_FORCES,Params);
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      puts(" Second step : Reset nodal values");
      puts(" \t WORKING ...");
      FreeMat(Nodal_MASS_MOMENTUM);
      FreeMat(Nodal_VELOCITY);
      FreeMat(Nodal_TOT_FORCES);
      puts(" DONE !!!");
      puts("*************************************************");
    }
    
    puts(" Third step : Get the nodal mass and the momentum");
    puts(" \t WORKING ...");
    Nodal_MASS_MOMENTUM = GetNodalMassMomentum(MPM_Mesh,FEM_Mesh);
    Nodal_MASS.nV = Nodal_MASS_MOMENTUM.nM[0];
    Nodal_MOMENTUM.nM[0] = Nodal_MASS_MOMENTUM.nM[1];
    Nodal_MOMENTUM.nM[1] = Nodal_MASS_MOMENTUM.nM[2];
    BCC_Nod_VALUE(FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    puts(" \t DONE !!! \n");
    puts("*************************************************");
    
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_MOMENTUM,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }
    
    puts("*************************************************");
    puts(" Four step : Update the particle stress state");
    puts(" \t a) Get the grid nodal velocity ... WORKING");
    Nodal_VELOCITY = GetNodalVelocity(FEM_Mesh, Nodal_MOMENTUM, Nodal_MASS);
    puts(" \t DONE !!!");
    puts(" \t b) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(MPM_Mesh, FEM_Mesh, Nodal_VELOCITY);
    puts(" \t DONE !!!");
    puts(" \t c) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(MPM_Mesh);
    puts(" \t DONE !!!");

    puts("*************************************************");
    puts(" Five step : Calculate total forces forces");
    puts(" \t WORKING ...");
    Nodal_TOT_FORCES = GetNodalForces(MPM_Mesh,FEM_Mesh,TimeStep);
    BCC_Nod_VALUE(FEM_Mesh,Nodal_TOT_FORCES,TimeStep);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Six step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    UpdateGridNodalMomentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_TOT_FORCES);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Seven step : Corrector stage");
    puts(" \t WORKING ...");
    PCE_Corrector(MPM_Mesh, FEM_Mesh, Nodal_MASS, Nodal_TOT_FORCES, Params);
    puts(" DONE !!!");

  } /* End of temporal integration */
}
