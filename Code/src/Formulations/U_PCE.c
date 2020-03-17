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
  Matrix Nodal_MASS;
  Matrix Nodal_VELOCITY = initialize_NodalVelocity(MPM_Mesh, FEM_Mesh);
  strcpy(Nodal_VELOCITY.Info,"VELOCITY");
  Matrix Nodal_TOT_FORCES = MatAllocZ(N_dim,FEM_Mesh.NumNodesMesh);
  
  /*********************************************************************/
  /*********************************************************************/

  for(int TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){


    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_VELOCITY,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);
    
    puts("*************************************************");
    puts(" First step : Predictor stage");
    puts(" \t WORKING ...");
    Nodal_MASS = GetNodalMass(MPM_Mesh, FEM_Mesh);
    Nodal_VELOCITY = PredictorNodalVelocity(FEM_Mesh, Nodal_VELOCITY,
					    Nodal_TOT_FORCES,
					    Nodal_MASS, Params,
					    DeltaTimeStep);
    BCC_Nod_VALUE(FEM_Mesh,Nodal_VELOCITY,TimeStep);
    FreeMat(Nodal_TOT_FORCES);
    
    puts(" DONE !!!");
    puts("*************************************************");
        
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
    puts(" Seven step : Corrector stage");
    puts(" \t WORKING ...");
    Nodal_VELOCITY = CorrectorNodalVelocity(FEM_Mesh,Nodal_VELOCITY,
					    Nodal_TOT_FORCES,
					    Nodal_MASS,Params,
					    DeltaTimeStep);    
    Update_Lagrangian_PCE(MPM_Mesh, FEM_Mesh,
			  Nodal_MASS, Nodal_VELOCITY,
			  Nodal_TOT_FORCES);
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    FreeMat(Nodal_MASS);
    puts(" DONE !!!");

  } /* End of temporal integration */
}
