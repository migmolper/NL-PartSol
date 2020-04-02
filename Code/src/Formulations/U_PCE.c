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
  Matrix Nodal_Mass;
  Matrix Nodal_Velocity;
  Matrix Nodal_Forces = MatAllocZ(N_dim,FEM_Mesh.NumNodesMesh);
  
  /*********************************************************************/
  /*********************************************************************/

  for(int TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);
    
    puts("*****************************************");
    puts(" First step : Predictor stage ... WORKING");
    Nodal_Mass = GetNodalMass(MPM_Mesh, FEM_Mesh);
    Nodal_Velocity = PCE_Predictor(MPM_Mesh,FEM_Mesh,
				   Nodal_Velocity, Nodal_Mass,
				   Params, DeltaTimeStep);
    BCC_Nod_VALUE(FEM_Mesh,Nodal_Velocity,TimeStep);
   
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_Velocity,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }    
    puts(" DONE !!!");
    
    puts("*********************************************************");        
    puts(" Second step : Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(MPM_Mesh, FEM_Mesh, Nodal_Velocity);
    puts(" DONE !!!");

    puts("**********************************************************");        
    puts(" Third step : Update the particle stress state ... WORKING");
    UpdateGaussPointStress(MPM_Mesh);
    ComputeDamage(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");

    puts("******************************************************");
    puts(" Four step : Calculate total forces forces ... WORKING");
    Nodal_Forces = GetNodalForces(MPM_Mesh,FEM_Mesh,TimeStep);
    BCC_Nod_VALUE(FEM_Mesh,Nodal_Forces,TimeStep);
    puts(" DONE !!!");
    
    puts("****************************************");
    puts(" Five step : Corrector stage ... WORKING");
    Nodal_Velocity = PCE_Corrector(FEM_Mesh,Nodal_Velocity,
				   Nodal_Forces, Nodal_Mass,
				   Params, DeltaTimeStep);
    puts(" DONE !!!"); 

    puts("***************************************************");
    puts(" Six step : Update particles lagrangian ... WORKING");
    PCE_Update_Lagrangian(MPM_Mesh, FEM_Mesh,
   			  Nodal_Mass, Nodal_Velocity,
    			  Nodal_Forces,DeltaTimeStep);
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");

    puts("********************************************************");
    puts(" Seven step : Reset nodal values of the mesh ... WORKING");
    FreeMat(Nodal_Mass);
    FreeMat(Nodal_Velocity);
    FreeMat(Nodal_Forces);
    puts(" DONE !!!");

  } /* End of temporal integration */
}
