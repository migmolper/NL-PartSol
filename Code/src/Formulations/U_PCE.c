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
  /* Time step */
  int TimeStep;

  /* Control parameters of the generalized-alpha algorithm 
   all the parameters are controled by a simple parameter :
   SpectralRadius */
  Time_Int_Params Params;
  Params.GA_gamma = 0.5;
  Params.GA_beta = 0.25;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the nodal kinetics
     Nodal_Kinetics = {mass, a0, a1, v}
   */
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = NumberDimensions;
  Matrix Nodal_Velocity =
    MatAssign(N_dim,N_Nodes,NAN,NULL,
	      (double**)malloc(NumberDimensions*sizeof(double *)));
  Matrix Nodal_Kinetics;

  /* Nodal forces for the balance */
  Matrix Nodal_Forces = MatAssign(N_dim,N_Nodes,NAN,NULL,NULL);

  puts("*************************************************");
  puts(" Get the nodal kinetics");
  puts(" \t WORKING ...");
  Nodal_Kinetics = GetNodalVelocityDisplacement(MPM_Mesh, FEM_Mesh);
  for(int i = 0 ; i<NumberDimensions ; i++){
    Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+N_dim+i];
  }
  puts(" \t DONE !!! \n");

  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);

    /* Print nodal and Gps values */
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_Kinetics,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }
    
    puts("*************************************************");
    puts("First step : Set the essential BCC (over P)");
    puts("\t WORKING ...");
    BCC_Nod_VALUE(FEM_Mesh,Nodal_Velocity,TimeStep);
    puts(" \t DONE !!!");

    puts("*************************************************");
    puts("Second  step : Predictor step");
    puts("\t WORKING ...");
    PCE_Predictor(FEM_Mesh, Nodal_Kinetics, Params);
    
    puts("*************************************************");
    puts("Third step : Calculate total forces forces");
    puts("\t a) Calculate the strain increment ... WORKING");
    puts("\t WORKING ...");
    UpdateGaussPointStrain(MPM_Mesh, FEM_Mesh, Nodal_Velocity);
    puts("\t DONE !!!");
    puts("\t b) Update the particle stress state ... WORKING");
    puts("\t WORKING ...");
    UpdateGaussPointStress(MPM_Mesh);
    puts("\t DONE !!!");
    puts("\t c) Evaluate total forces vector ... WORKING");
    puts("\t WORKING ...");
    Nodal_Forces = GetNodalForces(MPM_Mesh, FEM_Mesh, TimeStep);
    puts("DONE !!!");
    
    puts("*************************************************");
    puts("Four step : Corrector step");
    puts("\t WORKING ...");
    PCE_Corrector(FEM_Mesh, Nodal_Kinetics, Nodal_Forces, Params);
    BCC_Nod_VALUE(FEM_Mesh, Nodal_Velocity, TimeStep);
    puts("DONE !!!");
    
    puts("*************************************************");
    puts("Five step : Update the particle kinetics");
    puts("\t WORKING ...");
    PCE_AdvectionKinetics(MPM_Mesh, FEM_Mesh, Nodal_Kinetics);
    puts("DONE !!!");
    
    puts("*************************************************");
    puts("Six step : Search the GP in the mesh");
    puts("\t WORKING ...");
    LocalSearchGaussPoints(MPM_Mesh, FEM_Mesh);
    puts("DONE !!!");
    
    puts("*************************************************");
    puts("Seven step : Reset nodal forces");
    puts("\t WORKING ...");
    FreeMat(Nodal_Forces);
    puts("DONE !!!");

  } /* End of temporal integration */
  
}
