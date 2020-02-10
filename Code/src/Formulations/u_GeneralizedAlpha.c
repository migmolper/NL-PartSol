#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

void u_GeneralizedAlpha(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
/*!
 * The generalized-alpha algorithm here implemented is analogous
 * to the one described in "Temporal and null-space filter for the
 * material point method". DOI : 10.1002/nme.6138.
 * Tran and Solowski
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
  Params.GA_alpha =
    (2*SpectralRadius-1)/(1+SpectralRadius);
  Params.GA_beta = (5-3*SpectralRadius)/
    (pow((1+SpectralRadius),2)*(2-SpectralRadius));
  Params.GA_gamma =
    3/2 - Params.GA_alpha;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the nodal kinetics
     Nodal_Kinetics = {mass, a0, a1, v}
   */
  int N_Nodes = FEM_Mesh.NumNodesMesh;
  int N_dim = NumberDimensions;
  Matrix Nodal_Velocity;
  Matrix Nodal_Kinetics;

  /* Nodal forces for the balance */
  Matrix Nodal_Forces = MatAssign(N_dim,N_Nodes,NAN,NULL,NULL);

  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);

    puts("*************************************************");
    puts(" First step : Get the nodal kinetics");
    puts(" \t WORKING ...");
    Nodal_Kinetics = GetNodalKinetics(MPM_Mesh,FEM_Mesh);
    Nodal_Velocity =
      MatAssign(N_dim,N_Nodes,NAN,NULL,
		(double**)malloc(NumberDimensions*sizeof(double *)));
    for(int i = 0 ; i<NumberDimensions ; i++){
      Nodal_Velocity.nM[i] = Nodal_Kinetics.nM[1+2*N_dim+i];
    }
    puts(" \t DONE !!! \n");
    
    puts("*************************************************");
    puts(" Second step : Set the essential BCC (over P)");
    puts(" \t WORKING ...");
    BCC_Nod_VALUE(FEM_Mesh,Nodal_Velocity,TimeStep);
    puts(" \t DONE !!!");
    
    /* Print nodal and Gps values */
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_Velocity,
      		   (int)TimeStep/ResultsTimeStep);
      /* Print GPs results */
      WriteVtk_MPM("MPM_VALUES",MPM_Mesh,List_Fields,
      		   (int)TimeStep/ResultsTimeStep);
    }
    
    puts("*************************************************");
    puts(" Third step : Update the particle stress state");
    puts(" \t DONE !!!");
    puts(" \t a) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(MPM_Mesh, FEM_Mesh, Nodal_Velocity);
    puts(" \t DONE !!!");
    puts(" \t b) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(MPM_Mesh);
    puts(" \t DONE !!!");
    
    puts("*************************************************");
    puts(" Four step : Calculate total forces forces");
    puts(" \t WORKING ...");
    Nodal_Forces = GetNodalForces(MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Five step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    GA_UpdateNodalKinetics(FEM_Mesh, Nodal_Kinetics, Nodal_Forces, Params);
    BCC_Nod_VALUE(FEM_Mesh, Nodal_Velocity, TimeStep);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Six step : Update the particle kinetics");
    puts(" \t WORKING ...");
    GA_AdvectionKinetics(MPM_Mesh, FEM_Mesh, Nodal_Kinetics, Params);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Seven step : Search the GP in the mesh");
    puts(" \t WORKING ...");
    LocalSearchGaussPoints(MPM_Mesh, FEM_Mesh);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Nine step : Reset nodal forces");
    puts(" \t WORKING ...");
    FreeMat(Nodal_Forces);
    free(Nodal_Velocity.nM);
    FreeMat(Nodal_Kinetics);
    puts(" DONE !!!");

  } /* End of temporal integration */
  
}
