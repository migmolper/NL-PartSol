#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grams.h"

void U_GA(Mesh FEM_Mesh, GaussPoint MPM_Mesh)
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
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Ndim = NumberDimensions;
  Matrix V_I;
  Matrix Nodal_Kinetics;

  /* Nodal forces for the balance */
  Matrix F_I = MatAssign(Ndim,Nnodes,NAN,NULL,NULL);

  
  puts("*************************************************");
  puts(" First step : Get the nodal kinetics");
  puts(" \t WORKING ...");
  Nodal_Kinetics = GetNodalKinetics(MPM_Mesh,FEM_Mesh);
  V_I =  MatAssign(Nnodes,Ndim,NAN,NULL,
		   (double**)malloc(Nnodes*sizeof(double *)));
  for(int i = 0 ; i<Ndim ; i++){
    V_I.nM[i] = Nodal_Kinetics.nM[1+2*Ndim+i];
  }
  puts(" \t DONE !!! \n");
  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);
    
    puts("*************************************************");
    puts(" First step : Essential BCC in P ... WORKING");
    /* BCC_Nod_VALUE(FEM_Mesh,V_I,TimeStep); */
    puts(" \t DONE !!!");

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
    puts(" Second step : Compute internal forces ... WORKING");
    F_I = MatAllocZ(Nnodes,Ndim);
    puts(" \t Compute internal forces");
    F_I = compute_InternalForces(F_I, V_I, MPM_Mesh, FEM_Mesh, DeltaTimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute body forces");
    F_I = compute_BodyForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute contact forces");
    F_I = compute_ContacForces(F_I, MPM_Mesh, FEM_Mesh, TimeStep);
    puts(" \t DONE !!!");
    
    puts("*************************************************");
    puts(" Third step : Update nodal momentum ... WORKING");
    GA_UpdateNodalKinetics(FEM_Mesh, Nodal_Kinetics, F_I, Params);
    /* BCC_Nod_VALUE(FEM_Mesh, V_I, TimeStep); */
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Four step : Update lagrangian ... WORKING");
    update_Particles_GA(MPM_Mesh, FEM_Mesh, Nodal_Kinetics, Params);
    LocalSearchGaussPoints(MPM_Mesh, FEM_Mesh);
    puts(" DONE !!!");
        
    puts("*************************************************");
    puts(" Five step : Reset nodal forces ... WORKING");
    FreeMat(F_I);
    puts(" DONE !!!");

  } /* End of temporal integration */
  
}
