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

  int Ndim = NumberDimensions;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /* Control parameters of the generalized-alpha algorithm 
     all the parameters are controled by a simple parameter :
     SpectralRadius */
  Time_Int_Params Params;
  Params.GA_gamma = 0.5;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix Phi_I;

  Matrix M_I = MatAssign(1, Nnodes, NAN, NULL, NULL);
  strcpy(M_I.Info,"MASS");

  Matrix P_I =
    MatAssign(Ndim,Nnodes,
	      NAN, NULL, (double **)malloc((unsigned)Ndim*sizeof(double*)));
  strcpy(P_I.Info,"MOMENTUM");

  Matrix V_I;
  Matrix F_I;
  
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
      PCE_Predictor(MPM_Mesh, FEM_Mesh, M_I,
		    P_I, F_I,Params);
      LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
      puts(" Second step : Reset nodal values");
      puts(" \t WORKING ...");
      FreeMat(Phi_I);
      FreeMat(V_I);
      FreeMat(F_I);
      puts(" DONE !!!");
      puts("*************************************************");
    }
    
    puts(" Third step : Get the nodal mass and the momentum");
    puts(" \t WORKING ...");
    Phi_I = compute_NodalFields(MPM_Mesh,FEM_Mesh);
    imposse_NodalMomentum(FEM_Mesh,Phi_I,TimeStep);
    puts(" \t DONE !!! \n");
    puts("*************************************************");
    
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,P_I,
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
    imposse_NodalMomentum(FEM_Mesh,F_I,TimeStep);
    puts(" \t DONE !!!");
    
    puts("*************************************************");
    puts(" Third step : Update nodal momentum ... WORKING");
    update_NodalMomentum(FEM_Mesh,Phi_I,F_I);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Seven step : Corrector stage");
    puts(" \t WORKING ...");
    PCE_Corrector(MPM_Mesh, FEM_Mesh, M_I, F_I, Params);
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");

  } /* End of temporal integration */
}
