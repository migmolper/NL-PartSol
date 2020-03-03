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
  Matrix Phi_I;

  Matrix M_I = MatAssign(1, FEM_Mesh.NumNodesMesh, NAN, NULL, NULL);
  strcpy(M_I.Info,"MASS");

  Matrix P_I =
    MatAssign(N_dim,FEM_Mesh.NumNodesMesh,
	      NAN, NULL, (double **)malloc((unsigned)N_dim*sizeof(double*)));
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
    Phi_I = GetNodalMassMomentum(MPM_Mesh,FEM_Mesh);
    M_I.nV = Phi_I.nM[0];
    P_I.nM[0] = Phi_I.nM[1];
    P_I.nM[1] = Phi_I.nM[2];
    BCC_Nod_VALUE(FEM_Mesh,P_I,TimeStep);
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
    puts(" Four step : Update the particle stress state");
    puts(" \t a) Get the grid nodal velocity ... WORKING");
    V_I = GetNodalVelocity(FEM_Mesh, P_I, M_I);
    puts(" \t DONE !!!");
    puts(" \t b) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(MPM_Mesh, FEM_Mesh, V_I);
    puts(" \t DONE !!!");
    puts(" \t c) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(MPM_Mesh);
    puts(" \t DONE !!!");

    puts("*************************************************");
    puts(" Five step : Calculate total forces forces");
    puts(" \t WORKING ...");
    F_I = GetNodalForces(MPM_Mesh,FEM_Mesh,TimeStep);
    BCC_Nod_VALUE(FEM_Mesh,F_I,TimeStep);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Six step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    UpdateGridNodalMomentum(FEM_Mesh,P_I,F_I);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Seven step : Corrector stage");
    puts(" \t WORKING ...");
    PCE_Corrector(MPM_Mesh, FEM_Mesh, M_I, F_I, Params);
    puts(" DONE !!!");

  } /* End of temporal integration */
}
