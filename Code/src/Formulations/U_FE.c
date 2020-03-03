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

  int Ndim = 3;
  int Nnodes = FEM_Mesh.NumNodesMesh;

  /*********************************************************************/
  /***** INITIALIZE AUXILIAR STRUCTURES TO STORE NODAL INFORMATION *****/
  /*********************************************************************/

  /* Auxiliar variable for the mass and momentum */
  Matrix Phi_I;
  strcpy(Phi_I.Info,"MOMENTUM;MASS");
  Matrix V_I;
  Matrix F_I;
  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    DeltaTimeStep = DeltaT_CFL(MPM_Mesh, FEM_Mesh.DeltaX);
    printf("***************** STEP : %i , DeltaT : %f \n",
	   TimeStep,DeltaTimeStep);

    puts("*************************************************");
    puts(" First step : Get the nodal fields ... WORKING");
    puts(" \t Nodal mass and momentum");    
    Phi_I = compute_NodalFields(MPM_Mesh,FEM_Mesh);
    puts(" \t DONE !!!");
    puts(" \t Essential boundary conditions over P");
    imposse_NodalMomentum(FEM_Mesh,Phi_I,TimeStep);
    puts(" \t DONE !!!");
    puts(" \t Compute nodal velocity");
    V_I = compute_NodalVelocity(FEM_Mesh, Phi_I);
    puts(" \t DONE !!!");

    
    if(TimeStep % ResultsTimeStep == 0){
      /* Print Nodal values after appling the BCCs */
      WriteVtk_FEM("Mesh",FEM_Mesh,Phi_I,
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

    /* if(MPM_Mesh.Mat[0].Fracture){ */
    /*   printf(" \t d) %s %i %s \n", */
    /* 	     "Update particle status for material",0, */
    /* 	     "... WORKING"); */
    /*   UpdateBeps(MPM_Mesh,FEM_Mesh); */
    /*   MPM_Mesh.Phi.ji = */
    /* 	ComputeDamage(MPM_Mesh.Phi.ji, MPM_Mesh.Phi.W, MPM_Mesh.Phi.mass, */
    /* 		      MPM_Mesh.Phi.Stress, MPM_Mesh.MatIdx,MPM_Mesh.Mat, */
    /* 		      MPM_Mesh.Beps,FEM_Mesh.DeltaX); */
    /*   puts(" \t DONE !!!"); */
    /* } */

    
    puts("*************************************************");
    puts(" Third step : Update nodal momentum ... WORKING");
    update_NodalMomentum(FEM_Mesh,Phi_I,F_I);
    BCC_Nod_VALUE(FEM_Mesh,F_I,TimeStep);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Four step : Update lagrangian ... WORKING");
    UpdateVelocityAndPositionGP(MPM_Mesh, FEM_Mesh, Phi_I, F_I);
    LocalSearchGaussPoints(MPM_Mesh,FEM_Mesh);
    puts(" DONE !!!");
    
    puts("*************************************************");
    puts(" Five step : Reset nodal values ... WORKING");
    FreeMat(Phi_I);
    FreeMat(V_I);
    FreeMat(F_I);
    puts(" DONE !!!");

  } /* End of temporal integration */

}
