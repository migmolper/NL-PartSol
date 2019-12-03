#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../InOutFun/InOutFun.h"
#include "../MeshTools/MeshTools.h"
#include "../MPM_Subroutines/MPM_Subroutines.h"


void u_ForwardEuler(Mesh FEM_Mesh, GaussPoint GP_Mesh)
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

  /* Nodal mass */
  Matrix Nodal_MASS; 
  Nodal_MASS.N_rows = 1;
  Nodal_MASS.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MASS.nM = NULL;
  strcpy(Nodal_MASS.Info,"MASS");

  /* Nodal momentum */
  Matrix Nodal_MOMENTUM; 
  Nodal_MOMENTUM.N_rows = NumberDimensions;
  Nodal_MOMENTUM.N_cols = FEM_Mesh.NumNodesMesh;
  Nodal_MOMENTUM.nV = NULL;
  Nodal_MOMENTUM.nM =
    (double **)malloc((unsigned)NumberDimensions*sizeof(double*));
  strcpy(Nodal_MOMENTUM.Info,"MOMENTUM");

  /* Nodal velocities */
  Matrix Nodal_VELOCITY;
  
  /* Nodal total forces */
  Matrix Nodal_TOT_FORCES;
  
  /*********************************************************************/
  /*********************************************************************/

  for(TimeStep = 0 ; TimeStep<NumTimeStep ; TimeStep++ ){

    puts("*************************************************");
    printf("********************** STEP : %i \n",TimeStep);
    puts("*************************************************");
    
    /* First step : Output Gauss-Points values to Paraview */
    if(TimeStep % ResultsTimeStep == 0){
      puts("*************************************************");
      puts(" First step : Output Gauss-Points values to Paraview");
      puts(" \t WORKING ...");
      WriteVtk_MPM("MPM_VALUES",GP_Mesh,List_Fields,
		   (int)TimeStep/ResultsTimeStep);
      printf(" \t DONE !!! \n");
    }
    
    /* Second step : Get the nodal mass and the momentum */
    puts("*************************************************");
    puts(" Second step : Get the nodal mass and the momentum");
    puts(" \t WORKING ...");
    Nodal_MASS_MOMENTUM = GetNodalMassMomentum(GP_Mesh,FEM_Mesh);
    Nodal_MASS.nV = Nodal_MASS_MOMENTUM.nM[0];
    Nodal_MOMENTUM.nM[0] = Nodal_MASS_MOMENTUM.nM[1];
    Nodal_MOMENTUM.nM[1] = Nodal_MASS_MOMENTUM.nM[2];
    printf(" \t DONE !!! \n");

    /* Third step : Set the essential boundary conditions (over p)*/
    puts("*************************************************");
    puts(" Third step : Set the essential BCC (over P)");
    puts(" \t WORKING ...");
    BCC_Nod_VALUE(FEM_Mesh,Nodal_MOMENTUM,TimeStep);
    puts(" DONE !!!");

    /* Four step : Update the particle stress state */
    puts("*************************************************");
    puts(" Four step : Update the particle stress state");
  
    /* a) Get the grid nodal velocity */
    puts(" \t a) Get the grid nodal velocity ... WORKING");
    Nodal_VELOCITY = GetNodalVelocity(FEM_Mesh,
				      Nodal_MOMENTUM,
				      Nodal_MASS);
    puts(" \t DONE !!!");
    /* b) Calculate the strain increment */
    puts(" \t b) Calculate the strain increment ... WORKING");
    UpdateGaussPointStrain(GP_Mesh,
			   FEM_Mesh,
			   Nodal_VELOCITY);
    puts(" \t DONE !!!");
    /* c) Update the particle stress state */
    puts(" \t c) Update the particle stress state ... WORKING");
    UpdateGaussPointStress(GP_Mesh);
    puts(" \t DONE !!!");

    /* Five step : Calculate total forces */
    puts("*************************************************");
    puts(" Five step : Calculate total forces forces");
    puts(" \t WORKING ...");
    /* BCC_GP_Forces(GP_Mesh, BCC_Loads, TimeStep); */
    Nodal_TOT_FORCES = GetNodalForces(GP_Mesh,FEM_Mesh,TimeStep);
    puts(" DONE !!!");    

    /* Six step : Integrate the grid nodal momentum equation */
    puts("*************************************************");
    puts(" Six step : Integrate the grid nodal momentum equation");
    puts(" \t WORKING ...");
    UpdateGridNodalMomentum(FEM_Mesh,Nodal_MOMENTUM,Nodal_TOT_FORCES);
    BCC_Nod_VALUE(FEM_Mesh,Nodal_TOT_FORCES,TimeStep);
    puts(" DONE !!!");

    /* Seven step : Update the particle velocity and position */
    puts("*************************************************");
    puts(" Seven step : Update the particle velocity and position");
    puts(" \t WORKING ...");
    UpdateVelocityAndPositionGP(GP_Mesh,FEM_Mesh,Nodal_MASS,
				Nodal_MOMENTUM,Nodal_TOT_FORCES);
    puts(" DONE !!!");

    
    /* Eight step : Search the GP in the mesh */
    puts("*************************************************");
    puts(" Eight step : Search the GP in the mesh");
    puts(" \t WORKING ...");
    LocalSearchGaussPoints(GP_Mesh,FEM_Mesh);
    puts(" DONE !!!");

    /* Nine step : Print nodal values */
    if(TimeStep % ResultsTimeStep == 0){
      puts("*************************************************");
      puts(" Nine step : Print nodal values");
      puts(" \t WORKING ...");
      WriteVtk_FEM("Mesh",FEM_Mesh,Nodal_MOMENTUM,(int)TimeStep/ResultsTimeStep);
      puts(" DONE !!!");
    }
    
    /* Ten step : Store all the material properties in the particles
       so that the deformed grid can be discarted */
    puts("*************************************************");
    puts(" Ten step : Reset nodal values of the mesh");
    puts(" \t WORKING ...");
    /* Reset nodal values */
    FreeMat(Nodal_MASS_MOMENTUM);
    FreeMat(Nodal_VELOCITY);
    FreeMat(Nodal_TOT_FORCES);
    puts(" DONE !!!");

  } /* End of temporal integration */

}
