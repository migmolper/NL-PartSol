#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/GlobalVariables.h"
#include "InOutFun/InOutFun.h"
#include "Constitutive/Constitutive.h"
#include "ElementsFunctions/ElementTools.h"
#include "GaussPointsFunctions/GaussPointsTools.h"
#include "Solvers/Solvers.h"
#include "Num_Schemes/Num_Schemes.h"

void main(int argc, char *argv[])
/*
  Inputs parameters :
  * Mesh file
  * Data file
*/
{
  /* Check command-line arguments */
  if(argc == 1){
    perror("Error in main(), insuficient number of input files !");
    exit(0);
  }

  /* Read the .dat file */
  ReadDatFile(argv[1]);

  /* Read mesh data and initialize the element mesh */
  Element ElementMesh = ReadGidMesh(MeshFileName);
  ElementMesh.NumberDOF = 2; /* Change this */
    
  /* Read the initial conditions fields as a CSV */
  Matrix InputFields = Read_CSV(InitCondFileName,9);

  /* Read the boundary conditions */
  ReadBCC(BounCondFileName);
  
  /* Get material Constitutive matrix */
  Matrix D = LinearElastic1D(ElasticModulus);

  /* Define Gauss point mesh */
  GaussPoint GP_Mesh  = Initialize_GP_Mesh(InputFields,D,
					   ElementMesh);

  /* Define the fields of the sigma-v */
  Matrix Phi_n_GP; 
  /* Auxiliar pointer to the fields in t = n, in this special case, to avoid
     a second malloc that waste memory, we create a two rows table of pointer,
     so we have to fill the rest of the Matrix type field, in order to allow a 
     nide behaviour of the linear algebra functions */
  Phi_n_GP.nM =  (double **)malloc((unsigned)2*sizeof(double *));
  Phi_n_GP.nV = NULL; /* Set to NULL the (double *) pointer */
  Phi_n_GP.n = -999; /* Set to -999 the scalar variable */
  Phi_n_GP.N_rows = 2; /* Number of rows */
  Phi_n_GP.N_cols = GP_Mesh.NumGP; /* Number of columns */
  Phi_n_GP.nM[0] = GP_Mesh.Phi.Stress.nV; /* Asign to the first row the stress field */
  Phi_n_GP.nM[1] = GP_Mesh.Phi.vel.nV; /* Asign to the first row the velocity field */

  /* Define nodal values of the sigma-v */
  Matrix Phi_n_Nod = MatAlloc(2,ElementMesh.NumNodesMesh);
  
  Matrix A = MatAllocZ(2,2);
  A.nM[0][1] = - D.n;
  A.nM[1][0] = - (double)1/Density;
  printf("************************************************* \n");
  printf("Coupling matrix A : \n");
  printf("\t [ %f %f ] \n",A.nM[0][0],A.nM[0][1]);
  printf("\t [ %f %f ] \n",A.nM[1][0],A.nM[1][1]);

  /* printf("Phi_n_Nod \n %f ; %f \n %f ; %f \n", */
  /* 	   Phi_n_Nod.nM[0][0],Phi_n_Nod.nM[1][0], */
  /* 	   Phi_n_Nod.nM[0][1],Phi_n_Nod.nM[1][1]); */

  /* Get the lumped mass matrix, in this case, we dont have to update it */
  /* Matrix M = Get_Geom_Mass_Matrix(GP_Mesh,ElementMesh); */
  /* Matrix M_l = Get_Lumped_Matrix(M); */
  Matrix M_l = GetMassMatrix_L(ElementMesh,GP_Mesh);

  /* Transfer information from the GP to the mesh */
  GaussPointsToMesh(ElementMesh,GP_Mesh,Phi_n_GP,Phi_n_Nod,M_l);

  /* Run simulation  */
  printf("************************************************* \n");
  printf("Run simulation :  !!! \n");

  for(int t_i = 0 ; t_i<25; t_i++){
    printf("************************************************* \n");
    printf("Time step : %i \n",t_i);

    /* Print the solution in the GP */    
    GnuplotOutput1D(GP_Mesh.Phi.x_GC,
		    GP_Mesh.Phi.vel,
		    0.0, 8.0,
		    t_i,GP_Mesh.NumGP,
		    "Velocity");
    GnuplotOutput1D(GP_Mesh.Phi.x_GC,
		    GP_Mesh.Phi.Stress,
		    0.0, 8.0,
		    t_i,GP_Mesh.NumGP,
		    "Stress");
       
    /* Solve the Eulerian part of the problem in the mesh and get the results in the nodes */
    Two_Steps_TG(ElementMesh,GP_Mesh,Phi_n_Nod,Phi_n_GP,A,M_l,t_i);

    /* Transfer information from the mesh to the GP for MPM */
    MeshToGaussPoints(ElementMesh,GP_Mesh,Phi_n_Nod,Phi_n_GP,M_l);
  
  }
  
  exit(EXIT_SUCCESS);
  
  
}


