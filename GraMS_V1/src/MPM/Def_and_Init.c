#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*********************************************************************/

Mesh InitializeMesh(char * GDF){

  Mesh Back_Mesh;

  puts("*************************************************");
  puts(" Generate the MPM mesh");
  puts(" \t Defining background FEM mesh ...");
  Back_Mesh = ReadGidMesh(FEM_MeshFileName);
  puts(" \t DONE !!!");
  puts(" \t Searching neighbours elements for each node ...");
  GetNodalConnectivity(Back_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading FEM boundary conditions ...");
  Back_Mesh.Bounds = Set_FEM_BCC(GDF, Back_Mesh);
  puts(" \t DONE !!! ");

  return Back_Mesh;
}

/*********************************************************************/

GaussPoint Define_GP_Mesh(char * MPM_GID_MeshName,
			  double Density)
/*
  
*/
{
  /* Material point mesh (Gauss-Points) */
  Mesh MPM_GID_Mesh;
  GaussPoint MPM_Mesh;
  Matrix Poligon_Coordinates;
  ChainPtr Poligon_Connectivity;
  ChainPtr Vertex;
  int NumVertex;
  Matrix Poligon_Centroid;

  /* Screen message */
  printf("Begin of initialize the Gauss-Points mesh !!! \n");

  /* Read GP mesh */
  MPM_GID_Mesh = ReadGidMesh(MPM_GID_MeshName);

  /* The number of Gauss-Points is the same as the number of elements
   in the input mesh, because we set a GP in the middle of each element */
  MPM_Mesh.NumGP = MPM_GID_Mesh.NumElemMesh;
 
  /* Allocate fields */
  
  /* Index of the Element */
  MPM_Mesh.Element_id =
    (int *)Allocate_ArrayZ(MPM_Mesh.NumGP,sizeof(int));  

  /* A list with the number of tributary nodes of the GP */
  MPM_Mesh.NumberNodes =
    (int *)Allocate_ArrayZ(MPM_Mesh.NumGP,sizeof(int));
  /* A table of chains with the nodal connectivity of the GP */
  MPM_Mesh.ListNodes =
    (ChainPtr *)malloc(MPM_Mesh.NumGP*sizeof(ChainPtr));
  if(MPM_Mesh.ListNodes == NULL){
    puts("Error in Chain declaration");
    exit(0);
  }
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    MPM_Mesh.ListNodes[i] = NULL;  
  }
  
  /* Define the shapefunction employed by the GP */
  strcpy(MPM_Mesh.ShapeFunctionGP,ShapeFunctionGP);

  /* Coordinates of the GP (Global/Local)*/
  MPM_Mesh.Phi.x_GC = MatAllocZ(MPM_Mesh.NumGP,3);
  strcpy(MPM_Mesh.Phi.x_GC.Info,"Global Coordinates");

  /* Material parameters */
  /* Mass */
  MPM_Mesh.Mat.mass = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Mat.mass.Info,"Mass GP");
  /* Density */
  MPM_Mesh.Mat.rho = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Mat.rho.Info,"Density GP");
  /* Elastic Modulus */
  MPM_Mesh.Mat.E = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Mat.E.Info,"Elastic modulus GP");  
  /* Poisson ratio */
  MPM_Mesh.Mat.mu = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Mat.mu.Info,"Poisson ratio GP");
  /* Damage parameter (fracture) */
  MPM_Mesh.Mat.ji = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Mat.ji.Info,"Damage parameter GP");
  /* Normalizing constant (fracture) */
  MPM_Mesh.Mat.Ceps = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Mat.Ceps.Info,"Ceps fracture GP");
  
  /* Allocate vectorial/tensorial fields */
  switch(NumberDimensions){
  case 1 :
    /* Natural coordinates (Vectorial) */
    MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.x_EC.Info,"Element Coordinates GP");
    /* Displacement field (Vectorial) */
    MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.dis.Info,"Displacement field GP");
    /* Velocity field (Vectorial) */
    MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.vel.Info,"Velocity field GP");
    /* Acceleration field (Vectorial) */
    MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.acc.Info,"Acceleration field GP");
    /* Strain field (Tensor) */
    MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.Strain.Info,"Strain field GP");
    /* Stress field (Tensor) */
    MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.Stress.Info,"Stress field GP");
    /* Deformation Energy (Scalar) */
    MPM_Mesh.Phi.W = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.W.Info,"Deformation Energy GP");
    /* Lenght of the Voxel (Only GIMP) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      MPM_Mesh.lp = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
    }
    /* Lagrange Multipliers (Only LME ) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
      MPM_Mesh.lambda = MatAllocZ(MPM_Mesh.NumGP,1);
      strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
    }
    break;
  case 2 :
    /* Natural coordinates (Vectorial) */
    MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,2);
    strcpy(MPM_Mesh.Phi.x_EC.Info,"Element Coordinates GP");    
    /* Displacement field (Vectorial) */
    MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,2);
    strcpy(MPM_Mesh.Phi.dis.Info,"Displacement field GP");
    /* Velocity field (Vectorial) */
    MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,2);
    strcpy(MPM_Mesh.Phi.vel.Info,"Velocity field GP");
    /* Acceleration field (Vectorial) */
    MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,2);
    strcpy(MPM_Mesh.Phi.acc.Info,"Acceleration field GP");
    /* Strain field (Tensor) */
    MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.Strain.Info,"Strain field GP");
    /* Stress field (Tensor) */
    MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.Stress.Info,"Stress field GP");
    /* Deformation Energy (Scalar) */
    MPM_Mesh.Phi.W = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.W.Info,"Deformation Energy GP");
    /* Lenght of the Voxel (Only GIMP) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      MPM_Mesh.lp = MatAllocZ(MPM_Mesh.NumGP,2);
      strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
    }  
    /* Lagrange Multipliers (Only LME ) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
      MPM_Mesh.lambda = MatAllocZ(MPM_Mesh.NumGP,2);
      strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
    }
    break;
  case 3:
    /* Natural coordinates (Vectorial) */
    MPM_Mesh.Phi.x_EC = MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.x_EC.Info,"Element Coordinates GP");   
    /* Displacement field (Vectorial) */
    MPM_Mesh.Phi.dis = MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.dis.Info,"Displacement field GP");
    /* Velocity field (Vectorial) */
    MPM_Mesh.Phi.vel = MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.vel.Info,"Velocity field GP");
    /* Acceleration field (Vectorial) */
    MPM_Mesh.Phi.acc = MatAllocZ(MPM_Mesh.NumGP,3);
    strcpy(MPM_Mesh.Phi.acc.Info,"Acceleration field GP");
    /* Strain field (Tensor) */
    MPM_Mesh.Phi.Strain = MatAllocZ(MPM_Mesh.NumGP,9);
    strcpy(MPM_Mesh.Phi.Strain.Info,"Strain field GP");
    /* Stress field (Tensor) */
    MPM_Mesh.Phi.Stress = MatAllocZ(MPM_Mesh.NumGP,9);
    strcpy(MPM_Mesh.Phi.Stress.Info,"Stress field GP");
    /* Deformation Energy (Scalar) */
    MPM_Mesh.Phi.W = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.W.Info,"Deformation Energy GP");
    /* Lenght of the Voxel (Only GIMP) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      MPM_Mesh.lp = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
    }   
    /* Lagrange Multipliers (Only LME ) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
      MPM_Mesh.lambda = MatAllocZ(MPM_Mesh.NumGP,3);
      strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
    }
    break;
  default:
    puts("Error in Initialize_GP_Mesh() : Wrong number of dimensions");
    exit(0);
  }

  /* Fill geometrical properties of the GP mesh */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* Get the connectivity of the elements vertex */
    Poligon_Connectivity = MPM_GID_Mesh.Connectivity[i];
    
    /* Generate a matrix with the poligon coordinates */
    NumVertex = MPM_GID_Mesh.NumNodesElem[i];
    Poligon_Coordinates = MatAllocZ(NumVertex,NumberDimensions); 
    /* Initialize chain interator */
    Vertex = Poligon_Connectivity;
    /* Loop in the chain to fill the poligon */
    for(int k = 0, I_Vertex = 0;
	(k<NumVertex) || (Vertex != NULL);
	k++, Vertex = Vertex->next){
      I_Vertex = Vertex->I;
      for(int l = 0 ; l<NumberDimensions ; l++){
	Poligon_Coordinates.nM[k][l] = MPM_GID_Mesh.Coordinates.nM[I_Vertex][l];
      }
    }

    /* Free data */
    FreeChain(MPM_GID_Mesh.Connectivity[i]);

    /* Get the area (Poligon_Centroid.n) 
       and the position of the centroid (Poligon_Centroid.nV) */
    Poligon_Centroid = Centroid_Poligon(Poligon_Coordinates);
    /* Free data */
    FreeMat(Poligon_Coordinates);
    
    /* Assign the mass parameter */
    MPM_Mesh.Mat.mass.nV[i] = Poligon_Centroid.n*Density;
    /* Set the initial density */
    MPM_Mesh.Mat.rho.nV[i] = Density;
    /* Set the elastic modulus */
    MPM_Mesh.Mat.E.nV[i] = ElasticModulus;
    /* Set the poisson ratio */
    MPM_Mesh.Mat.mu.nV[i] = PoissonModulus;
    /* Damage parameter */
    MPM_Mesh.Mat.ji.nV[i] = 0.0;

    
    
    /* Get the coordinates of the centre */
    MPM_Mesh.Phi.x_GC.nM[i][0] = Poligon_Centroid.nV[0];
    MPM_Mesh.Phi.x_GC.nM[i][1] = Poligon_Centroid.nV[1];
    MPM_Mesh.Phi.x_GC.nM[i][2] = 0.0;         
    /* Local coordinates of the element */
    MPM_Mesh.Element_id[i] = -999;
    MPM_Mesh.NumberNodes[i] = 4;
    /* Location in the natural coordinates
       of the element (Init to zero) */
    MPM_Mesh.Phi.x_EC.nM[i][0] = 0.0;
    MPM_Mesh.Phi.x_EC.nM[i][1] = 0.0;
    /* Initial displacement */
    MPM_Mesh.Phi.dis.nM[i][0] = 0.0;
    MPM_Mesh.Phi.dis.nM[i][1] = 0.0;
    /* Initial Velocity */
    MPM_Mesh.Phi.vel.nM[i][0] = 0.0;
    MPM_Mesh.Phi.vel.nM[i][1] = 0.0;
    /* Initial Acceleration */
    MPM_Mesh.Phi.acc.nM[i][0] = 0.0;
    MPM_Mesh.Phi.acc.nM[i][1] = 0.0;
    /* Initial Stress */
    MPM_Mesh.Phi.Stress.nM[i][0] = 0.0;
    MPM_Mesh.Phi.Stress.nM[i][1] = 0.0;
    MPM_Mesh.Phi.Stress.nM[i][2] = 0.0;
    /* Initial Strain */
    MPM_Mesh.Phi.Strain.nM[i][0] = 0.0;
    MPM_Mesh.Phi.Strain.nM[i][1] = 0.0;
    MPM_Mesh.Phi.Strain.nM[i][2] = 0.0;
    /* Initial Energy of deformation */
    MPM_Mesh.Phi.W.nV[i] = 0.0;

    /* Voxel lenght (GIMP) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"uGIMP2D") == 0){
      for(int j =0 ; j<NumberDimensions ; j++){
	MPM_Mesh.lp.nM[i][j] =
	  0.5*pow(Poligon_Centroid.n,(double)1/NumberDimensions);
      }
    }
    /* Lagrange Multipliers (Only LME ) */
    if(strcmp(MPM_Mesh.ShapeFunctionGP,"LME") == 0){
      for(int j =0 ; j<NumberDimensions ; j++){
	MPM_Mesh.lambda.nM[i][j] = 0.0;
      }
    }
    MPM_Mesh.Gamma = 1.8;
    /* Free data */
    FreeMat(Poligon_Centroid);
    
  }

  /* Free the input data */
  FreeMat(MPM_GID_Mesh.Coordinates);
  free(MPM_GID_Mesh.Connectivity);
  free(MPM_GID_Mesh.ActiveNode);
  free(MPM_GID_Mesh.NumNeighbour);
  free(MPM_GID_Mesh.NodeNeighbour);
  
  /* Final messages */
  printf("End of initialize the Gauss-Points mesh !!! \n");
  
  return MPM_Mesh;
}


/*********************************************************************/

GaussPoint InitializeGP(char * GDF, Mesh FEM_Mesh){

  GaussPoint GP_Mesh;

  puts("*************************************************");
  puts(" Generate the MPM mesh");
  puts(" \t Defining MPM mesh of GPs ...");
  GP_Mesh = Define_GP_Mesh(MPM_MeshFileName,Density);
  puts(" \t DONE !!!");
  puts(" \t Constitutive library for GPs ...");
  GP_Mesh.D = Contitutive();
  puts(" \t DONE !!!");
  puts(" \t Searching GPs in the FEM mesh ...");
  GlobalSearchGaussPoints(GP_Mesh,FEM_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading MPM load cases ...");
  GP_Mesh.F = Read_MPM_LoadCase_ExtForces(GDF,GP_Mesh);
  GP_Mesh.B = Read_MPM_LoadCase_BodyForces(GDF,GP_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading MPM initial conditions ...");  
  Read_MPM_InitVal(GDF,GP_Mesh);
  puts(" \t DONE !!!");

  return GP_Mesh;

}