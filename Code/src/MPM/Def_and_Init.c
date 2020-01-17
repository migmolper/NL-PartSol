#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/*********************************************************************/

GaussPoint GramsSolid2D(char * GDF, char * FileName)
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
  MPM_GID_Mesh = ReadGidMesh(FileName);

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
    printf("%s : %s \n",
	   "Define_GP_Mesh",
	   "Memory error for ListNodes");
    exit(0);
  }
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    MPM_Mesh.ListNodes[i] = NULL;  
  }
  /* Table with */
  MPM_Mesh.Beps =
    (ChainPtr *)malloc(MPM_Mesh.NumGP*sizeof(ChainPtr));
  if(MPM_Mesh.Beps == NULL){
    printf("%s : %s \n",
	   "Define_GP_Mesh",
	   "Memory error for Beps");
    exit(0);
  }
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    MPM_Mesh.Beps[i] = NULL;  
  }
  
  /* Coordinates of the GP (Global/Local)*/
  MPM_Mesh.Phi.x_GC = MatAllocZ(MPM_Mesh.NumGP,3);
  strcpy(MPM_Mesh.Phi.x_GC.Info,"Global Coordinates");

  /* Read Material parameters */
  MPM_Mesh.Mat = GramsMaterials(GDF,MPM_Mesh);
  MPM_Mesh.MatIdx = (int *)malloc(MPM_Mesh.NumGP*sizeof(int));
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    MPM_Mesh.MatIdx[i] = 0;
  }

  /* Read Shape functions parameters */
  GramsShapeFun(GDF);

  /* Lenght of the Voxel (Only GIMP) */
  if(strcmp(ShapeFunctionGP,"uGIMP") == 0){
    MPM_Mesh.lp = MatAllocZ(MPM_Mesh.NumGP,NumberDimensions);
    strcpy(MPM_Mesh.lp.Info,"Voxel lenght GP");
  }
  /* Lagrange Multipliers / Beta (Only LME ) */
  if(strcmp(ShapeFunctionGP,"LME") == 0){
    MPM_Mesh.lambda = MatAllocZ(MPM_Mesh.NumGP,NumberDimensions);
    strcpy(MPM_Mesh.lambda.Info,"Lagrange Multiplier");
    MPM_Mesh.Beta = MatAllocZ(MPM_Mesh.NumGP,NumberDimensions);
    strcpy(MPM_Mesh.Beta.Info,"Beta parameter");
  }
    
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
    /* Damage parameter (fracture) */
    MPM_Mesh.Phi.ji = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.ji.Info,"Damage parameter GP");
    /* Mass */
    MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.mass.Info,"Mass GP");
    /* Density */
    MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.rho.Info,"Density GP");
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
    /* Damage parameter (fracture) */
    MPM_Mesh.Phi.ji = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.ji.Info,"Damage parameter GP");
    /* Mass */
    MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.mass.Info,"Mass GP");
    /* Density */
    MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.rho.Info,"Density GP");
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
    /* Damage parameter (fracture) */
    MPM_Mesh.Phi.ji = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.ji.Info,"Damage parameter GP");
    /* Mass */
    MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.mass.Info,"Mass GP");
    /* Density */
    MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
    strcpy(MPM_Mesh.Phi.rho.Info,"Density GP");
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
    MPM_Mesh.Phi.mass.nV[i] = Poligon_Centroid.n*MPM_Mesh.Mat[0].rho;
    /* Set the initial density */
    MPM_Mesh.Phi.rho.nV[i] = MPM_Mesh.Mat[0].rho;
    
    /* Get the coordinates of the centre */
    MPM_Mesh.Phi.x_GC.nM[i][0] = Poligon_Centroid.nV[0];
    MPM_Mesh.Phi.x_GC.nM[i][1] = Poligon_Centroid.nV[1];
    MPM_Mesh.Phi.x_GC.nM[i][2] = 0.0;         
    /* Local coordinates of the element */
    MPM_Mesh.Element_id[i] = -999;
    MPM_Mesh.NumberNodes[i] = 4;
   
    /* Free data */
    FreeMat(Poligon_Centroid);    
  }


  /* Read initial values */
  Read_MPM_InitVal(GDF,MPM_Mesh);


  /* Read external and internal forces */
  MPM_Mesh.F = Read_MPM_LoadCase_ExtForces(GDF,MPM_Mesh);
  MPM_Mesh.B = Read_MPM_LoadCase_BodyForces(GDF,MPM_Mesh);


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

/* /\*********************************************************************\/ */

/* GaussPoint InitializeGP(char * GDF, Mesh FEM_Mesh){ */

/*   GaussPoint MPM_Mesh; */
  
/*   puts("*************************************************"); */
/*   puts(" Generate the MPM mesh"); */
/*   puts(" \t Defining MPM mesh of GPs ..."); */
/*   MPM_Mesh = GramsSolid2D(GDF,MPM_MeshFileName); */
/*   puts(" \t DONE !!!"); */
/*   /\* puts(" \t Reading MPM load cases ..."); *\/ */
/*   /\* MPM_Mesh.F = Read_MPM_LoadCase_ExtForces(GDF,MPM_Mesh); *\/ */
/*   /\* MPM_Mesh.B = Read_MPM_LoadCase_BodyForces(GDF,MPM_Mesh); *\/ */
/*   /\* puts(" \t DONE !!!"); *\/ */
/*   /\* puts(" \t Reading MPM initial conditions ...");   *\/ */
/*   /\* Read_MPM_InitVal(GDF,MPM_Mesh); *\/ */
/*   /\* puts(" \t DONE !!!"); *\/ */
/*   puts(" \t Initialize shape-functions parameters ..."); */
/*   if(strcmp(ShapeFunctionGP,"MPMQ4") == 0){ */
/*     Q4_Initialize(MPM_Mesh, FEM_Mesh); */
/*   } */
/*   if(strcmp(ShapeFunctionGP,"uGIMP") == 0){ */
/*     GIMP_Initialize(MPM_Mesh,FEM_Mesh); */
/*   } */
/*   if(strcmp(ShapeFunctionGP,"LME") == 0){ */
/*     LME_Initialize(MPM_Mesh,FEM_Mesh); */
/*   } */
/*   puts(" \t DONE !!!"); */
  
/*   return MPM_Mesh; */

/* } */
