#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../MathTools/MathTools.h"
#include "../InOutFun/InOutFun.h"
#include "../MeshTools/MeshTools.h"
#include "../Constitutive/Constitutive.h"
#include "../MPM_Subroutines/MPM_Subroutines.h"
#include "GaussPointsTools.h"

/*********************************************************************/

GaussPoint Define_GP_Mesh(char * MPM_GID_MeshName,
			  double Density)
/*
  
*/
{
  /* Material point mesh (Gauss-Points) */
  Mesh MPM_GID_Mesh;
  GaussPoint MPM_Mesh;

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

  /* Coordinates of the GP (Global/Local)*/
  MPM_Mesh.Phi.x_GC = MatAllocZ(MPM_Mesh.NumGP,3);
  strcpy(MPM_Mesh.Phi.x_GC.Info,"Global Coordinates");

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
    break;
  default:
    puts("Error in Initialize_GP_Mesh() : Wrong number of dimensions");
    exit(0);
  }  

  /* Mass field (Scalar) */
  MPM_Mesh.Phi.mass = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Phi.mass.Info,"Mass field GP");
 
  /* Density field (Scalar) */
  MPM_Mesh.Phi.rho = MatAllocZ(MPM_Mesh.NumGP,1);
  strcpy(MPM_Mesh.Phi.rho.Info,"Density field GP");

  
  Matrix Poligon = MatAllocZ(MPM_GID_Mesh.NumNodesElem,NumberDimensions);
  int * Poligon_Connectivity;  
  Matrix Poligon_Centroid;

  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* Get the connectivity of the elements vertex */
    Poligon_Connectivity = MPM_GID_Mesh.Connectivity[i];

    /* Get the coordiantes of the element vertex */
    for(int j = 0 ; j<MPM_GID_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumberDimensions ; k++){
	Poligon.nM[j][k] =
	  MPM_GID_Mesh.Coordinates.nM[Poligon_Connectivity[j]][k];
      }
    }

    /* Get the area (Poligon_Centroid.n) 
       and the position of the centroid (Poligon_Centroid.nV) */
    Poligon_Centroid = Centroid_Poligon(Poligon);
           
    /* Assign the mass parameter */
    MPM_Mesh.Phi.mass.nV[i] = Poligon_Centroid.n*Density;
    /* Set the initial density */
    MPM_Mesh.Phi.rho.nV[i] = Density;
    /* Get the coordinates of the centre */
    MPM_Mesh.Phi.x_GC.nM[i][0] = Poligon_Centroid.nV[0];
    MPM_Mesh.Phi.x_GC.nM[i][1] = Poligon_Centroid.nV[1];
    MPM_Mesh.Phi.x_GC.nM[i][2] = 0.0;

    /* Free data */
    FreeMat(Poligon_Centroid);
      
    /* Local coordinates of the element */
    MPM_Mesh.Element_id[i] = -999;
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

    
  }
  /* Initialize all the fields :
     - Mass of the material point
     - Position field (Vectorial) in global coordiantes 
     and in element coordinates :
     - Displacement, Velocity and acceleration field (Vectorial)
     - Stress and Strain fields (Tensorial)
     Note : It is not necessary to allocate memory...think about it ;)
  */

  /* Free the input data */
  FreeMat(MPM_GID_Mesh.Coordinates);
  free(MPM_GID_Mesh.Connectivity);
  free(MPM_GID_Mesh.ActiveNode);

  /* Final messages */
  printf("End of initialize the Gauss-Points mesh !!! \n");
  
  return MPM_Mesh;
}


/*********************************************************************/

GaussPoint InitializeGP(char * GDF, Mesh FEM_Mesh, Matrix D){

  GaussPoint GP_Mesh;

  puts("*************************************************");
  puts(" Generate the MPM mesh");
  puts(" \t Defining MPM mesh of GPs ...");
  GP_Mesh = Define_GP_Mesh(MPM_MeshFileName,Density);
  puts(" \t DONE !!!");
  puts(" \t Constitutive law for the GPs ...");
  GP_Mesh.D = D;
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
