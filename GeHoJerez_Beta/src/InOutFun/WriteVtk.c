
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../ToolsLib/Utils.h"
#include "InOutFun.h"

void WriteVtk_MPM(char * Name_File, GaussPoint MPM_Mesh,
		  Matrix List_Fields, int TimeStep_i){

  FILE * Vtk_file;

  char Name_file_t[80];

  sprintf(Name_file_t,"%s/MPM_%s_%i.vtk",OutputDir,Name_File,TimeStep_i);
  
  Vtk_file = fopen(Name_file_t,"w");

  /* Header */
  fprintf(Vtk_file,"# vtk DataFile Version 3.0 \n");
  fprintf(Vtk_file,"vtk output \n");
  fprintf(Vtk_file,"ASCII \n");
  fprintf(Vtk_file,"DATASET UNSTRUCTURED_GRID \n");

  /* Coordinates */
  fprintf(Vtk_file,"POINTS %i float \n",MPM_Mesh.NumGP);
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%f %f %f \n",
	    MPM_Mesh.Phi.x_GC.nM[i][0],
	    MPM_Mesh.Phi.x_GC.nM[i][1],
	    MPM_Mesh.Phi.x_GC.nM[i][2]);
  }

  /* Connectivity */
  fprintf(Vtk_file,"CELLS %i %i \n",MPM_Mesh.NumGP,2*MPM_Mesh.NumGP);
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i %i \n",1,i);
  }

  /* Type of element */
  fprintf(Vtk_file,"CELL_TYPES %i \n",MPM_Mesh.NumGP);
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i \n",1);
  }

  /* Point data */
  fprintf(Vtk_file,"POINT_DATA %i \n",MPM_Mesh.NumGP);
  fprintf(Vtk_file,"VECTORS %s float \n","X_GC");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%f %f %f \n",
	    MPM_Mesh.Phi.x_GC.nM[i][0],
	    MPM_Mesh.Phi.x_GC.nM[i][1],
	    MPM_Mesh.Phi.x_GC.nM[i][2]);
  }
  fprintf(Vtk_file,"VECTORS %s float \n","X_EC");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%f %f %f \n",
	    MPM_Mesh.Phi.x_EC.nM[i][0],
	    MPM_Mesh.Phi.x_EC.nM[i][1],
	    MPM_Mesh.Phi.x_EC.nM[i][2]);
  } 

  /* Cell data */  
  fprintf(Vtk_file,"CELL_DATA %i \n",MPM_Mesh.NumGP);

  fprintf(Vtk_file,"SCALARS MASS float \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%f \n",MPM_Mesh.Phi.mass.nV[i]);
  }

  fprintf(Vtk_file,"SCALARS MASS float \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%f \n",MPM_Mesh.Phi.mass.nV[i]);
  }

  fprintf(Vtk_file,"SCALARS DENSITY float \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%f \n",MPM_Mesh.Phi.rho.nV[i]);
    
  }
  
  fprintf(Vtk_file,"SCALARS ELEM_i float \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<MPM_Mesh.NumGP ; i++){
    fprintf(Vtk_file,"%i \n",MPM_Mesh.Element_id[i]);
    
  }
  
}

/*********************************************************************/

void WriteVtk_FEM(char * Name_File, Mesh ElementMesh,
		  Matrix List_Nod_Fields, int TimeStep_i){

  FILE * Vtk_file;

  char Name_file_t[80];

  sprintf(Name_file_t,"%s/FEM_%s_%i.vtk",OutputDir,Name_File,TimeStep_i);
  
  Vtk_file = fopen(Name_file_t,"w");

  /* Header */
  fprintf(Vtk_file,"# vtk DataFile Version 3.0 \n");
  fprintf(Vtk_file,"vtk output \n");
  fprintf(Vtk_file,"ASCII \n");
  fprintf(Vtk_file,"DATASET UNSTRUCTURED_GRID \n");

   /* Coordinates */
  fprintf(Vtk_file,"POINTS %i float \n",ElementMesh.NumNodesMesh);
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    fprintf(Vtk_file,"%f %f %f \n",
	    ElementMesh.Coordinates.nM[i][0],
	    ElementMesh.Coordinates.nM[i][1],
	    ElementMesh.Coordinates.nM[i][2]);
  }

  /* Connectivity */
  fprintf(Vtk_file,"CELLS %i %i \n",
	  ElementMesh.NumElemMesh,
	  ElementMesh.NumElemMesh*(1+ElementMesh.NumNodesElem));
  for(int i = 0 ; i<ElementMesh.NumElemMesh ; i++){
    fprintf(Vtk_file,"%i",ElementMesh.NumNodesElem);
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      fprintf(Vtk_file," %i ",ElementMesh.Connectivity[i][j]);
    }
    fprintf(Vtk_file,"\n");
  }

  /* Type of element */
  fprintf(Vtk_file,"CELL_TYPES %i \n",ElementMesh.NumElemMesh);
  for(int i = 0 ; i<ElementMesh.NumElemMesh ; i++){
    fprintf(Vtk_file,"%i \n",9);
  }

  /* Point data */
  fprintf(Vtk_file,"POINT_DATA %i \n",ElementMesh.NumNodesMesh);
  fprintf(Vtk_file,"VECTORS %s float \n","X_GC");
  for(int i =  0 ; i<ElementMesh.NumNodesMesh ; i++){
    fprintf(Vtk_file,"%f %f %f \n",
	    ElementMesh.Coordinates.nM[i][0],
	    ElementMesh.Coordinates.nM[i][1],
	    ElementMesh.Coordinates.nM[i][2]);
  }

 
  fprintf(Vtk_file,"SCALARS Nod_Mass float \n");
  fprintf(Vtk_file,"LOOKUP_TABLE default \n");
  for(int i =  0 ; i<ElementMesh.NumNodesMesh ; i++){
    fprintf(Vtk_file,"%f \n",List_Nod_Fields.nM[0][i]);
  }
  
  /* Cell data */  
  fprintf(Vtk_file,"CELL_DATA %i \n",ElementMesh.NumElemMesh);

  fprintf(Vtk_file,"SCALARS Active_Elem int 1 \n");
  fprintf(Vtk_file,"LOOKUP_TABLE table_Elem \n");
  for(int i = 0 ; i<ElementMesh.NumElemMesh ; i++){
    fprintf(Vtk_file,"%i\n",
	    ElementMesh.ActiveElem[i]);
  } 

}
