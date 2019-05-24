#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../ToolsLib/Utils.h"

Mesh RectangularMesh(double X0, double Y0,
		     double X1, double Y1,
		     double Dx, double Dy,
		     char * ElementType){
  Mesh OutMesh;
  int node_i = 0;

  if( strcmp(ElementType,"Linear") == 0 ){
    strcpy(OutMesh.TypeElem,ElementType);
    OutMesh.NumNodesElem = 4;
    OutMesh.Dimension = 2;
    OutMesh.N_ref = Q4;
    OutMesh.dNdX_ref = dQ4;
    OutMesh.NumNodesMesh = (int)(((X1-X0)/Dx + 1) *
				 ((Y1-Y0)/Dy + 1));
    OutMesh.NumElemMesh = (int)(((X1-X0)/Dx) *
				((Y1-Y0)/Dy));
    OutMesh.NumNodesBound = (int)(((X1-X0)/Dx - 1)*2 +
				  ((Y1-Y0)/Dy - 1)*2) + 4;

    /* Allocate matrix with the coordinates of the nodes */
    OutMesh.Coordinates = MatAlloc(OutMesh.NumNodesMesh,2);
    OutMesh.Connectivity = (int **)Allocate_Matrix(OutMesh.NumElemMesh,4,sizeof(int));
    OutMesh.ActiveElem = (int *)Allocate_ArrayZ(OutMesh.NumElemMesh,sizeof(int));

    /* Fill matrix with coordinates */
    for(int i = 0 ; i<(int)((X1-X0)/Dx + 1) ; i++){
      for(int j = 0 ; j<(int)((Y1-Y0)/Dy + 1) ; j++){
	OutMesh.Coordinates.nM[node_i][0] = Dx*i;
	OutMesh.Coordinates.nM[node_i][1] = Dy*j;
	node_i++;
      }
    }
    
    /* Fill conectivity mesh */
        
  }
 
  
  return OutMesh;
}

