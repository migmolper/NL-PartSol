#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../InOutFun/InOutFun.h"
#include "MeshTools.h"

Mesh InitializeMesh(char * GDF){

  Mesh Back_Mesh;

  puts("*************************************************");
  puts(" Generate the MPM mesh");
  puts(" \t Defining background FEM mesh ...");
  Back_Mesh = ReadGidMesh(FEM_MeshFileName);
  puts(" \t DONE !!!");
  puts(" \t Searching neighbours elements for each node ...");
  Back_Mesh.NodeNeighbour = GetNodalConnectivity(Back_Mesh);
  puts(" \t DONE !!!");
  puts(" \t Reading FEM boundary conditions ...");
  Back_Mesh.Bounds = Set_FEM_BCC(GDF, Back_Mesh);
  puts(" \t DONE !!! ");

  return Back_Mesh;
}
