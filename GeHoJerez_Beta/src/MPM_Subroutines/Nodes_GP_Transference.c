#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../InOutFun/InOutFun.h"

Matrix GetNodalValuesFromGP(GaussPoint MPM_Mesh, Mesh FEM_Mesh, char LisOfFields[MAXC])
{

  /* Variable declaration */
  Matrix Nod_Fields; /* Output list of fields */
  Matrix GP_Fields; /* Auxiliar structure to do the selection */
  int Elem_GP; /* Index of the element */
  int * Elem_Nods;
  Matrix X_EC_GP;
  Matrix N_Ref_GP, dNdX_Ref_GP;
  Matrix F_Ref_GP, F_Ref_GP_T, F_Ref_GP_Tm1;
  Matrix dNdx_XG;
  char * FieldsList[MAXW] = {NULL};
  int NumberFields;
  Matrix Elem_Coords = MatAllocZ(FEM_Mesh.NumNodesElem,FEM_Mesh.Dimension);
  int CalcGradient = 0; /* Boolean variable, set to zero by defaul */


  /* 0º Get those fields to transfeer */
  NumberFields = parse (FieldsList,LisOfFields,";\n");
  if(NumberFields<=0){
    puts("Error in GaussPointToMesh() : Wrong number of fields !!!");
    exit(0);
  }
  
  /* 1º Allocate the output list of fields */
  Nod_Fields.N_rows = NumberFields;
  Nod_Fields.N_cols = FEM_Mesh.NumNodesMesh;
  Nod_Fields.n = NAN;
  Nod_Fields.nV = NULL;
  Nod_Fields.nM = (double **)Allocate_MatrixZ(NumberFields,
					      FEM_Mesh.NumNodesMesh,
					      sizeof(double));

  
  /* 2º Iterate over the GP to get the nodal values */
  for(int i = 0 ; i<MPM_Mesh.NumGP ; i++){

    /* 3º Get the coordinates of the GP in the Element */
    X_EC_GP.nV = MPM_Mesh.Phi.x_EC.nM[i];

    /* 4º Get the index of the element */
    Elem_GP = MPM_Mesh.Element_id[i];

    /* 5º Nodes of the element */
    Elem_Nods = FEM_Mesh.Connectivity[Elem_GP];

    /* 6º Coordinates of the element */
    for(int j = 0 ; j<FEM_Mesh.NumNodesElem ; j++){
      for(int k = 0 ; k<FEM_Mesh.Dimension ; k++){
	Elem_Coords.nM[j][k] = FEM_Mesh.Coordinates.nM[Elem_Nods[j]][k];
      }
    }

    /* 7º Evaluate the shape function and it derivarive in the GP */
    N_Ref_GP = FEM_Mesh.N_ref(X_EC_GP);

    /* Do this only when we have to evaluate a gradient, for example to get the
     internal forces */
    if(CalcGradient){
      /* 8º Evaluate the gradient of the shape function in the GP */
      dNdX_Ref_GP = FEM_Mesh.dNdX_ref(X_EC_GP);

      /* 9º Get the reference deformation gradient in the GP */
      F_Ref_GP = Get_RefDeformation_Gradient_Q4(X_EC_GP,Elem_Coords);

      /* 10º Get the deformation gradient (dNdx_XG) : Only in some cases */
      F_Ref_GP_T = Transpose_Mat(F_Ref_GP),
	free(F_Ref_GP.nM);
      F_Ref_GP_Tm1 = Get_Inverse(F_Ref_GP_T),
	free(F_Ref_GP_T.nM);
      dNdx_XG = Scalar_prod(F_Ref_GP_Tm1,dNdX_Ref_GP),
	free(F_Ref_GP_Tm1.nM),
	free(dNdX_Ref_GP.nM);
    }

    /* 11º For each nodal field that we want to calc : */
    for(int j = 0 ; j<NumberFields ; j++){

      /* 12aº Asign GP mass to the nodes */
      if(strcmp(FieldsList[j],"MASS") == 0){
      	for(int k = 0 ; k<FEM_Mesh.NumNodesElem ; k++){
      	  Nod_Fields.nM[j][Elem_Nods[k]] += MPM_Mesh.Phi.mass.nV[i]*N_Ref_GP.nV[k];
      	}
      }

      /* 12bº Asign GP momentum to the nodes */
      if(strcmp(FieldsList[j],"MOMENTUM") == 0){
	for(int k = 0 ; j<FEM_Mesh.NumNodesElem ; k++){
	  for(int l = 0 ; l<NumberDimensions ; l++){
	    Nod_Fields.nM[j][Elem_Nods[k]] +=
	      MPM_Mesh.Phi.mass.nV[i]*
	      MPM_Mesh.Phi.vel.nM[i][l]*
	      N_Ref_GP.nV[j];
	  }
	}
      }

      /* 12cº Asign GP internal forces to the nodes */
      if(strcmp(FieldsList[j],"F_INT") == 0){
	
	puts("Functionality not implemented yet");
	exit(0);
      }

      /* 12dº Asign GP body forces to the nodes */
      if(strcmp(FieldsList[j],"F_B") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }

      /* 12eº Asign GP local forces to the nodes */
      if(strcmp(FieldsList[j],"F_T") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }

      /* 12fº Asign GP total forces to the nodes */
      if(strcmp(FieldsList[j],"F_TOT") == 0){
	puts("Functionality not implemented yet");
	exit(0);
      }
    
    }
    
  }

  return Nod_Fields;
  
}
