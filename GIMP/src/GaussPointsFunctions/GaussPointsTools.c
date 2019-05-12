#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "../InOutFun/InOutFun.h"
#include "../ElementsFunctions/ShapeFunctions.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../Constitutive/Constitutive.h"
#include "GaussPointsTools.h"

/*********************************************************************/

GaussPoint Initialize_GP_Mesh(Matrix InputFields,
			      Matrix D,
			      Element Elem)
/*
  
*/
{
  /* Material point mesh (Gauss-Points) */
  GaussPoint GP_Mesh;
  int Init_Num_GP_Elem;
  int Size_GP_Mesh;
  int NumFields;
  char * Field[MAXW] = {NULL};
  /* Initialize parser to read files */
  ParserDictionary Dict = InitParserDictionary();
  char * delims = Dict.sep[4];

  /* Screen message */
  printf("************************************************* \n");
  printf("Begin of initialize the Gauss-Points mesh !!! \n");

  /* do this with a swich...*/
  if (strcmp(Elem.TypeElem,"Linear")==0){
    if(Elem.NumNodesElem == 2){ /* 1D Linear mesh */
      Init_Num_GP_Elem = 1;
      Size_GP_Mesh = Init_Num_GP_Elem*(Elem.NumElemMesh);
      GP_Mesh.Phi.x_EC = MatAlloc(1,Size_GP_Mesh);
    }
  }
  
  /* Initialize id of the GP */
  GP_Mesh.NumGP = Size_GP_Mesh;
  /* Initialize all the fields :
     - Mass of the material point  
       - Position field (Vectorial) in global coordiantes and in element coordinates : 
       - Displacement, Velocity and acceleration field (Vectorial)  
       - Stress and Strain fields (Tensorial) 
       Note : It is not necessary to allocate memory...think about it ;) 
  */

  NumFields = parse (Field, InputFields.Info, ";\n");
    
  for(int i = 0; i<NumFields ; i++){
    if(strcmp(Field[i],"X_GP")==0)
      GP_Mesh.Phi.x_GC.nV = InputFields.nM[i];

    if(strcmp(Field[i],"V_X")==0)
      GP_Mesh.Phi.vel.nV = InputFields.nM[i];

    if(strcmp(Field[i],"SIGMA_X")==0)
      GP_Mesh.Phi.Stress.nV = InputFields.nM[i];

    if(strcmp(Field[i],"MASS")==0)
      GP_Mesh.Phi.mass.nV = InputFields.nM[i];
  }
  
  /* Id of the element set to a NaN value */ 
  GP_Mesh.Element_id = (int *)Allocate_Array(Size_GP_Mesh,sizeof(int));
  for(int i = 0 ; i<Size_GP_Mesh ; i++){
    GP_Mesh.Element_id[i] = -999;
  }
  
  /* Allocate and Initialize the constitutive response */
  GP_Mesh.D = D;

  /* Localize all the Gauss-Points */
  UpdateElementLocationGP(GP_Mesh,Elem);

  /* Initial conditions */

  for(int i = 0 ; i<GP_Mesh.NumGP ;i++){
    printf("* Set the GP %i initial conditions : \n",i+1);
    printf("\t -> Element : %i \n",
	   GP_Mesh.Element_id[i]+1);
    printf("\t -> X Glob coordinate : %f \n",
	   GP_Mesh.Phi.x_GC.nV[i]);
    printf("\t -> X Ref coordiante : %f \n",
	   GP_Mesh.Phi.x_EC.nV[i]);
    printf("\t -> Velocity : %f \n",
	   GP_Mesh.Phi.vel.nV[i]);
    printf("\t -> Stress : %f \n",
	   GP_Mesh.Phi.Stress.nV[i]);
    printf("\t -> Mass : %f \n",
	   GP_Mesh.Phi.mass.nV[i]);
  }

  /* Final messages */
  printf("End of initialize the Gauss-Points mesh !!! \n");
  
  return GP_Mesh;

}

/*********************************************************************/

void UpdateElementLocationGP(GaussPoint GP_Mesh,Element ElementMesh){

  double GP_x;
  int Elem_1D_0,Elem_1D_1;
  
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){
    if(GP_Mesh.Element_id[i] == -999){
      GP_x = GP_Mesh.Phi.x_GC.nV[i];
      for(int j = 0 ; j<ElementMesh.NumElemMesh ; j++){
	Elem_1D_0 = ElementMesh.Connectivity[j][0];
	Elem_1D_1 = ElementMesh.Connectivity[j][1];
	if( (GP_x >= ElementMesh.Coordinates[Elem_1D_0-1][0]) &&
	    (GP_x <= ElementMesh.Coordinates[Elem_1D_1-1][0]) ){
	  GP_Mesh.Element_id[i] = j;
	  GP_Mesh.Phi.x_EC.nV[i] = 0.5;
	  /* Activate element */
	  ElementMesh.ActiveElem[i] = 1;
	}
      }   
    }
  }
  
}

/*********************************************************************/

/* double GetJacobian(){ */
/* } */

/*********************************************************************/


/* Matrix Get_Lagrangian_CG_Tensor(GaussPoint * GP) */
/* /\*  */
/*    The Right Cauchy Green Deformation Tensor : */
/*    C = F^{T} \cdot F  */

/*    Inputs : Gauss point */
/*    Output : C */
/* *\/ */
/* { */
/*   /\* Check if we have a null matrix *\/ */
/*   if (GP->F.nM == NULL){ */
/*     puts("Error in Get_Lagrangian_CG_Tensor : GP->F tensor = NULL"); */
/*     exit(0); */
/*   } */

/*   /\* Output, we dont need to allocate it because the memory reserve is done  */
/*      in the function Scalar_prod() *\/ */
/*   Matrix C; */
/*   /\* Make a temporal copy of F because Scalar_prod() destroy the input *\/ */
/*   Matrix F = CopyMat(GP->F); */
/*   /\* Get C *\/ */
/*   C = Scalar_prod(Transpose_Mat(F),F); */

/*   return C; */
  
/* } */

/*********************************************************************/

/* Matrix Get_Eulerian_CG_Tensor(GaussPoint * GP)  */
/* /\* */
/*  The Left Cauchy Green Deformation Tensor : */
/*  B = F \cdot F^{T}  */

/*  Inputs : Gauss point  */
/*  Output : B */
/* *\/ */
/* {   */
/*   /\* Check if we have a null matrix *\/ */
/*   if (GP->F.nM == NULL){ */
/*     puts("Error in Get_Eulerian_CG_Tensor : GP->F tensor = NULL"); */
/*     exit(0); */
/*   } */

/*   /\* Output, we dont need to allocate it because the memory reserve is done  */
/*    in the function Scalar_prod() *\/ */
/*   Matrix B; */
/*   /\* Make a temporal copy of F because Scalar_prod() destroy the input *\/ */
/*   Matrix F = CopyMat(GP->F); */
  
/*   /\* Get B *\/ */
/*   B = Scalar_prod(F,Transpose_Mat(F)); */

/*   return B; */
  
/* } */
