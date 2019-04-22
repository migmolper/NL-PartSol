#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/ShapeFunctions.h"
#include "ToolsLib/Constitutive.h"

/*********************************************************************/

GaussPoint * AllocateGaussPoints(int NumberGP)
/* 
   Gives to you the total size of the material point mesh :
   Inputs : Total number of material point in the mesh (Gauss points)
*/
{
  GaussPoint * GP_Mesh = (GaussPoint *)Allocate_Array(NumberGP,
						       sizeof(GaussPoint));

  return GP_Mesh;
}


/*********************************************************************/

void Initialize_GP(GaussPoint * GP,
		   int GP_id,
		   Matrix RefCoords,
		   double PoissonRatio,
		   double YoungModulus)
/* In this function we Initialize and allocate all the variables of the Gauss Points,
   if they are necessary */
{

  /* Initialize id of the GP */
  GP->id = GP_id;

  /* Id of the element where it is */
  GP->Element_id = 0; 

  /* Initialize kind of material */
  strcpy(GP->Material,"Elastic");

  /* Initialize position field (Vectorial) in global coordiantes
   and in element coordinates : 
   It is not necessary to allocate memory...think about it ;) */
  
  /* GP.x_GC.nV = ... Initialize */
  GP->x_EC = RefCoords;

  /* Allocate and Initialize the Velocity and acceleration field (Vectorial) */
  GP->v = MatAlloc(2,1);
  GP->a = MatAlloc(2,1);
  /*  GP.v.nV = ... Initialize */
  /*  GP.v.nV = ... Initialize */

  /* Allocate and Initialize the Stress and Strain fields (Tensorial) */
  GP->Stress = MatAlloc(3,1);
  GP->Strain = MatAlloc(3,1);
  /* GP.Strain.nM = ...  Initialize */
  /* GP.Stress.nM = ... Initialize */

  /* Allocate and Initialize the constitutive response */
  GP->D = LinearElastic(PoissonRatio,YoungModulus);

}

/*********************************************************************/

Matrix Get_Lagrangian_CG_Tensor(GaussPoint * GP)
/* 
   The Right Cauchy Green Deformation Tensor :
   C = F^{T} \cdot F 

   Inputs : Gauss point
   Output : C
*/
{
  /* Check if we have a null matrix */
  if (GP->F.nM == NULL){
    puts("Error in Get_Lagrangian_CG_Tensor : GP->F tensor = NULL");
    exit(0);
  }

  /* Output, we dont need to allocate it because the memory reserve is done 
     in the function Scalar_prod() */
  Matrix C;
  /* Make a temporal copy of F because Scalar_prod() destroy the input */
  Matrix F = CopyMat(GP->F);
  /* Get C */
  C = Scalar_prod(Transpose_Mat(F),F);

  return C;
  
}

/*********************************************************************/

Matrix Get_Eulerian_CG_Tensor(GaussPoint * GP) 
/*
 The Left Cauchy Green Deformation Tensor :
 B = F \cdot F^{T} 

 Inputs : Gauss point 
 Output : B
*/
{  
  /* Check if we have a null matrix */
  if (GP->F.nM == NULL){
    puts("Error in Get_Eulerian_CG_Tensor : GP->F tensor = NULL");
    exit(0);
  }

  /* Output, we dont need to allocate it because the memory reserve is done 
   in the function Scalar_prod() */
  Matrix B;
  /* Make a temporal copy of F because Scalar_prod() destroy the input */
  Matrix F = CopyMat(GP->F);
  
  /* Get B */
  B = Scalar_prod(F,Transpose_Mat(F));

  return B;
  
}
