#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/ShapeFunctions.h"

/*********************************************************************/

void Initialize_GP(GaussPoint GP,Matrix RefCoords)
/* In this function we Initialize and allocate all the variables of the Gauss Points,
   if they are necessary */
{

  /* Initialize_Element id of the GP */
  GP.id = 0;

  /* Id of the element where it is */
  GP.Element_id = 0; 

  /* Initialize kind of material */
  strcpy(GP.Material,"Elastic");

  /* Initialize and allocate position field (Vectorial) in global coordiantes
   and in element coordinates */
  GP.x_GC = MatAlloc(1,2);
  GP.x_EC = MatAlloc(1,2);
  /* GP.x_GC.nV = ... Initialize */
  GP.x_EC.nV = RefCoords.nV;

  /* Allocate and Initialize the Velocity and acceleration field (Vectorial) */
  GP.v = MatAlloc(1,2);
  GP.a = MatAlloc(1,2);
  /*  GP.v.nV = ... Initialize */
  /*  GP.v.nV = ... Initialize */

  /* Allocate and Initialize the Stress and Strain fields (Tensorial) */
  GP.Stress = MatAlloc(2,2);
  GP.Strain = MatAlloc(2,2);
  /* GP.Strain.nM = ...  Initialize */
  /* GP.Stress.nM = ... Initialize */
}

/*********************************************************************/


/* Function for defining the element properties */
void Initialize_Element(Element Elem,
			GaussPoint * GP_e,
			char * TypeElement,
			Matrix GlobalNodsCoords,
			int * GlobalNodsId){
  
  if (strcmp(TypeElement,"Q4")==0){

    /* Identification number of the element */
    Elem.id = 0;

    /* Number of nodes in this element */
    Elem.NumberNodes = 4;

    /* Number of degree of freedom for each node */
    Elem.NumberDOF = 2;

    /* Global coordiantes of the nodes */
    Elem.X_g = GlobalNodsCoords;

    /* Conectivity mesh */
    Elem.N_id = GlobalNodsId;

    /* Shape functions and its derivatives */
    Elem.N_ref = Q4;
    Elem.dNdX_ref = dQ4;

    /* List of GP inside of the element */
    Elem.GP_e = GP_e;

    /* Auxiliar matrix to do the operations in K matrix */
    Elem.B = MatAlloc(4,2);
  }
  
}



/*********************************************************************/

void Get_RefDeformation_Gradient(GaussPoint GP,
				 Element Elem_GP)
/*
  Get the deformation gradient of the reference element:

  F = grad(N_{\alpha}) \otimes x_{\alpha}
  
  Inputs : 
  - GP -> GP structure
  - Elem -> Structure with the properties of the element

  Output :
  - F_ref -> Deformation gradient of the reference element
  evaluated in the GP
  - B_ref -> Ma
*/
{
  GP.F_ref = Scalar_prod(Elem_GP.dNdX_ref(GP.x_EC), /* grad(N_{\alpha}) */
			 Elem_GP.X_g);  /* x_{\alpha} */
}


/*********************************************************************/

/* void Get_Deformation_Gradient(Element * Elem, /\* Gauss point *\/ */
/* 			      GaussPoint * GP) /\* Element *\/ */
/* { */
  
/*   /\* Deformation reference gradient F_ref *\/ */
/*   Get_RefDeformation_Gradient(GP,Elem); */

/*   /\* Get the inverse of the deformation gradient *\/ */
/*   Tensor F_ref_m1; */
/*   F_ref_m1 = Get_Inverse(GP[0].F_ref); */

/*   /\* Iterate over the nodes in the Element to get B *\/ */
/*   for(int i = 0 ; i<Element[0].NumberNodes ; i++){ */

    
    
/*   }/\* Nodes loop *\/ */

/*   /\* Once we have calculated the deformation gradient */
/*      matrix we dont need anymore the inverse of the */
/*      reference deformation gradient so we free memory *\/ */
/*   free(F_ref_m1.n); */
  
/* } */



/*********************************************************************/

Matrix Get_dNdi_matrix(Matrix dN_di){

  Matrix dNdi_matrix;

  if(dNdi_matrix.nM != NULL){
    puts("Error in Get_dNdi_matrix() : The input must be an array !");
    exit(0);
  }

  switch(dNdi_matrix.N_cols*dNdi_matrix.N_rows){
  case 1:
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
  case 2:
    dNdi_matrix = MatAlloc(3,2);
    dNdi_matrix.nM[0][0] = dN_di[0];
    dNdi_matrix.nM[1][0] = 0; 
    dNdi_matrix.nM[2][0] = dN_di[1]; 
    dNdi_matrix.nM[0][1] = 0; 
    dNdi_matrix.nM[1][1] = dN_di[1]; 
    dNdi_matrix.nM[2][1] = dN_di[0];
    break;
  case 3:
    puts("Error in Get_dNdi_matrix() : 3D cases not implemented yet");
    exit(0);
  default :
    puts("Error in Get_dNdi_matrix() : Wrong case select");
    exit(0);
  }
  
  return dNdi_matrix;

}

/*********************************************************************/

Matrix Get_dNdx_matrix(Element Elem_GP,
		       GaussPoint GP)
/* Get the derivative matrix of the shape functions in global coordiantes, it is
   the so called B matrix in the classical formulation of the finite element method 
   Inputs: Element, GP where to evaluate it
   Outputs : Matrix dNdX
*/
{
  /* Variables */
  Matrix dNdX;
  Matrix dNdX_ref;
  Matrix dNdi_matrix;
  Matrix F_ref_m1T;
  
  /* Allocate the output */
  dNdX = MatAlloc(3,4);

  /* Get the inverse of the reference deformation gradient and transpose it */
  F_ref_m1T = Transpose_Mat(Get_Inverse(GP.F_ref));

  /* Get the values of the shape functions derivatives in each node */  
  dNdX_ref = Elem_GP.dNdX_ref(GP.x_EC);

  for(int i = 0 ; i < Elem_GP.NumberNodes ; i++){  

    /* Ojo, cuidado con las dimensiones del producto !!!!!!!!!!!!!!!*/
    dNdi_matrix = Get_dNdi_matrix(Scalar_prod(F_ref_m1T,dNdX_ref));

  }
  
  free(dNdX_ref.nM);
  free(F_ref_m1T.nM)
  
  return dNdX;  
}


/*********************************************************************/

Matrix Get_Lagrangian_CG_Tensor(GaussPoint GP)
/* 
   The Right Cauchy Green Deformation Tensor :
   C = F^{T} \cdot F 

   Inputs : Gauss point
   Output : C
*/
{
  /* Check if we have a null matrix */
  if (GP[0].F.n == NULL){
    puts("Error in Get_Lagrangian_CG_Tensor : GP.F tensor = NULL");
    exit(0);
  }

  Matrix C;
  
  /* Get C */
  C = Scalar_prod(Transpose_Mat(GP.F), /* F^T */
		  GP.F); /* F */

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
  if (GP[0].F.n == NULL){
    puts("Error in Get_Eulerian_CG_Tensor : GP.F tensor = NULL");
    exit(0);
  }

  Matrix B;
  
  /* Get B */
  B = Scalar_prod(GP.F, /* F */
		  Transpose_Mat(GP.F)); /* F^T */

  return B;
  
}


/*********************************************************************/


/* Matrix Get_Stiffness_Matrix(Matrix D, /\* Constitutive matrix *\/ */
/* 			    Element * Elem) /\* Element to analyze *\/ */
/* /\* Calcule the element stiffness matrix K *\/ */
/* { */

/*   Element Elem_e; */
/*   Matrix dNdX; */
/*   Matrix K; */
/*   double J; */
  
/*   /\* Allocate and initialize K *\/ */
/*   K.N_rows = Elem.NumberNodes * Elem.NumberDOF; */
/*   K.N_cols = Elem.NumberNodes * Elem.NumberDOF; */
/*   K.n = (double **)Allocate_Matrix(K.N_rows,K.N_cols,sizeof(double **)); */
  
/*   /\* Iterate over the Gauss Points *\/ */
/*   for(int i_gp = 0 ; i_gp<Elem.NumberGP ; i_gp++){ */
    
/*     /\* Get the derivatives of the shape function in global coordiantes *\/ */
/*     dNdx = Get_dNdx_matrix(Elem_e,Elem_e.GP_e[i_gp]); */
    
/*     /\* Get the Jacobian of the transformation *\/ */
/*     /\* Check if we have a null matrix *\/ */
/*     J = Get_Determinant(Elem_e.GP_e[i_gp].F_ref); */
    
/*     /\* Iterate nodes *\/ */
/*     K_gp = Mat_Mul(D,dNdX); */
/*     K_gp = Mat_Mul(dNdX,K); */
/*     K_gp *= J; */
    
/*     K = Mat_Sum(K,K_gp); */
    
/*   } */
/*   /\* We dont need the derivatives *\/ */
/*   free(dNdX.n); */
  
/*   return K; */
/* } */

