
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "ShapeFunctions.h"
#include "../Constitutive/Constitutive.h"

/*********************************************************************/

/* Function for defining the element properties */
Element Initialize_Element(GaussPoint * GP_e,
			   char * TypeElement,
			   Matrix GlobalNodsCoords,
			   int * GlobalNodsId){

  
  Element Elem;

  /* Linear quadrilateral of four nodes */
  if (strcmp(TypeElement,"Q4")==0){

    /* Identification number of the element */
    Elem.id = 0;

    /* Number of nodes in this element */
    Elem.NumberNodes = 4;

    /* Number of degree of freedom for each node */
    Elem.NumberDOF = 2;
    
    /* Number of Gauss points inside of the element */
    Elem.NumberGP = 4;

    /* Conectivity mesh :
       Here we only pass by reference the connectivity mesh, the 
       allocation of it occurs ouside of the Initialize_Element() */
    Elem.N_id = GlobalNodsId;
    
    /* Global coordiantes of the nodes :
       x and y coordiantes of each node of the element (Nnodes X Ndof) */
    Elem.X_g = GlobalNodsCoords;
    /*
      ^
      |_____ Future updates : remove this field and use only the connectivity mesh
     */
    
    /* Shape functions and its derivatives:
       Here we only pass by reference the function, the output Matrix 
       is allocated insede of N_ref() and dNdX_ref() */
    Elem.N_ref = Q4;
    Elem.dNdX_ref = dQ4;

    /* List of GP inside of the element :
       This list is generated outside of the element and passed by reference */
    Elem.GP_e = GP_e;
    /* Auxiliar matrix to do the operations in K matrix */
    Elem.B = MatAlloc(4,2);
  }

  return Elem;
  
}

/*********************************************************************/

void Get_RefDeformation_Gradient(GaussPoint * GP,
				 Element * Elem_GP)
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
  Matrix X_g = CopyMat(Elem_GP->X_g);
  
  GP->F_ref = Scalar_prod(Elem_GP->dNdX_ref(GP->x_EC), /* grad(N_{\alpha}) */
			 X_g);  /* x_{\alpha} */
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

Matrix Get_dNdx_matrix(Element * Elem_GP,
		       GaussPoint * GP)
/* 
   Get the derivative matrix of the shape functions in global coordiantes, it is
   the so called B matrix in the classical formulation of the finite element method
   Inputs: Element, GP where to evaluate it
   Outputs : Matrix dNdX
*/
{
  /* Variables */
  Matrix dNdX;
  Matrix dNdX_ref;
  Matrix F_ref;
  Matrix F_ref_m1T;
  Matrix dN_di;
  Matrix dN_di_ref;
  Matrix aux;

  /* Get the values of the shape functions derivatives in each node */
  dNdX_ref = Elem_GP->dNdX_ref(GP->x_EC);
  
  
  switch(Elem_GP->NumberDOF){
    
  case 1:
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
    
  case 2:
    
    /* Allocate the output */
    dNdX = MatAlloc(3,2*Elem_GP->NumberNodes);
    
    /* Fill the matrix */
    for(int i = 0 ; i<Elem_GP->NumberNodes ; i++){

      /* Fill an array with the nodal derivatives in the reference element */
      /* Allocate the auxiliar array */
      dN_di_ref = MatAlloc(2,1);
      dN_di_ref.nV[0] = dNdX_ref.nM[0][i];
      dN_di_ref.nV[1] = dNdX_ref.nM[1][i];

      /* Get the inverse of the reference deformation gradient and transpose it */
      F_ref = CopyMat(GP->F_ref);
      F_ref_m1T = Transpose_Mat(F_ref);
	    
      /* Obtain the nodal derivarives in the real element */
      dN_di = Scalar_prod(F_ref_m1T,dN_di_ref);
      
      /* Fill the array with the nodal partial derivation of the reference element */     
      dNdX.nM[0][2*i] = dN_di.nV[0];
      dNdX.nM[1][2*i] = 0; 
      dNdX.nM[2][2*i] = dN_di.nV[1]; 

      dNdX.nM[0][2*i + 1] = 0; 
      dNdX.nM[1][2*i + 1] = dN_di.nV[1]; 
      dNdX.nM[2][2*i + 1] = dN_di.nV[0];
      
    }

    break;
    
  case 3:
    puts("Error in Get_dNdi_matrix() : 3D cases not implemented yet");
    exit(0);
    
  default :
    puts("Error in Get_dNdi_matrix() : Wrong case select");
    exit(0);
  }
  
  /* Free memory */
  free(dNdX_ref.nM);
  free(dN_di.nV);  
  
  return dNdX;
}

/*********************************************************************/


Matrix Get_Stiffness_Matrix(Element * Elem) /* Element to analyze */
/* 
   Calcule the element-stiffness matrix K :
   Inputs : Constitutive matrix
   Outputs : Element
*/
{

  Matrix dNdX;
  Matrix dNdX_T;
  Matrix K;
  Matrix K_GP;
  Matrix D;
  Matrix J;

  /* Auxiliar pointer */
  GaussPoint * GP_ei;
  
  /* Allocate and initialize K to zero */
  K = MatAllocZ(Elem->NumberNodes * Elem->NumberDOF,
		 Elem->NumberNodes * Elem->NumberDOF);
    
  /* Iterate over the Gauss Points */
  for(int i_GP = 0 ; i_GP<Elem->NumberGP ; i_GP++){

    /* Get the derivative Matrix and it transposse evaluated in the gauss point */

    /* Pointer to the GP in the element */
    GP_ei = &Elem->GP_e[i_GP];
    
    dNdX = Get_dNdx_matrix(Elem,GP_ei);
    dNdX_T= Transpose_Mat(Get_dNdx_matrix(Elem,GP_ei));

    /* Get the jacobian of the transformation and store it as a scalar */
    J.n = Get_Determinant(GP_ei->F_ref);
    J.nM = NULL;
    J.nV = NULL;
    J.N_cols = 1;
    J.N_rows = 1;

    /* Multiply the constitutive matrix by the derivative matrix */
    D = CopyMat(GP_ei->D);    
    K_GP = Scalar_prod(D,dNdX);

    /* Multiply the result by the derivative matrix transpose */
    K_GP = Scalar_prod(dNdX_T,K_GP);
    
    /* Multiply by the jacobian*/
    K_GP = Scalar_prod(K_GP,J);

    /* Acumulate for the integral rule */
    K = Add_Mat(K,K_GP);
    
  }
  
  return K;
}

