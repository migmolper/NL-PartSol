#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/Utils.h"
#include "ShapeFunctions.h"
#include "../Constitutive/Constitutive.h"
#include "ElementTools.h"

/*********************************************************************/

Matrix Get_RefDeformation_Gradient(Matrix X_g,
				   Matrix dNdX_Ref_GP)
/*
  Get the deformation gradient of the reference element:

  F = grad(N_{\alpha}) \otimes x_{\alpha}
  
  Inputs :
  - X_g -> This are the coordinates of the nodes
  - dNdX_Ref_GP -> Derivative gradient evaluated in the GP

  Output :
  - F_Ref -> Deformation gradient of the reference element
  evaluated in the GP
*/
{
  /* Output declaration */
  Matrix F_Ref;
  
  F_Ref = Scalar_prod(dNdX_Ref_GP, /* grad(N_{\alpha}) */
		      X_g);  /* x_{\alpha} */
    

  return F_Ref;
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

Matrix Get_dNdx_matrix(Matrix X_g,
		       Matrix dNdX_Ref_GP,
		       int NumNodesElem,
		       int Ndim)
/*
   Get the derivative matrix of the shape functions in global coordiantes, it is
   the so called B matrix in the classical formulation of the finite element method
   Inputs: 
   - dNdX_ref_GP : Values of the shape functions derivatives in each node
   evaluate in the GP 

   Outputs : Matrix dNdX_GP
*/
{
  /* Decalaration of the output matrix */
  Matrix dNdX_GP;
  /* Inverse and transpose of the reference deformation gradient */
  Matrix F_Ref;
  Matrix F_Ref_m1;
  Matrix F_Ref_m1T;
  /* Nodal values of the reference shape functions derivatives in each node
     evaluate in the GP */
  Matrix dN_di_Ref_GP;
  /* Nodal derivative matrix evaluated in the GP */
  Matrix dN_di_GP;
  
  switch(Ndim){
    
  case 1:
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
    
  case 2:
    
    /* Allocate the output */
    dNdX_GP = MatAlloc(3,2*NumNodesElem);

    /* Allocate the auxiliar array */
    dN_di_Ref_GP = MatAlloc(2,1);

    /* Get the inverse of the reference deformation gradient and transpose it */
    F_Ref = Get_RefDeformation_Gradient(X_g,dNdX_Ref_GP);
    F_Ref_m1 = Get_Inverse(F_Ref);
    free(F_Ref.nM);
    F_Ref_m1T = Transpose_Mat(F_Ref_m1);
    free(F_Ref_m1.nM);
    
    /* Fill the matrix */
    for(int i = 0 ; i<NumNodesElem ; i++){

      /* Fill an array with the nodal derivatives in the reference element */
      dN_di_Ref_GP.nV[0] = dNdX_Ref_GP.nM[0][i];
      dN_di_Ref_GP.nV[1] = dNdX_Ref_GP.nM[1][i];
	    
      /* Obtain the nodal derivarives in the real element */
      dN_di_GP = Scalar_prod(F_Ref_m1T,dN_di_Ref_GP);
      
      /* Fill the array with the nodal partial derivation of the reference element */
      dNdX_GP.nM[0][2*i] = dN_di_GP.nV[0];
      dNdX_GP.nM[1][2*i] = 0;
      dNdX_GP.nM[2][2*i] = dN_di_GP.nV[1];

      dNdX_GP.nM[0][2*i + 1] = 0;
      dNdX_GP.nM[1][2*i + 1] = dN_di_GP.nV[1];
      dNdX_GP.nM[2][2*i + 1] = dN_di_GP.nV[0];
      
    }

    /* Free memory */
    free(dN_di_GP.nV);
    free(dN_di_Ref_GP.nV);
    free(F_Ref_m1T.nM);

    break;
    
  case 3:
    puts("Error in Get_dNdi_matrix() : 3D cases not implemented yet");
    exit(0);
    
  default :
    puts("Error in Get_dNdi_matrix() : Wrong case select");
    exit(0);
  }
  
  return dNdX_GP;
}

/*********************************************************************/


/* Matrix Get_Stiffness_Matrix(Element * Elem) /\* Element to analyze *\/ */
/* /\*  */
/*    Calcule the element-stiffness matrix K : */
/*    Inputs : Constitutive matrix */
/*    Outputs : Element */
/* *\/ */
/* { */

/*   Matrix dNdX; */
/*   Matrix dNdX_T; */
/*   Matrix K; */
/*   Matrix K_GP; */
/*   Matrix D; */
/*   Matrix J; */

/*   /\* Auxiliar pointer *\/ */
/*   GaussPoint * GP_ei; */
  
/*   /\* Allocate and initialize K to zero *\/ */
/*   K = MatAllocZ(Elem->NumberNodes * Elem->NumberDOF, */
/* 		 Elem->NumberNodes * Elem->NumberDOF); */
    
/*   /\* Iterate over the Gauss Points *\/ */
/*   for(int i_GP = 0 ; i_GP<Elem->NumberGP ; i_GP++){ */

/*     /\* Get the derivative Matrix and it transposse evaluated in the gauss point *\/ */

/*     /\* Pointer to the GP in the element *\/ */
/*     GP_ei = &Elem->GP_e[i_GP]; */
    
/*     dNdX = Get_dNdx_matrix(Elem,GP_ei); */
/*     dNdX_T= Transpose_Mat(Get_dNdx_matrix(Elem,GP_ei)); */

/*     /\* Get the jacobian of the transformation and store it as a scalar *\/ */
/*     J.n = Get_Determinant(GP_ei->F_ref); */
/*     J.nM = NULL; */
/*     J.nV = NULL; */
/*     J.N_cols = 1; */
/*     J.N_rows = 1; */

/*     /\* Multiply the constitutive matrix by the derivative matrix *\/ */
/*     D = CopyMat(GP_ei->D);     */
/*     K_GP = Scalar_prod(D,dNdX); */

/*     /\* Multiply the result by the derivative matrix transpose *\/ */
/*     K_GP = Scalar_prod(dNdX_T,K_GP); */
    
/*     /\* Multiply by the jacobian*\/ */
/*     K_GP = Scalar_prod(K_GP,J); */

/*     /\* Acumulate for the integral rule *\/ */
/*     K = Add_Mat(K,K_GP); */
    
/*   } */
  
/*   return K; */
/* } */

/*********************************************************************/

Matrix Get_Geom_Mass_Matrix(GaussPoint GP_Mesh,
			    Element ElementMesh) 
/*
   Calcule the mass matrix M :
   Inputs : Element
   Outputs : Mass matrix (M)
*/
{
  /* Mass matrix */
  Matrix M;
  /* Index of the element where it is the GP */
  int i_GP_Elem;
  
  /* Coordenates of the element nodes */
  Matrix X_ElemNod;
  /* Pointer to the Connectivity of the element */
  int * Id_ElemNod;
  /* Derivative Matrix */
  Matrix dNdX_Ref_GP;

  /* Jacobian */
  Matrix J_GP;
  
  /* Allocate and initialize M to zero */
  M = MatAllocZ(ElementMesh.NumNodesMesh,
		ElementMesh.NumNodesMesh);

  /* Allocate the mesh with the nodal coordinates */
  X_ElemNod = MatAlloc(ElementMesh.NumNodesElem,
		       ElementMesh.Dimension);

  /* Iterate over the Gauss Points */
  for(int i_GP = 0 ; i_GP<GP_Mesh.NumGP ; i_GP++){
    
    i_GP_Elem = GP_Mesh.Element_id[i_GP];

    /* Calcule the Jacobian of the element evaluated in the GP : */
    /* Get the element connectivity and take in to account that 
       in the C programming languaje the index starts in 0 */
    Id_ElemNod = ElementMesh.Connectivity[i_GP_Elem];
    /* Get the coordinates of the nodes */
    X_ElemNod.nV[0] = ElementMesh.Coordinates[Id_ElemNod[0] - 1][0];
    X_ElemNod.nV[1] = ElementMesh.Coordinates[Id_ElemNod[1] - 1][0];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(X_ElemNod);
    J_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);

    M.nM[Id_ElemNod[0]-1][Id_ElemNod[0]-1] += (double)1/3 * J_GP.n;
    M.nM[Id_ElemNod[0]-1][Id_ElemNod[1]-1] += (double)1/6 * J_GP.n;
    M.nM[Id_ElemNod[1]-1][Id_ElemNod[0]-1] += (double)1/6 * J_GP.n;
    M.nM[Id_ElemNod[1]-1][Id_ElemNod[1]-1] += (double)1/3 * J_GP.n;
    
  }

  return M;


}
