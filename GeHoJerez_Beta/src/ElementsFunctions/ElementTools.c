#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "ShapeFunctions.h"
#include "../Constitutive/Constitutive.h"
#include "ElementTools.h"


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

Matrix Get_B_GP(Matrix X_EC_GP,Matrix Element)
/*
   Get the B matrix (Usual in the classical formulation of 
   the finite element method )
   Inputs:
   - Matrix X_NC_GP : Element coordinates
   - Matrix Element : Coordinates of the element 

   Outputs : Matrix B
*/
{

  /* 0º Define variables */
  Matrix B_GP; /* Declaration of the output matrix (NdimVecStrain x Nnodes*Ndim) */
  Matrix dNdX_Ref_GP; /* Derivative of the shape function evaluated in the GP (Ndim x Nnodes) */
  Matrix F_Ref_GP; /* Reference deformation gradient evaluated in the GP (Ndim x Ndim) */
  Matrix F_Ref_GP_m1; /* Inverse of the reference deformation gradient */
  Matrix F_Ref_GP_Tm1; /* Transpose of the inverse reference deformation gradient */
  Matrix dNdx_XG_GP; /* Derivatives of the shape function evaluates in the GP (Ndim x Ndim) */

  /* 1º Select the case to solve */
  switch(NumberDimensions){
    
  case 1:  /* 1D stress tensor */
    
    puts("Error in Get_dNdi_matrix() : 1D cases not implemented yet");
    exit(0);
    
  case 2: /* 2D stress tensor */
    
    /* 2º Allocate the output */
    B_GP = MatAlloc(3,2*Element.N_rows);

    if(strcmp(Element.Info,"Quadrilateral") == 0){      
      /* 3º Evaluate the gradient of the shape function in the GP */
      dNdX_Ref_GP = dQ4(X_EC_GP);
      /* 4º Get the reference deformation gradient in the GP */
      F_Ref_GP = Get_RefDeformation_Gradient_Q4(X_EC_GP,Element);
    }

    /* 5º Get the deformation gradient (dNdx_XG) : Only in some cases */
    /* 6aº Get the inverse of the deformation gradient */
    F_Ref_GP_m1 = Get_Inverse(F_Ref_GP), 
      free(F_Ref_GP.nM);
    /* 6bº Get the transpose of the inverse of the deformation gradient */
    F_Ref_GP_Tm1 = Transpose_Mat(F_Ref_GP_m1), 
      free(F_Ref_GP_m1.nM);
    /* 6cº Get the gradient of the shape functions in global coordinates */
    dNdx_XG_GP = Scalar_prod(F_Ref_GP_Tm1,dNdX_Ref_GP), 
      free(F_Ref_GP_Tm1.nM),
      free(dNdX_Ref_GP.nM);
    
    /* 7º Fill the array with the nodal partial derivation of the reference element */    
    for(int i = 0 ; i<Element.N_rows ; i++){
      B_GP.nM[0][2*i] = dNdx_XG_GP.nM[0][i];
      B_GP.nM[1][2*i] = 0;
      B_GP.nM[2][2*i] = dNdx_XG_GP.nM[1][i];
      
      B_GP.nM[0][2*i + 1] = 0;
      B_GP.nM[1][2*i + 1] = dNdx_XG_GP.nM[1][i];
      B_GP.nM[2][2*i + 1] = dNdx_XG_GP.nM[0][i];      
    }

    /* 8º Free memory */
    free(dNdx_XG_GP.nM);

    break;
    
  case 3: /* 3D stress tensor */
    puts("Error in Get_dNdi_matrix() : 3D cases not implemented yet");
    exit(0);
    
  default :
    puts("Error in Get_dNdi_matrix() : Wrong case select");
    exit(0);
  }
  
  return B_GP;
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

/* Matrix Get_Geom_Mass_Matrix(GaussPoint GP_Mesh, */
/* 			    Mesh ElementMesh)  */
/* /\* */
/*    Calcule the mass matrix M : */
/*    Inputs : Element */
/*    Outputs : Mass matrix (M) */
/* *\/ */
/* { */
/*   /\* Mass matrix *\/ */
/*   Matrix M; */
/*   /\* Index of the element where it is the GP *\/ */
/*   int i_GP_Elem; */
  
/*   /\* Coordenates of the element nodes *\/ */
/*   Matrix X_ElemNod; */
/*   /\* Pointer to the Connectivity of the element *\/ */
/*   int * Id_ElemNod; */
/*   /\* Derivative Matrix *\/ */
/*   Matrix dNdX_Ref_GP; */

/*   /\* Jacobian *\/ */
/*   Matrix J_GP; */
  
/*   /\* Allocate and initialize M to zero *\/ */
/*   M = MatAllocZ(ElementMesh.NumNodesMesh, */
/* 		ElementMesh.NumNodesMesh); */

/*   /\* Allocate the mesh with the nodal coordinates *\/ */
/*   X_ElemNod = MatAlloc(ElementMesh.NumNodesElem, */
/* 		       ElementMesh.Dimension); */

/*   /\* Iterate over the Gauss Points *\/ */
/*   for(int i_GP = 0 ; i_GP<GP_Mesh.NumGP ; i_GP++){ */
    
/*     i_GP_Elem = GP_Mesh.Element_id[i_GP]; */

/*     /\* Calcule the Jacobian of the element evaluated in the GP : *\/ */
/*     /\* Get the element connectivity and take in to account that  */
/*        in the C programming languaje the index starts in 0 *\/ */
/*     Id_ElemNod = ElementMesh.Connectivity[i_GP_Elem]; */
/*     /\* Get the coordinates of the nodes *\/ */
/*     X_ElemNod.nV[0] = ElementMesh.Coordinates.nM[Id_ElemNod[0] - 1][0]; */
/*     X_ElemNod.nV[1] = ElementMesh.Coordinates.nM[Id_ElemNod[1] - 1][0]; */
/*     dNdX_Ref_GP = ElementMesh.dNdX_ref(X_ElemNod); */
/*     J_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP); */

/*     M.nM[Id_ElemNod[0]-1][Id_ElemNod[0]-1] += (double)1/3 * J_GP.n; */
/*     M.nM[Id_ElemNod[0]-1][Id_ElemNod[1]-1] += (double)1/6 * J_GP.n; */
/*     M.nM[Id_ElemNod[1]-1][Id_ElemNod[0]-1] += (double)1/6 * J_GP.n; */
/*     M.nM[Id_ElemNod[1]-1][Id_ElemNod[1]-1] += (double)1/3 * J_GP.n; */
    
/*   } */

/*   return M; */


/* } */


/*********************************************************************/


Matrix GetNaturalCoordinates(Matrix X_EC_GP,Matrix X_GC_GP,Matrix Element_GC_Nod)
/* 
   The function return the natural coordinates of a point inside of the element, 
   Inputs :
   - Coordinates of the element nodes
   - Initial coordinate of the point to start the search 
   - Derivative function of the element

   Depending of the kind of element, we employ differents types of shape functions
*/
{
  
  X_EC_GP = Newton_Rapson(Get_GlobalCoordinates_Q4,Element_GC_Nod,
			  Get_RefDeformation_Gradient_Q4,Element_GC_Nod,
			  X_GC_GP,X_EC_GP);

  return X_EC_GP;

}


