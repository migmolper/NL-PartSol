#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ToolsLib/TypeDefinitions.h"
#include "ToolsLib/Utils.h"
#include "ToolsLib/ShapeFunctions.h"

/*********************************************************************/

void Initialize_GP(GaussPoint * GP,Vector * RefCoords)
/* In this function we Initialize and allocate all the variables of the Gauss Points,
   if they are necessary */
{

  /* Initialize_Element id of the GP */
  GP[0].id = 0;

  /* Id of the element where it is */
  GP[0].Element_id = 0; 

  /* Initialize kind of material */
  strcpy(GP[0].Material,"Elastic");

  /* Initialize and allocate position field (Vectorial) in global coordiantes
   and in element coordinates */
  GP[0].x_GC.Size = 2;
  GP[0].x_GC.n = (double *)Allocate_Array(2,sizeof(double));
  /* GP[0].x_GC.n = ... Initialize */
  
  GP[0].x_EC.Size = 2;
  GP[0].x_EC.n = (double *)Allocate_Array(2,sizeof(double));
  GP[0].x_EC.n = RefCoords[0].n;

  /* Allocate and Initialize the Velocity and acceleration field (Vectorial) */
  GP[0].v.Size = 2;
  GP[0].v.n = (double *)Allocate_Array(2,sizeof(double));
  /*  GP[0].v.n = ... Initialize */

  GP[0].a.Size = 2;
  GP[0].a.n = (double *)Allocate_Array(2,sizeof(double));
  /*  GP[0].v.n = ... Initialize */

  /* Allocate and Initialize the Stress and Strain fields (Tensorial) */
  GP[0].Stress.Size = 2;
  GP[0].Stress.n = (double **)Allocate_Matrix(2,2,sizeof(double));
  /* GP[0].Stress.n = ... Initialize */

  GP[0].Strain.Size = 2;
  GP[0].Strain.n = (double **)Allocate_Matrix(2,2,sizeof(double));
  /*GP[0].Strain.n = ...  Initialize */
  
  /* Initialize and allocate the deformation gradient of the reference element 
     It is a tensor field */
  GP[0].F_ref.Size = 2;
  GP[0].F_ref.n = (double **)Allocate_Matrix(2,2,sizeof(double));
  /*  GP[0].F_ref.n = ... Initialize */
  
  /* Initialize and allocate the deformation gradient of the element 
     It is a tensor field */
  GP[0].F.Size = 2;
  GP[0].F.n = (double **)Allocate_Matrix(2,2,sizeof(double));
  /*  GP[0].F.n = ... Initialize */

  /* Initialize and allocate the Lagrangian Cauchy-Green tensor (right) 
   of the element, it is a tensor field */
  GP[0].C.Size = 2;
  GP[0].C.n = (double **)Allocate_Matrix(2,2,sizeof(double));
  /*  GP[0].C.n = ... Initialize */

  /* Initialize and allocate the Eulerian Cauchy-Green tensor (left) 
   of the element, it is a tensor field */
  GP[0].B.Size = 2;
  GP[0].B.n = (double **)Allocate_Matrix(2,2,sizeof(double));
  /*  GP[0].B.n = ... Initialize */
  
}

/*********************************************************************/


/* Function for defining the element properties */
void Initialize_Element(Element * Elem,
			GaussPoint * GP_e,
			char * TypeElement,
			double ** GlobalNodsCoords,
			int * GlobalNodsId){
  
  if (strcmp(TypeElement,"Q4")==0){

    /* Identification number of the element */
    Elem[0].id = 0;

    /* Number of nodes in this element */
    Elem[0].NumberNodes = 4;

    /* Number of degree of freedom for each node */
    Elem[0].NumberDOF = 2;

    /* Global coordiantes of the nodes */
    Elem[0].X_g = GlobalNodsCoords;

    /* Conectivity mesh */
    Elem[0].N_id = GlobalNodsId;

    /* Shape functions and its derivatives */
    Elem[0].N_ref = Q4;
    Elem[0].dNdX_ref = dQ4;

    /* List of GP inside of the element */
    Elem[0].GP_e = GP_e;

    /* Auxiliar matrix to do the operations in K matrix */
    Elem[0].B.N_cols = 2;
    Elem[0].B.N_rows = 4;
    Elem[0].B.n = (double **)Allocate_Matrix(4,2,sizeof(double));
  }
  
}



/*********************************************************************/

void Get_RefDeformation_Gradient(GaussPoint * GP,
				 Element * Elem)
/*
  Get the deformation gradient of the reference element:

  F = x_{\alpha} \otimes grad(N_{\alpha})
  
  Inputs : 
  - GP -> GP structure
  - Elem -> Structure with the properties of the element

  Output :
  - F_ref -> Deformation gradient of the reference element
  evaluated in the GP
  - B_ref -> Ma
*/
{
  Element Elem_GP = Elem[GP[0].Element_id];
  
  /* Get the values of the shape functions derivatives */
  double ** dNdX_ref = (double **)Allocate_Matrix(4,2,sizeof(double));
  dNdX_ref = Elem_GP.dNdX_ref(&(GP[0].x_EC));

  for(int n = 0 ; n<Elem_GP.NumberNodes ; n++){
    /* Loop over the nodes of the element */    
    for(int i = 0 ; i<Elem_GP.NumberDOF ; i++){
      /* Loop over the degree of freedoms (Reference Element) */
      for(int j = 0 ; j<Elem_GP.NumberDOF ; j++){
	/* Loop over the degree of freedoms (Element) */

	/* In this case we use isoparametric formulation, */
	if(n==0){/* Initialize and increase */
	  GP[0].F_ref.nM[i][j] = 0;
	  GP[0].F_ref.nM[i][j] += Elem_GP.X_g.nM[n][j]*dNdX_ref[n][i];
	}
	else{ /* Increase */
	  GP[0].F_ref.nM[i][j] += Elem_GP.X_g.nM[n][j]*dNdX_ref[n][i];
	}
	
      }/* for j */
    }/* for i */
  }/* for n */

  /* Free the gradient of the shape functions */
  free(dNdX_ref);
  
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

Matrix Get_dNdx_matrix(Element * Elem,
		       GaussPoint * GP)
/* Get the derivative matrix of the shape functions in global coordiantes, it is
   the so called B matrix in the classical formulation of the finite element method 
   Inputs: Element, GP where to evaluate it
   Outputs : Matrix dNdX
*/
{

  /* Store the GP information in a auxiliary element */
  Element Elem_GP = Elem[GP[0].Element_id];
  
  /* Allocate the output */
  Matrix dNdX;
  dNdX.N_rows = Elem_GP.NumberNodes;
  dNdX.N_cols = Elem_GP.NumberDOF;
  dNdX.nM = (double **)Allocate_Matrix(dNdX.N_rows,
				      dNdX.N_cols,
				      sizeof(double));
  
  /* Get the values of the shape functions derivatives */
  double ** dNdX_ref = (double **)Allocate_Matrix(dNdX.N_rows,
						  dNdX.N_cols,
						  sizeof(double));
  dNdX_ref = Elem_GP.dNdX_ref(GP[0].x_EC);
  
  /* Get the inverse of the reference deformation gradient and transpose it */
  Matrix F_ref_m1T;
  double aux;
  F_ref_m1T = Get_Inverse(GP[0].F_ref);
  for(int i = 0 ; i < F_ref_m1T.N_rows ; i++ ){
    for(int j = i ; j < F_ref_m1T.N_cols ; j++){
      aux = F_ref_m1T.n[i][j];
      F_ref_m1T.nM[i][j] = F_ref_m1T.nM[j][i];
      F_ref_m1T.nM[j][i] = aux;
    }
  }

  /* Fill the matrix */
  for(int i = 0 ; i < dNdX.N_rows ; i++){ /* Iterate over the nodes */
    for(int j = 0 ; j < dNdX.N_cols ; j++){ /* Iterate over the degree of freedom */
      dNdX.nM[i][j] = 0; /* Fill first with a zero */
      
      for(int k = 0 ; k < dNdX.N_cols ; k++){
	dNdX.nM[i][j] += F_ref_m1T.nM[j][k]*dNdX_ref[i][j+k];
      }
      
    } /* for i */
  } /* for j */

  free(dNdX_ref);
  
  return dNdX;  
}


/*********************************************************************/

void Get_Lagrangian_CG_Tensor(GaussPoint * GP)
/* 
   The Right Cauchy Green Deformation Tensor :
   C = F^{T} \cdot F 

   Inputs : Gauss point
*/
{
  /* Check if we have a null matrix */
  if (GP[0].F.n == NULL){
    puts("Error in Get_Lagrangian_CG_Tensor : GP.F tensor = NULL");
    exit(0);
  }
  
  /* Get C */
  GP[0].C = Scalar_prod(Transpose_Mat(GP[0].F), /* F^T */
			GP[0].F); /* F */
  
}



/*********************************************************************/

void Get_Eulerian_CG_Tensor(GaussPoint * GP) 
/*
 The Left Cauchy Green Deformation Tensor :
 B = F \cdot F^{T} 

 Inputs : Gauss point 
*/
{  
  /* Check if we have a null matrix */
  if (GP[0].F.n == NULL){
    puts("Error in Get_Eulerian_CG_Tensor : GP.F tensor = NULL");
    exit(0);
  }
  
  /* Get B */
  GP[0].B = Scalar_prod(GP[0].F, /* F */
			Transpose_Mat(GP[0].F)); /* F^T */
  
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

