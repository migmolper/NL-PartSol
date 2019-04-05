#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ShapeFunctions.h"
//#include "Constitutive.h"



/* Function for defining the element properties */
void Initialize_Element(Element * Elem,
			char * TypeElement,
			double ** GlobalNodsCoords,
			int * GlobalNodsId){
  
  if (strcmp(TypeElement,"Q4")==0){
    Elem->id = 0;
    Elem->NumberNodes = 4;
    Elem->NumberDOF = 2;
    Elem->X_g = GlobalNodsCoords;
    Elem->N_id = GlobalNodsId;
    Elem->n = Q4;
    Elem->dn = dQ4;
  }
}

void Initialize_GP(GaussPoint * GP,Vector * RefCoords){
  GP->x_EC = RefCoords;
}

void Get_Deformation_Gradient(GaussPoint * GP,
			      Element * Elem,
			      double ** B_ref,
			      double ** F)
/*
  Get the deformation gradient of an element:

  F = x_{\alpha} \otimes grad(N_{\alpha})
  
  Inputs : 
  - X_e -> Coordenates of a GP in the reference element
  - Elem -> Structure with the properties of the element

  Output :
  - F -> Deformation gradient evaluated in the GP
*/
{
  /* Get the values of the shape functions derivatives */
  B_ref = Elem->dn(GP->x_EC);

  for(int n = 0 ; n<Elem->NumberNodes ; n++){
    /* Loop over the nodes of the element */    
    for(int i = 0 ; i<Elem->NumberDOF ; i++){
      /* Loop over the degree of freedoms (Reference Element) */
      for(int j = 0 ; j<Elem->NumberDOF ; j++){
	/* Loop over the degree of freedoms (Element) */

	/* In this case we use isoparametric formulation, */
	if(n==0){/* Initialize and increase */
	  F[i][j] = 0;
	  F[i][j] += Elem->X_g[n][j]*B_ref[n][i];
	}
	else{ /* Increase */
	  F[i][j] += Elem->X_g[n][j]*B_ref[n][i];
	}
	
      }/* if j */
    }/* if i */
  }/* if n */
  
}


void Get_Lagrangian_CG_Tensor(GaussPoint * GP,Element * Elem, double ** C)
  /* The Right Cauchy Green Deformation Tensor :
   C = F^{T} \cdot F */
{
  double ** F = (double **)Allocate_Matrix(Elem->NumberDOF,
					   Elem->NumberDOF,
					   sizeof(double));

  /* Calcule the deformation gradient F to get C */
  Get_Deformation_Gradient(GP,Elem,F);

  /* Get C */
  for(int i = 0; i<Elem->NumberDOF ; i++){
    for(int j = 0; j<Elem->NumberDOF; j++){
      C[i][j] = F[j][i] * F[i][j];
      /* C    =  F^{T}  *   F  */
    }
  }

  /* Free memory */
  free(F);
}



void Get_Eulerian_CG_Tensor(GaussPoint * GP,Element * Elem,double ** B)
  /* The Left Cauchy Green Deformation Tensor :
   B = F \cdot F^{T} */
{
  double ** F = (double **)Allocate_Matrix(Elem->NumberDOF,
					   Elem->NumberDOF,
					   sizeof(double));
  
  /* Calcule the deformation gradient to get B */
  Get_Deformation_Gradient(GP,Elem,F);

  /* Get B */
  for(int i = 0; i<Elem->NumberDOF ; i++){
    for(int j = 0; j<Elem->NumberDOF; j++){
      B[i][j] = F[i][j] * F[j][i];
      /* B    =    F    *  F^{T}  */
    }
  }
  
  /* Free memory */
  free(F);
}

void Get_Deformation_Gradient_GP(Element * Elem,
				GaussPoint * GP,
				double ** B)
{
  double ** F = (double **)Allocate_Matrix(2,2,sizeof(double));
  double ** F_m1 = (double **)Allocate_Matrix(2,2,sizeof(double));

  Get_Deformation_Gradient(GP,Elem,F); /* Deformation gradient F */
  Get_Inverse(F,F_m1,2); /* Get the inverse of the deformation gradient */
    
  for(int i = 0 ; i<Element.NumberNodes ; i++){/* Iterate over the nodes in the Element to get B */

    
    
  }/* Nodes loop */

  /* Once we have calculated the deformation gradient 
     matrix we dont need anymore the deformation gradient
     and it inverse, so we free memory */
  free(F);
  free(F_m1);
  
}


double ** Get_Stiffness_Matrix(double ** D, /* Constitutive matrix */
			       GaussPoint * GP, /* Gauss points of the element */
			       Element * Elem) /* Element to analyze */
/* Calcule the element stiffness matrix K */
{

  /* We need the deformation grafient F and later it inverse */
  
  double ** B = (double **)Allocate_Matrix(Element->NumberNodes,
					   Element->NumberDOF,
					   sizeof(double));

  for(int i = 0 ; i<Element.NumberGP ; i++){/* Iterate over the Gauss Points */

    
    for(int j = 0 ; j<Element.NumberDOF ; j++){/* Iterate over the degree of freedom */
      
      
    }/* GP loop */
    
  }/* DOF loop */


  free(F);
  free(F_m1);
  
  return K;
}
