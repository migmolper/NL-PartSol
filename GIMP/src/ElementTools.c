#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "include/Utils.h"
#include "include/ShapeFunctions.h"


/* Function for defining the element properties */
void Initialize_Element(Element * Elem,
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
    Elem[0].n = Q4;
    Elem[0].dn = dQ4;

    /* Auxiliar matrix to do the operations in K matrix */
    Elem[0].B.N_cols = 2;
    Elem[0].B.N_rows = 4;
    Elem[0].B.n = (double **)Allocate_Matrix(4,2,sizeof(double));
  }
  
}

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
  double ** deriv_SHPF = (double **)Allocate_Matrix(4,2,sizeof(double));
  deriv_SHPF = Elem_GP.dn(&(GP[0].x_EC));

  for(int n = 0 ; n<Elem_GP.NumberNodes ; n++){
    /* Loop over the nodes of the element */    
    for(int i = 0 ; i<Elem_GP.NumberDOF ; i++){
      /* Loop over the degree of freedoms (Reference Element) */
      for(int j = 0 ; j<Elem_GP.NumberDOF ; j++){
	/* Loop over the degree of freedoms (Element) */

	/* In this case we use isoparametric formulation, */
	if(n==0){/* Initialize and increase */
	  GP[0].F_ref.n[i][j] = 0;
	  GP[0].F_ref.n[i][j] += Elem_GP.X_g[n][j]*deriv_SHPF[n][i];
	}
	else{ /* Increase */
	  GP[0].F_ref.n[i][j] += Elem_GP.X_g[n][j]*deriv_SHPF[n][i];
	}
	
      }/* for j */
    }/* for i */
  }/* for n */

  /* Free the gradient of the shape functions */
  free(deriv_SHPF);
  
}


/*********************************************************************/

/* void Get_Deformation_Gradient(Element * Elem, /\* Gauss point *\/ */
/* 			      GaussPoint * GP) /\* Element *\/ */
/* { */
/*   double ** F_ref_m1 = (double **)Allocate_Matrix(2,2,sizeof(double)); */

/*   /\* Deformation reference gradient F_ref *\/ */
/*   Get_RefDeformation_Gradient(GP,Elem); */

/*   /\* Get the inverse of the deformation gradient *\/ */
/*   Get_Inverse(GP[0].F_ref.n,F_ref_m1,2); */

/*   /\* Iterate over the nodes in the Element to get B *\/ */
/*   for(int i = 0 ; i<Element[0].NumberNodes ; i++){ */

    
    
/*   }/\* Nodes loop *\/ */

/*   /\* Once we have calculated the deformation gradient */
/*      matrix we dont need anymore the inverse of the */
/*      reference deformation gradient so we free memory *\/   */
/*   free(F_ref_m1); */
  
/* } */


/*********************************************************************/

void Get_Lagrangian_CG_Tensor(GaussPoint * GP) /* Gauss point */
			        
/* The Right Cauchy Green Deformation Tensor :
   C = F^{T} \cdot F 
*/
{
  /* Get C */
  for(int i = 0; i<GP[0].F.Size ; i++){
    for(int j = 0; j<GP[0].F.Size; j++){
      GP[0].C.n[i][j] = GP[0].F.n[j][i] * GP[0].F.n[i][j];
      /*    C         =       F^{T}     *       F  */
    } /* for j */
  } /* for i */
  
}



/*********************************************************************/

void Get_Eulerian_CG_Tensor(GaussPoint * GP) /* Gauss point */
  
  /* The Left Cauchy Green Deformation Tensor :
   B = F \cdot F^{T} */
{
  /* Get B */
  for(int i = 0; i<GP[0].F.Size ; i++){
    for(int j = 0; j<GP[0].F.Size; j++){
      GP[0].B.n[i][j] = GP[0].F.n[i][j] * GP[0].F.n[j][i];
      /*    B         =       F         *       F^{T}  */
    } /* for j */
  } /* for i */
  
}


/*********************************************************************/


/* double ** Get_Stiffness_Matrix(double ** D, /\* Constitutive matrix *\/ */
/* 			       GaussPoint * GP, /\* Gauss points of the element *\/ */
/* 			       Element * Elem) /\* Element to analyze *\/ */
/* /\* Calcule the element stiffness matrix K *\/ */
/* { */

/*   /\* We need the deformation grafient F and later it inverse *\/ */
  
/*   double ** B = (double **)Allocate_Matrix(Element->NumberNodes, */
/* 					   Element->NumberDOF, */
/* 					   sizeof(double)); */

/*   for(int i = 0 ; i<Element.NumberGP ; i++){/\* Iterate over the Gauss Points *\/ */

    
/*     for(int j = 0 ; j<Element.NumberDOF ; j++){/\* Iterate over the degree of freedom *\/ */
      
      
/*     }/\* GP loop *\/ */
    
/*   }/\* DOF loop *\/ */


/*   free(F); */
/*   free(F_m1); */
  
/*   return K; */
/* } */

