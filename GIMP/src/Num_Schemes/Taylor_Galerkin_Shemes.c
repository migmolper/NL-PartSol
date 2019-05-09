#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../ElementsFunctions/ElementTools.h"

void Two_Steps_TG_Mc(Element ElementMesh, GaussPoint GP_Mesh, Matrix A)
/* 

   This sheme pretends to solve :
   Ui_t + Aij*Ui_x = 0, (1D)
   
   Where : 
   - Ui : array with the fields
   - ()_t : partial time derivative
   - ()_x : partial x derivative
   - Aij : Matrix with the physical properties, 
   depends of the case that we are solving

   Input parameters :
   - Element mesh
   - Gauss points mesh
   - "A" matrix

   Two steps Taylor-Galerkin using the consistent mass matrix :

   Employ linear element for the discretization via Bubnov-Galerkin 
   --> First step (advective predictor) :
   U_n12 = U_n - 0.5*DeltaT*Flux(U^a)_n
   --> Second step (corrector therm) :
   U_n1 - U_n = - DeltaT*Flux(U^a)_n12 + DeltaT*Flux(U^v)_n

   Solve the resulting system  with a Jacobi conjugate gradient with 
   at least three iterations :
   --> Solver :
   Mc * (U_n1 - U_n) = RHS_n
   Where :
   ------> Mc is the consistent mass matrix
   ------> RHS_n is the right-hand therm

   See on : J.Peraire, "A Finite Element Method for Convection dominated Flows",
   Ph.D. Thesis. University of Wales, Swansea (1986).
*/
{
  /* Variable definitions */
  Matrix Phi_n; /* Fields values in t = n */
  Matrix Phi_n12; /* Fields values in t = n + 1/2 */
  Matrix Flux_n; /* Flux in t = n */
  Matrix Flux_n12; /* Flux in t = n + 1/2 */
  Matrix J_GP; /* Jacobian of the element */
  Matrix dNdX_Ref_GP; /* Derivative Matrix */
  Matrix X_ElemNod = MatAlloc(2,1); /* Coordenates of the element nodes */
  int * Id_ElemNod; /* Pointer to the Connectivity of the element */
  Matrix M; /* Geometrical mass matrix */
  Matrix M_l; /* Lumped geometrical mass matrix */
  Matrix RHS; /* Array with the right-hand side */
  int NumEquations;
  /* Auxiliar pointer to the fields in t = n, in this special case, to avoid
   a second malloc that waste memory, we create a two rows table of pointer,
  so we have to fill the rest of the Matrix type field, in order to allow a 
  nide behaviour of the linear algebra functions */
  Phi_n.nM =  (double **)malloc((unsigned)2*sizeof(double *));
  Phi_n.nV = NULL; /* Set to NULL the (double *) pointer */
  Phi_n.n = -999; /* Set to -999 the scalar variable */
  Phi_n.N_rows = 2; /* Number of rows */
  Phi_n.N_cols = GP_Mesh.NumGP; /* Number of columns */
  Phi_n.nM[0] = GP_Mesh.Phi.Stress.nV; /* Asign to the first row the stress field */
  Phi_n.nM[1] = GP_Mesh.Phi.vel.nV; /* Asign to the first row the velocity field */

  /* Set the number of equations to solve */
  if(A.N_cols != A.N_rows){
    printf("Error in Two_Steps_TG_Mc() : The matrix A is not square \n");
    exit(0);
  }
  NumEquations = A.N_cols;
  
  /* Allocate an auxiar variable with the fiels of analysis : 
   change this in the future for a clever fields selector...*/
  Phi_n12 = MatAllocZ(2,GP_Mesh.NumGP);

  /* Get the flux for the GP in t = n */
  Flux_n = Scalar_prod(A,Phi_n);

  /* First step : Get U_n12 */
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){

    /* Calcule the Jacobian of the element evaluated in the GP : */
    /* Get the element connectivity and take in to account that 
       in the C programming languaje the index starts in 0 */
    Id_ElemNod = ElementMesh.Connectivity[i];

    /* Get the coordinates of the nodes */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates[Id_ElemNod[j] - 1][0];
    }
    
    dNdX_Ref_GP = ElementMesh.dNdX_ref(X_ElemNod);
    J_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);

    /* Get the value of U in t = n+1/2 */
    for(int j = 0 ; j<NumEquations ; j++){
      Phi_n12.nM[j][i] = Phi_n.nM[j][i] - 0.5*DeltaTimeStep*J_GP.n*
	Flux_n.nM[j][i]*(dNdX_Ref_GP.nV[0] + dNdX_Ref_GP.nV[1]); 
    }
  }

  /* Get the flux for the GP in t = n + 1/2 */
  Flux_n12 = Scalar_prod(A,Phi_n12); 
  
  /* Second step : Get RHS and solve */

  /* Allocate the right hand side */
  RHS = MatAllocZ(2,ElementMesh.NumNodesMesh);

  /* Include the flux therms */
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){ /* Iterate over the GP */

    /* Calcule the Jacobian of the element evaluated in the GP : */
    /* Get the element connectivity and take in to account that 
       in the C programming languaje the index starts in 0 */
    Id_ElemNod = ElementMesh.Connectivity[i];

    /* Get the coordinates of the nodes */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates[Id_ElemNod[j] - 1][0];
    }
    
    dNdX_Ref_GP = ElementMesh.dNdX_ref(X_ElemNod);
    J_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);

    /* Add the flux contribution to the RHS, note that each 
     GP contributes to each node of the element where it is */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumEquations ; k++){
	RHS.nM[k][Id_ElemNod[j]-1] -= dNdX_Ref_GP.nV[j]*
	(double)1/J_GP.n*
	fabs(J_GP.n)*
	Flux_n12.nM[j][i];
      }      
    }
    
  }

  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    printf("%f ; %f \n",RHS.nM[0][i],RHS.nM[1][i]);
  }

  /* Include the boundary conditions */
  
  /* Get the mass matrix */
  M = Get_Geom_Mass_Matrix(GP_Mesh,ElementMesh);
  /* Get the Lumped-Mass matrix */
  M_l = Get_Lumped_Matrix(M);

  
  
}
