#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ToolsLib/TypeDefinitions.h"
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
  Matrix J_GP; /* Jacobian of the element */
  Matrix dNdX_Ref_GP; /* Derivative Matrix */
  Matrix X_ElemNod; /* Coordenates of the element nodes */
  int * Id_ElemNod; /* Pointer to the Connectivity of the element */
  
  /* Auxiliar pointer to the fields */
  //Phi_n =  MatAlloc(2,GP_Mesh.NumGP);
  Phi_n.nM[0] = GP_Mesh.Phi.Stress.nV;
  Phi_n.nM[1] = GP_Mesh.Phi.vel.nV;

  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){
    printf("%f ; %f \n",GP_Mesh.Phi.Stress.nV[i],GP_Mesh.Phi.vel.nV[i]);
  }
  /* Allocate an auxiar variable with the fiels of analysis : 
   change this in the future for a clever fields selector...*/
  Phi_n12 = MatAlloc(2,GP_Mesh.NumGP);

  exit(0);
  /* Get the flux for the GP */
  Flux_n = Scalar_prod(A,Phi_n);

  /* First step : Get U_n12 */
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){

    /* Calcule the Jacobian of the element evaluated in the GP : */
    /* Get the element connectivity and take in to account that 
       in the C programming languaje the index starts in 0 */
    Id_ElemNod = ElementMesh.Connectivity[i];
    /* Get the coordinates of the nodes */
    X_ElemNod.nV[0] = ElementMesh.Coordinates[Id_ElemNod[0] - 1][0];
    X_ElemNod.nV[1] = ElementMesh.Coordinates[Id_ElemNod[1] - 1][0];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(X_ElemNod);
    J_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);
    
    Phi_n12.nM[0][i] += Phi_n.nM[0][i]; 
    Phi_n12.nM[1][i] += Phi_n.nM[1][i];

    Phi_n12.nM[0][i] -= 0.5*J_GP.n*Flux_n.nM[0][i]*(dNdX_Ref_GP.nV[0] +
						  dNdX_Ref_GP.nV[1]); 
    Phi_n12.nM[1][i] -= 0.5*J_GP.n*Flux_n.nM[1][i]*(dNdX_Ref_GP.nV[0] +
						  dNdX_Ref_GP.nV[1]); 
  }

  /* Second step : Get RHS and solve */

  
  
}
