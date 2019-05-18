#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../Solvers/Solvers.h"

void Two_Steps_TG(Element ElementMesh, GaussPoint GP_Mesh,
		  Matrix Phi_n_Nod, Matrix Phi_n_GP,
		  Matrix A, int TimeStep)
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
   - ElementMesh -> Finite element mesh properties
   - GP_Mesh - > Gauss points mesh
   - A -> Matrix of the flux
   - Phi_n_GP -> Matrix with the values in the GP
   - TimeStep -> Time step

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

  /* Set the number of equations to solve */
  if(A.N_cols != A.N_rows){
    printf("Error in Two_Steps_TG_Mc() : The matrix A is not square !!! \n");
    exit(0);
  }
  if (Phi_n_GP.N_rows != A.N_cols){
    printf("Error in Two_Steps_TG_Mc() : The Phi and A, has not the same fields !!! \n");
  }
  int NumDOF = Phi_n_GP.N_rows;

  /* Fields defined in the Gauss-Points */
   /* Fields values in t = n */
  Matrix Phi_n12; /* Fields values in t = n + 1/2 */
  Matrix Phi_n1; /* Fields values in t = n + 1 */
  
  /* Flux fields defined inthe Gauss-Points*/
  Matrix Flux_n; /* Flux in t = n */
  Matrix Flux_n12; /* Flux in t = n + 1/2 */

  /* Elements properties */
  Matrix GP_i_XE; /* Element coordinates of the Gauss Points */
  Matrix F_Ref_GP; /* Reference deformation gradient evaluated in the GP */
  double J_GP; /* Jacobian evaluated in the Gauss-Point */
  Matrix dNdX_Ref_GP; /* Derivative Matrix */
  Matrix X_ElemNod = MatAlloc(2,1); /* Coordenates of the element nodes */
  int * Id_ElemNod; /* Pointer to the Connectivity of the element */
  Matrix M; /* Geometrical mass matrix */
  Matrix M_l; /* Geometrical mass matrix lumped */

  /* Nodal variables */
  Matrix Flux_n_Nod; /* Nodal values of the fluxes */
  Matrix RHS; /* Array with the right-hand side */
  Matrix RHS_i; /* Matrix with pointer functions to each field for the global RHS */
  Matrix DeltaPhiNod; /* Increment of the field solution in the nodes */
  Matrix DeltaPhiNod_i; /* Matrix of pointers to each field of the global DeltaPhinod */

  /* Boundary conditions variables */
  double BCC_val;
  int i_BCC_GP;

  /* Pointers fields for the input solver */
  DeltaPhiNod_i.nM = NULL;
  DeltaPhiNod_i.n = -999;
  DeltaPhiNod_i.N_rows = ElementMesh.NumNodesMesh; 
  DeltaPhiNod_i.N_cols = 1;  
  RHS_i.nM = NULL;
  RHS_i.n = -999;
  RHS_i.N_rows = ElementMesh.NumNodesMesh;
  RHS_i.N_cols = 1;

  printf(" * Get Gauss-Points-Fields values in t = n + 1/2 \n ");
  /**********************************************************************/
  /********* First step : Get GP-Fields values in t = n + 1/2  **********/
  /**********************************************************************/
  
  printf("\t --> Get the flux in the mesh for t = n \n");
  /****************** Get the flux in the mesh for t = n ****************/
  /**********************************************************************/
  /* F(Phi_n0)             F(Phi_n1)             F(Phi_n2)              */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n      */
  /*  Phi_n0     Phi_e0     Phi_n1     Phi_e1     Phi_n2                */
  /**********************************************************************/
  Flux_n_Nod = Scalar_prod(A,Phi_n_Nod);

  /* printf("Flux_n_Nod \n %f ; %f \n %f ; %f \n", */
  /* 	 Flux_n_Nod.nM[0][0],Flux_n_Nod.nM[1][0], */
  /* 	 Flux_n_Nod.nM[0][1],Flux_n_Nod.nM[1][1]); */
  

  printf("\t --> Get the fields values in the Gauss-Points for t = n + 1/2 \n");
  /*** Get the value of the fields in the Gauss points for t = n + 1/2 ***/
  /***********************************************************************/
  /*           F(Phi_e0)             F(Phi_e1)                           */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1/2 */
  /*           --->^<---             --->^<---                           */
  /*         __|   |   |__         __|   |   |__                         */
  /* F(Phi_n0)     |      F(Phi_n1)      |      F(Phi_n2)                */
  /*   [N0]------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n        */
  /*            Phi_e0                Phi_e1                             */
  /***********************************************************************/
  
  /* 1º Allocate an auxiar variable with the fiels of analysis */
  Phi_n12 = MatAllocZ(NumDOF,GP_Mesh.NumGP);  
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){

    /* 2º Get the element connectivity */
    Id_ElemNod = ElementMesh.Connectivity[i];
    
    /* 3º Get the nodal coordinates (take in to account that 
       in the C programming languaje the index starts in 0) */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates[Id_ElemNod[j] - 1][0];
    }
    
    /* 4º Get the gradient of the element shape functions evaluated 
       in the Gauss-Point position */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    F_Ref_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);

    /* 5º Iterate over the fields */
    for(int j = 0 ; j<NumDOF ; j++){
      /* 6º Add the therm of the previous step */
      Phi_n12.nM[j][i] = Phi_n_GP.nM[j][i];
      /* 7º Add the flux contributions iterating over the elements nodes */
      for(int k = 0 ; k<ElementMesh.NumNodesElem ; k++){
	Phi_n12.nM[j][i] -= 0.5*DeltaTimeStep*(double)1/F_Ref_GP.n*
	  Flux_n_Nod.nM[j][Id_ElemNod[k]-1]*dNdX_Ref_GP.nV[k];
      }      
    }    
  }

  /* printf("Phi_n12 \n %f ; %f \n %f ; %f \n", */
  /* 	 Phi_n12.nM[0][0],Phi_n12.nM[1][0], */
  /* 	 Phi_n12.nM[0][1],Phi_n12.nM[1][1]); */
  
  /* Once we have use it, we free the memory */
  free(Flux_n_Nod.nM);

  printf(" * Get Nodal-Fields values in t = n + 1 \n ");
  /**********************************************************************/
  /******** Second step : Get Nodal-Fields values in t = n + 1 **********/
  /**********************************************************************/

  printf("\t --> Get the flux for the Gauss-Points in t = n + 1/2 \n");
  /************** Get the flux for the GP for t = n + 1/2 ****************/
  /***********************************************************************/
  /*            F(Phi_e0)             F(Phi_e1)                          */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1/2 */
  /*  Phi_n0     Phi_e0     Phi_n1     Phi_e1     Phi_n2                 */
  /***********************************************************************/
  Flux_n12 = Scalar_prod(A,Phi_n12);
  
  /* printf("Flux_n12 \n %f ; %f \n %f ; %f \n", */
  /* 	 Flux_n12.nM[0][0],Flux_n12.nM[1][0], */
  /* 	 Flux_n12.nM[0][1],Flux_n12.nM[1][1]); */
  

  /* Once we have use it, we free the memory */
  free(Phi_n12.nM);

  printf("\t --> Get the RHS in the mesh in t = n + 1/2 \n");
  /************ Calcule the RHS in the mesh for t = n + 1/2 **************/
  /***********************************************************************/
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1   */
  /*         ··· SOLVER ···      ··· SOLVER ···                          */
  /*  RHS_n0               RHS_n1                RHS_n2                  */
  /*     ^                    ^                     ^                    */
  /*     |                    |                     |                    */
  /* F(Phi_n0) <F(Phi_e0)> F(Phi_n1) <F(Phi_e1)> F(Phi_n2)               */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1/2 */
  /***********************************************************************/
  
  /* 1º Allocate the right-hand side and initialize to zero */
  RHS = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh);
  
  /* Include the flux therms and the boundary conditions */
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){ /* Iterate over the GP */
    
    /* 2º Get the element connectivity */
    Id_ElemNod = ElementMesh.Connectivity[i];
    
    /* 3º Get the nodal coordinates (take in to account that 
	in the C programming languaje the index starts in 0)  */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates[Id_ElemNod[j] - 1][0];
    }
    
    /* 4º Get the GP coordiantes in the element */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    
    /* 5º Get the shape functions derivatives of the element evaluated in the GP */
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    
    /* 6º Get the deformation gradient of the element in the GP */
    F_Ref_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);
    F_Ref_GP.N_cols = 1, F_Ref_GP.N_rows = 1;
    
    /* 7º Get the Jacobian of the element in the GP */
    J_GP = Get_Determinant(F_Ref_GP);

    /* 8º Add the flux contribution to the RHS, note that each 
       GP contributes to each node of the element where it is */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumDOF ; k++){
	/* 9º If Get calcule the flux and add to the RHS */
	RHS.nM[k][Id_ElemNod[j]-1] +=
	  - dNdX_Ref_GP.nV[j]*(double)1/F_Ref_GP.n*fabs(J_GP)*Flux_n12.nM[k][i];
      }      
    }    
  }
    
  /* 10º Add the flux contributions in the boundary */
  for(int i = 0 ; i<ElementMesh.NumNodesBound ; i++){
    for(int j = 0 ; j<NumDOF ; j++){
      
      BCC_val = GetBoundaryCondition(ElementMesh.NodesBound[i]-1,j,TimeStep);

      if( (ElementMesh.NodesBound[i] - 1) == 0) i_BCC_GP = 0;
      if( (ElementMesh.NodesBound[i] - 1) == 20) i_BCC_GP = 19;

      if ( BCC_val != BCC_val ){
	RHS.nM[j][i] += -0.5*Flux_n12.nM[j][i_BCC_GP];
      }
      else {
	RHS.nM[j][i] += -1*BCC_val;
      }   
      
    }
  }
  
  /* 11º Multiply by the time step */
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumDOF ; j++){
      RHS.nM[j][i] *= -DeltaTimeStep;
    }
  }
  
  
  printf("\t --> Solve the sistem of equations to get U_n1 in the nodes \n");
  /******** Solve the sistem of equations to get U_n1 in the nodes *******/
  
  /* 1º Get the mass matrix */
  M = Get_Geom_Mass_Matrix(GP_Mesh,ElementMesh);
  M_l = Get_Lumped_Matrix(M);

  /* 2º Allocate the output array solution */
  DeltaPhiNod = MatAllocZ(2,ElementMesh.NumNodesMesh);
  
  /* 3º Iterate over the fields */
  for(int i = 0 ; i<NumDOF ; i++){
    /* 4º Asign values to the pointer to solve the equations */
    DeltaPhiNod_i.nV = DeltaPhiNod.nM[i];
    RHS_i.nV = RHS.nM[i];
    /* 5º Run the solver */
    /* DeltaPhiNod_i = Jacobi_Conjugate_Gradient_Method(M,RHS_i,DeltaPhiNod_i); */
    DeltaPhiNod_i = One_Iteration_Lumped(M_l,RHS_i,DeltaPhiNod_i);
    /* Note that we dont have to return to DeltaPhinod,
       because the pointer will do for us */
  }

  /* Once we have use it, we free the memory */
  free(RHS.nM);
  free(M.nM);
  free(M_l.nV);

  /* 6º Update the solution */
  for(int i = 0; i<ElementMesh.NumNodesMesh ; i++){
    for(int k=0 ; k<NumDOF ; k++){
      Phi_n_Nod.nM[k][i] += DeltaPhiNod.nM[k][i];
    }    
  }
  
  /* printf("DeltaPhiNod \n %f ; %f \n %f ; %f \n", */
  /* 	 DeltaPhiNod.nM[0][0],DeltaPhiNod.nM[1][0], */
  /* 	 DeltaPhiNod.nM[0][1],DeltaPhiNod.nM[1][1]); */
  
  /* printf("Phi_n_Nod \n %f ; %f \n %f ; %f \n", */
  /* 	 Phi_n_Nod.nM[0][0],Phi_n_Nod.nM[1][0], */
  /* 	 Phi_n_Nod.nM[0][1],Phi_n_Nod.nM[1][1]); */
  
  
  /* Once we have use it, we free the memory */
  free(DeltaPhiNod.nM);
  
  
} /* End of Two_Steps_TG_Mc() */


void Flux_Corrected_Transport(Element ElementMesh, GaussPoint GP_Mesh,
			      Matrix Phi_n_Nod, Matrix Phi_n_GP,
			      Matrix A, int TimeStep)
/* Flux-corrected transport solver, based on :
   "Finite element Flux-Corrected transport (FEM-FCT) 
   for the euler and navier-stokes equations"
   Rainald Löhner, Ken Morgan, Jaime Peraire and Medhi Vahdati
 */
{
  
} 
