#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../GaussPointsFunctions/GaussPointsTools.h"
#include "../Solvers/Solvers.h"
#include "../Boundary_Conditions/BoundaryConditions.h"

void Two_Steps_TG_FEM(Mesh ElementMesh, GaussPoint GP_Mesh,
		  Matrix Phi_n_Nod, Matrix Phi_n_GP,
		  Matrix A, Matrix M, int TimeStep)
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

  /* Fields defined in the Gauss-Points in t = n + 1/2 */
  Matrix Phi_n12_GP; 
  
  /* Flux fields defined inthe Gauss-Points*/
  Matrix Flux_n12_GP; /* Flux in t = n + 1/2 */

  /* Elements properties */
  Matrix GP_i_XE; /* Element coordinates of the Gauss Points */
  Matrix F_Ref_GP; /* Reference deformation gradient evaluated in the GP */
  double J_GP; /* Jacobian evaluated in the Gauss-Point */
  Matrix dNdX_Ref_GP; /* Derivative Matrix */
  double dNdX_GP; /* Deformation gradient */
  Matrix X_ElemNod = MatAlloc(2,1); /* Coordenates of the element nodes */
  int Elem_GP_i;
  int * Id_ElemNod; /* Pointer to the Connectivity of the element */

  /* Nodal variables */
  Matrix Flux_n_Nod; /* Nodal values of the fluxes */
  Matrix RHS; /* Array with the right-hand side */
  Matrix RHS_i; /* Matrix with pointer functions to each field for the global RHS */
  Matrix DeltaPhiNod; /* Increment of the field solution in the nodes */
  Matrix DeltaPhiNod_i; /* Matrix of pointers to each field of the global DeltaPhinod */

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
  Phi_n12_GP = MatAllocZ(NumDOF,GP_Mesh.NumGP);  
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){

    /* 2º Index of the element where the G-P is located */
    Elem_GP_i = GP_Mesh.Element_id[i];

    /* 3º Get the element connectivity */
    Id_ElemNod = ElementMesh.Connectivity[Elem_GP_i];
    
    /* 4º Get the nodal coordinates (take in to account that 
       in the C programming languaje the index starts in 0) */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates.nM[Id_ElemNod[j] - 1][0];
    }
    
    /* 5º Get the gradient of the element shape functions evaluated 
       in the Gauss-Point position */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    F_Ref_GP = Get_RefDeformation_Gradient_L2(X_ElemNod,dNdX_Ref_GP);

    /* 6º Iterate over the fields */
    for(int j = 0 ; j<NumDOF ; j++){
      /* 7º Add the therm of the previous step */
      Phi_n12_GP.nM[j][i] = Phi_n_GP.nM[j][i];
      for(int k = 0 ; k<ElementMesh.NumNodesElem ; k++){
	/* 8º Get the deformation gradient */
	dNdX_GP = (double)1/F_Ref_GP.n*dNdX_Ref_GP.nV[k];
	/* 9º Add the flux contributions */
	Phi_n12_GP.nM[j][i] -= 0.5*DeltaTimeStep*
	  dNdX_GP*Flux_n_Nod.nM[j][Id_ElemNod[k]-1];
      }      
    }    
  }

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
  
  Flux_n12_GP = Scalar_prod(A,Phi_n12_GP);

  /* Once we have use it, we free the memory */
  free(Phi_n12_GP.nM);

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
      X_ElemNod.nV[j] = ElementMesh.Coordinates.nM[Id_ElemNod[j] - 1][0];
    }
    
    /* 4º Get the GP coordiantes in the element */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    
    /* 5º Get the reference shape functions derivatives of the element evaluated in the GP */
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    
    /* 6º Get the reference deformation gradient of the element in the GP */
    F_Ref_GP = Get_RefDeformation_Gradient_L2(X_ElemNod,dNdX_Ref_GP);
    F_Ref_GP.N_cols = 1, F_Ref_GP.N_rows = 1;
    
    /* 7º Get the Jacobian of the element in the GP */
    J_GP = Get_Determinant(F_Ref_GP);

    /* 8º Add the flux contribution to the RHS, note that each 
       GP contributes to each node of the element where it is */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      /* 9º Get the deformation gradient */
      dNdX_GP = dNdX_Ref_GP.nV[j]*(double)1/F_Ref_GP.n;
      for(int k = 0 ; k<NumDOF ; k++){
	/* 10º Calcule the flux and add to the RHS */
	RHS.nM[k][Id_ElemNod[j]-1] +=
	  - dNdX_GP*fabs(J_GP)*Flux_n12_GP.nM[k][i];
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
  
  /* 1º Allocate the output array solution */
  DeltaPhiNod = MatAllocZ(2,ElementMesh.NumNodesMesh);
  
  /* 2º Iterate over the fields */
  for(int i = 0 ; i<NumDOF ; i++){
    /* 3º Asign values to the pointer to solve the equations */
    DeltaPhiNod_i.nV = DeltaPhiNod.nM[i];
    RHS_i.nV = RHS.nM[i];
    /* 4º Run the solver */
    /* DeltaPhiNod_i = Jacobi_Conjugate_Gradient_Method(M,RHS_i,DeltaPhiNod_i); */
    DeltaPhiNod_i = One_Iteration_Lumped(M,RHS_i,DeltaPhiNod_i);
    /* Note that we dont have to return to DeltaPhinod,
       because the pointer will do for us */
  }
  
  /* Once we have use it, we free the memory */
  free(RHS.nM);

  /* 5º Update the solution */
  for(int i = 0; i<ElementMesh.NumNodesMesh ; i++){
    for(int k=0 ; k<NumDOF ; k++){
      Phi_n_Nod.nM[k][i] += DeltaPhiNod.nM[k][i];
    }    
  }

  /* Once we have use it, we free the memory */
  free(DeltaPhiNod.nM);

  /* 6º Add the Boundary conditions contributions in the boundary */
  /* ApplyBoundaryCondition_Nod(Phi_n_Nod,TimeStep);   */

  /* 7º Transfer information from the mesh to the GP for MPM */
  MeshToGaussPoints(ElementMesh,GP_Mesh,Phi_n_Nod,Phi_n_GP,M);
  
    
} /* End of Two_Steps_TG_Mc() */


void Two_Steps_TG_MPM(Mesh ElementMesh, GaussPoint GP_Mesh,
		      Matrix Phi_n_GP,Matrix M_l, Matrix A, int TimeStep){

  Matrix Nod_Int_Forces = MatAllocZ(2,ElementMesh.NumNodesMesh);
  Matrix Phi_n_Nod;
  int Elem_GP_i;
  int * Id_ElemNod;
  Matrix X_ElemNod = MatAlloc(2,1);
  Matrix GP_i_XE;
  Matrix dNdX_Ref_GP;
  Matrix F_Ref_GP;
  double dNdX_GP;
  Matrix Flux_n_Nod;
  Matrix Phi_n12_GP; 

  Phi_n_Nod = MatAllocZ(2,ElementMesh.NumNodesMesh);
  
  GaussPointsToMesh(ElementMesh,GP_Mesh,Phi_n_GP,Phi_n_Nod,M_l);

  Flux_n_Nod = Scalar_prod(A,Phi_n_Nod);

  Phi_n12_GP = MatAllocZ(2,GP_Mesh.NumGP);  
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){

    Elem_GP_i = GP_Mesh.Element_id[i];

    Id_ElemNod = ElementMesh.Connectivity[Elem_GP_i];
    
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates.nM[Id_ElemNod[j] - 1][0];
    }
    
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    F_Ref_GP = Get_RefDeformation_Gradient_L2(X_ElemNod,dNdX_Ref_GP);

    /* 6º Iterate over the fields */
    for(int j = 0 ; j<2 ; j++){
      /* 7º Add the therm of the previous step */
      Phi_n12_GP.nM[j][i] = Phi_n_GP.nM[j][i];
      for(int k = 0 ; k<ElementMesh.NumNodesElem ; k++){
	/* 8º Get the deformation gradient */
	dNdX_GP = (double)1/F_Ref_GP.n*dNdX_Ref_GP.nV[k];
	/* 9º Add the flux contributions */
	Phi_n12_GP.nM[j][i] -= 0.5*DeltaTimeStep*
	  dNdX_GP*Flux_n_Nod.nM[j][Id_ElemNod[k]-1];
      }      
    }    
  }

  /* Once we have use it, we free the memory */
  free(Flux_n_Nod.nM);  

  /* Get internal forces and external forces */
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){
    Elem_GP_i = GP_Mesh.Element_id[i];
    Id_ElemNod = ElementMesh.Connectivity[Elem_GP_i];
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates.nM[Id_ElemNod[j] - 1][0];
    }
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    F_Ref_GP = Get_RefDeformation_Gradient_L2(X_ElemNod,dNdX_Ref_GP);

    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      dNdX_GP = (double)1/F_Ref_GP.n*dNdX_Ref_GP.nV[j];
      Nod_Int_Forces.nM[1][Id_ElemNod[j]-1] += dNdX_GP*Phi_n12_GP.nM[0][Id_ElemNod[j]-1];
      Nod_Int_Forces.nM[0][Id_ElemNod[j]-1] += dNdX_GP*Phi_n12_GP.nM[1][Id_ElemNod[j]-1];
    }
      
  }
  
  free(Phi_n12_GP.nM);

  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    Phi_n_Nod.nM[1][i] += (double)1/M_l.nV[i]*DeltaTimeStep*(-Nod_Int_Forces.nM[1][i]);
    Phi_n_Nod.nM[0][i] += (double)1/M_l.nV[i]*DeltaTimeStep*(-Nod_Int_Forces.nM[0][i]);
  }

  /* ApplyBoundaryCondition_Nod(Phi_n_Nod,TimeStep); */

  MeshToGaussPoints(ElementMesh,GP_Mesh,Phi_n_Nod,Phi_n_GP,M_l);

  free(Phi_n_Nod.nM);
  free(Nod_Int_Forces.nM);


}
