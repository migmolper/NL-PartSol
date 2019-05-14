#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../ToolsLib/GlobalVariables.h"
#include "../ToolsLib/Utils.h"
#include "../ElementsFunctions/ElementTools.h"
#include "../Solvers/Solvers.h"

void Two_Steps_TG_Mc(Element ElementMesh, GaussPoint GP_Mesh,
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
/******************************************************************/
/*                                                                */
/*  Phi_n --> Phi_e <-- Phi_n --> Phi_e <-- Phi_n                 */
/*  [N0]------(e0)------[N1]------(e1)------[N2]-- ;  t = n + 1   */
/*     ^                ^  ^                ^                     */
/*     |___          ___|  |___          ___|                     */
/*          \       /          \        /                         */
/*           F(Phi)_e           F(Phi)_e                          */
/*  [N0]------(e0)------[N1]------(e1)------[N2]-- ;  t = n + 1/2 */
/*             ^                   ^                              */
/*             |                   |                              */
/*  [N0]------(e0)------[N1]------(e1)------[N2]-- ;  t = n       */
/*                                                                */
/******************************************************************/
  
{
  /* Variable definitions */

  /* Set the number of equations to solve */
  if(A.N_cols != A.N_rows){
    printf("Error in Two_Steps_TG_Mc() : The matrix A is not square \n");
    exit(0);
  }
  int NumDOF = A.N_cols;

  /* Fields defined in the Gauss-Points */
  Matrix Phi_n; /* Fields values in t = n */
  Matrix Phi_n12; /* Fields values in t = n + 1/2 */
  Matrix Phi_n1; /* Fields values in t = n + 1 */
  
  /* Flux fields defined inthe Gauss-Points*/
  Matrix Flux_n; /* Flux in t = n */
  Matrix Flux_n12; /* Flux in t = n + 1/2 */

  /* Elements properties */
  Matrix F_Ref_GP; /* Reference deformation gradient evaluated in the GP */
  double J_GP; /* Jacobian evaluated in the Gauss-Point */
  Matrix dNdX_Ref_GP; /* Derivative Matrix */
  Matrix X_ElemNod = MatAlloc(2,1); /* Coordenates of the element nodes */
  int * Id_ElemNod; /* Pointer to the Connectivity of the element */
  Matrix M; /* Geometrical mass matrix */

  /* Nodal variables */
  Matrix Phi_n_Nod; /* Nodal values of each field */
  Matrix Flux_n_Nod; /* Nodal values of the fluxes */
  Matrix RHS; /* Array with the right-hand side */
  Matrix RHS_i; /* Matrix with pointer functions to each field for the global RHS */
  Matrix DeltaPhiNod; /* Increment of the field solution in the nodes */
  Matrix DeltaPhiNod_i; /* Matrix with pointer functions to each field of the global DeltaPhinod */

  /* Boundary conditions variables */
  double BCC_GP_DOF;
  int i_BCC;  

  /* Auxiliar variable to the transference from the nodes to the Gauss Points */
  int Elem_GP_i; /* Element where the GP is placed */
  int * Nodes_Elem_GP_i; /* Nodes of the element where the GP is placed */
  Matrix GP_i_XE; /* Element coordinates of the Gauss Points */
  Matrix N_ref_XG; /* Element function evaluated in the Gauss-Points */
  
  /* Auxiliar pointer to the fields in t = n, in this special case, to avoid
   a second malloc that waste memory, we create a two rows table of pointer,
  so we have to fill the rest of the Matrix type field, in order to allow a 
  nide behaviour of the linear algebra functions */
  Phi_n.nM =  (double **)malloc((unsigned)NumDOF*sizeof(double *));
  Phi_n.nV = NULL; /* Set to NULL the (double *) pointer */
  Phi_n.n = -999; /* Set to -999 the scalar variable */
  Phi_n.N_rows = NumDOF; /* Number of rows */
  Phi_n.N_cols = GP_Mesh.NumGP; /* Number of columns */
  Phi_n.nM[0] = GP_Mesh.Phi.Stress.nV; /* Asign to the first row the stress field */
  Phi_n.nM[1] = GP_Mesh.Phi.vel.nV; /* Asign to the first row the velocity field */

  /* Pointers fields for the input solver */
  DeltaPhiNod_i.nM = NULL;
  DeltaPhiNod_i.n = -999;
  DeltaPhiNod_i.N_rows = ElementMesh.NumNodesMesh; 
  DeltaPhiNod_i.N_cols = 1;  
  RHS_i.nM = NULL;
  RHS_i.n = -999;
  RHS_i.N_rows = ElementMesh.NumNodesMesh; 
  RHS_i.N_cols = 1;


  printf(" * First step : Get U_n12 \n ");
  /**********************************************************************/
  /*********************** First step : Get U_n12 ***********************/
  /**********************************************************************/

  printf("\t --> Transfer the GP information to the mesh to calculate the gradients \n");
  /* Transfer the GP information to the mesh to calculate the gradients */
  /**********************************************************************/
  /*  [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n       */
  /* Phi_n0 <-- Phi_e0 --> Phi_n1 <-- Phi_e1 --> Phi_n2                 */
  /**********************************************************************/
  /* 1º Allocate the variable to store the nodal information */
  Phi_n_Nod = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh);
  for(int i=0 ; i<GP_Mesh.NumGP ; i++){
    
    /* 2º Index of the element where the G-P is located */
    Elem_GP_i = GP_Mesh.Element_id[i];
    
    /* 3º List of nodes of the element */
    Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i];
    
    /* 4º Element coordinates of the G-P */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    
    /* 5º Evaluate the shape functions of the element in the 
       coordinates of the G-P */
    N_ref_XG = ElementMesh.N_ref(GP_i_XE);
    
    /* 6º Iterate over the nodes of the element */
    for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){
      /* 7º Iterate over the fields */
      for(int k=0 ; k<NumDOF ; k++){
  	/* 8º Update the field value */;
  	Phi_n_Nod.nM[k][Nodes_Elem_GP_i[j]-1] += N_ref_XG.nV[j]*Phi_n.nM[k][i];
      }
    }
  }

  printf("\t --> Get the flux for the nodes in t = n \n");
  /****************** Get the flux for the nodes in t = n ***************/
  /**********************************************************************/
  /* F(Phi_n0)             F(Phi_n1)             F(Phi_n2)              */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n      */
  /*  Phi_n0     Phi_e0     Phi_n1     Phi_e1     Phi_n2                */
  /**********************************************************************/
  Flux_n_Nod = Scalar_prod(A,Phi_n_Nod);  

  printf("\t --> Get the value of the fields in the Gauss points in t = n + 1/2 \n");
  /*** Get the value of the fields in the Gauss points in t = n + 1/2 ****/
  /***********************************************************************/
  /*            F(Phi)_e              F(Phi)_e                           */
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
      Phi_n12.nM[j][i] = Phi_n.nM[j][i];
      /* 7º Check if this DOF of this Gauss-Point has some kind of restriction */
      BCC_GP_DOF = GetBoundaryCondition(i,j,TimeStep);
      if(BCC_GP_DOF != BCC_GP_DOF){
	/* 8aº Add the flux contributions iterating over the elements nodes */
	for(int k = 0 ; k<ElementMesh.NumNodesElem ; k++){
	  Phi_n12.nM[j][i] -= 0.5*DeltaTimeStep*(double)1/F_Ref_GP.n*
	    Flux_n_Nod.nM[j][Id_ElemNod[k]-1]*dNdX_Ref_GP.nV[k];
	}
      }
      else{
	/* 8bº Set the value */
	Phi_n12.nM[j][i] = BCC_GP_DOF;
      }
      
    }
    
  }
  
  /* Once we have use it, we free the memory */
  free(Flux_n_Nod.nM);

  printf(" * Second step : Get U_n1 \n ");
  /**********************************************************************/
  /*********************** Second step : Get U_n1 ***********************/
  /**********************************************************************/

  printf("\t --> Get the flux for the GP in t = n + 1/2 \n");
  /*************** Get the flux for the GP in t = n + 1/2 ****************/
  /***********************************************************************/
  /* F(Phi_n0)             F(Phi_n1)             F(Phi_n2)               */
  /*   [N0]-------(e0)-------[N1]-------(e1)-------[N2]-- ;  t = n + 1/2 */
  /*  Phi_n0     Phi_e0     Phi_n1     Phi_e1     Phi_n2                 */
  /***********************************************************************/
  Flux_n12 = Scalar_prod(A,Phi_n12);

  /* Once we have use it, we free the memory */
  free(Phi_n12.nM);
      
  /* Allocate the right-hand side and initialize to zero */
  RHS = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh);
  
  /* Include the flux therms and the boundary conditions */
  for(int i = 0 ; i<GP_Mesh.NumGP ; i++){ /* Iterate over the GP */
    /* Calcule the reference deformation gradient of the element evaluated 
       in the Gauss-Point :
        - Get the element connectivity 
	- Get the nodal coordinates (take in to account that 
	in the C programming languaje the index starts in 0) 
	- Get the gradient of the element shape functions evaluated 
	in the Gauss-Point position
    */
    Id_ElemNod = ElementMesh.Connectivity[i];
    /* Get the coordinates of the nodes */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      X_ElemNod.nV[j] = ElementMesh.Coordinates[Id_ElemNod[j] - 1][0];
    }
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    dNdX_Ref_GP = ElementMesh.dNdX_ref(GP_i_XE);
    F_Ref_GP = Get_RefDeformation_Gradient(X_ElemNod,dNdX_Ref_GP);
    F_Ref_GP.N_cols = 1;
    F_Ref_GP.N_rows = 1;
    /* Get the Jacobian */
    J_GP = Get_Determinant(F_Ref_GP);

    /* Add the flux contribution to the RHS, note that each 
     GP contributes to each node of the element where it is */
    for(int j = 0 ; j<ElementMesh.NumNodesElem ; j++){
      for(int k = 0 ; k<NumDOF ; k++){
	/* Check if this Gauss-Point has a set value in this DOF */
	BCC_GP_DOF = GetBoundaryCondition(i,k,TimeStep);
	if(BCC_GP_DOF != BCC_GP_DOF){
	  /* If not, get calcule the flux and add to the RHS */
	  RHS.nM[k][Id_ElemNod[j]-1] +=
	    dNdX_Ref_GP.nV[j]*(double)1/F_Ref_GP.n*fabs(J_GP)*Flux_n12.nM[k][i];
	}
	else{
	  /* If yes, calcule the RHS taking this into account */
	  RHS.nM[k][Id_ElemNod[j]-1] += BCC_GP_DOF;
	}
	
      }      
    }    
  }
  /* Multiply by the time step */
  for(int i = 0 ; i<ElementMesh.NumNodesMesh ; i++){
    for(int j = 0 ; j<NumDOF ; j++){
      RHS.nM[j][i] *= -DeltaTimeStep;
    }
  }

  /* Solve the sistem of equations */
  /* Get the mass matrix */
  M = Get_Geom_Mass_Matrix(GP_Mesh,ElementMesh); 
  /* Allocate the array with the solutions */
  DeltaPhiNod = MatAllocZ(NumDOF,ElementMesh.NumNodesMesh);
  /* Solve both all the systems of equations */
  for(int i = 0 ; i<NumDOF ; i++){
    /* Asign values to the pointer to solve the equations */
    DeltaPhiNod_i.nV = DeltaPhiNod.nM[i];
    RHS_i.nV = RHS.nM[i];
    /* Run the solver */
    /* DeltaPhiNod_i = Jacobi_Conjugate_Gradient_Method(M,RHS_i,DeltaPhiNod_i); */
    DeltaPhiNod_i = One_Iteration_Lumped(M,RHS_i,DeltaPhiNod_i);
    /* Note that we dont have to return to DeltaPhinod,
       because the pointer will do for us */
  }

  /* Once we have use it, we free the memory */
  free(RHS.nM);
  free(M.nM);

  for(int i = 0; i<ElementMesh.NumNodesMesh ; i++){
    for(int k=0 ; k<NumDOF ; k++){
      Phi_n_Nod.nM[k][i] += DeltaPhiNod.nM[k][i];
    }
  }

  
  /* Transfer the nodal values to the G-P */
  for(int i=0 ; i<GP_Mesh.NumGP ; i++){
    /* Index of the element where the G-P is located */
    Elem_GP_i = GP_Mesh.Element_id[i];
    /* List of nodes of the element */
    Nodes_Elem_GP_i = ElementMesh.Connectivity[Elem_GP_i];
    /* Element coordinates of the G-P */
    GP_i_XE.n = GP_Mesh.Phi.x_EC.nV[i];
    /* Evaluate the shape functions of the element in the coordinates of the G-P */
    N_ref_XG = ElementMesh.N_ref(GP_i_XE);
    /* Set to zero the field value */
    for(int k=0 ; k<NumDOF ; k++){
      /* Update the field value */
      Phi_n.nM[k][i] = 0;
    }    
    /* Iterate over the nodes of the element */
    for(int j=0 ; j<ElementMesh.NumNodesElem ; j++){
      /* Iterate over the fields */
      for(int k=0 ; k<NumDOF ; k++){
  	/* Update the field value */
  	Phi_n.nM[k][i] += N_ref_XG.nV[j]*Phi_n_Nod.nM[k][Nodes_Elem_GP_i[j]-1];
      }
    }
  }

  /* Free memory */
  free(Phi_n_Nod.nM);
  
} /* End of Two_Steps_TG_Mc() */
