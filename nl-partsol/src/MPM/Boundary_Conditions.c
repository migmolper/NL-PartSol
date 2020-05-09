#include "grams.h"

/**********************************************************************/

void imposse_NodalMomentum(Mesh FEM_Mesh, Matrix Phi_I, int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for(int j = 0 ; j<NumNodesBound ; j++){
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for(int k = 0 ; k<NumDimBound ; k++){
	/* 7º Apply only if the direction is active (1) */
	if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1){
	  /* 8º Check if the curve it is on time */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num)){
	    printf("%s : %s \n",
		   "Error in imposse_NodalMomentum()",
		   "The time step is out of the curve !!");
	    exit(0);
	  }
	  /* 9º Assign the boundary condition */
	  Phi_I.nM[Id_BCC][k] =
	    FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
	}
      }
    }    
  }

}

/**********************************************************************/

void imposse_NodalVelocity(Mesh FEM_Mesh, Matrix V_I, int TimeStep)
/*
  Apply the boundary conditions over the nodes 
*/
{

  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for(int j = 0 ; j<NumNodesBound ; j++){
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for(int k = 0 ; k<NumDimBound ; k++){
	/* 7º Apply only if the direction is active (1) */
	if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1){
	  /* 8º Check if the curve it is on time */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num)){
	    printf("%s : %s \n",
		   "Error in imposse_NodalMomentum()",
		   "The time step is out of the curve !!");
	    exit(0);
	  }
	  /* 9º Assign the boundary condition */
	  V_I.nM[Id_BCC][k] =
	    FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
	}
      }
    }    
  }

}

/*********************************************************************/

Matrix Eval_Body_Forces(Load * B, int NumLoads, int NumGP, int TimeStep)
/*
  Evaluate body forces in a time step.
*/
{
  Matrix Body_Forces_t =
    MatAllocZ(NumberDimensions,NumGP);
  int GP_Force;
  
  if(NumLoads>0){
    for(int i = 0 ; i<NumLoads; i++){
      for(int j = 0 ; j<B[i].NumNodes ; j++){
	GP_Force = B[i].Nodes[j];
	for(int k = 0 ; k<B[i].Dim ; k++){
	  if( (B[i].Dir[k] == 1) ||
	      (B[i].Dir[k] == -1)){
	    if( (TimeStep < 0) || (TimeStep > B[i].Value[k].Num)){
	      printf("%s : %s\n",
		     "Error in Eval_Body_Forces()",
		     "The time step is out of the curve !!");
	      exit(0);
	    }
	    Body_Forces_t.nM[k][GP_Force] +=
	      B[i].Value[k].Fx[TimeStep]*
	      (double)B[i].Dir[k];
	  }
	}
      }
    }
  }

  return Body_Forces_t;

}

/*********************************************************************/

Matrix Eval_Contact_Forces(Load * F, int NumLoads, int NumGP, int TimeStep)
/*
  Evaluate contact forces in a time step
 */
{
  int GP_Force;
  Matrix Contact_Forces_t = 
    MatAllocZ(NumberDimensions,NumGP);
  
  if(NumLoads>0){
    for(int i = 0 ; i<NumLoads; i++){
      for(int j = 0 ; j<F[i].NumNodes ; j++){
	GP_Force = F[i].Nodes[j];
	for(int k = 0 ; k<F[i].Dim ; k++){
	  if( (F[i].Dir[k] == 1) ||
	      (F[i].Dir[k] == -1)){
	    if( (TimeStep < 0) || (TimeStep > F[i].Value[k].Num)){
	      printf("%s : %s \n",
		     "Error in Eval_Contact_Forces()",
		     "The time step is out of the curve !!");
	      exit(0);
	    }
	    Contact_Forces_t.nM[k][GP_Force] +=
	      F[i].Value[k].Fx[TimeStep]*
	      (double)F[i].Dir[k];
	  }
	}
      }
    }
  }

  return Contact_Forces_t;
}

/**********************************************************************/

Matrix compute_Reactions(Mesh FEM_Mesh, Matrix F_I)
/*
  Compute the nodal reactions
*/
{
  /* 1º Define auxilar variables */
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound; /* Number of dimensions */
  int Id_BCC; /* Index of the node where we apply the BCC */
  int Ndim = NumberDimensions;

  Matrix R_I = MatAllocZ(FEM_Mesh.NumNodesMesh,Ndim);
  strcpy(R_I.Info,"REACTIONS");

  /* 2º Loop over the the boundaries */
  for(int i = 0 ; i<FEM_Mesh.Bounds.NumBounds ; i++){
    /* 3º Get the number of nodes of this boundarie */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;
    /* 4º Get the number of dimensions where the BCC it is applied */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;
    for(int j = 0 ; j<NumNodesBound ; j++){
      /* 5º Get the index of the node */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      /* 6º Loop over the dimensions of the boundary condition */
      for(int k = 0 ; k<NumDimBound ; k++){
	/* 7º Apply only if the direction is active (1) */
	if(FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1){
	  /* 8º Set to zero the forces in the nodes where velocity is fixed */
	  R_I.nM[Id_BCC][k] = F_I.nM[Id_BCC][k];
	  F_I.nM[Id_BCC][k] = 0;
	}
      }
    }    
  }

  return R_I;
}

/**********************************************************************/


void free_Load(Load F){

  free(F.Nodes);

}

/**********************************************************************/

void free_LoadTable(Load * Table){

}

/**********************************************************************/


void free_Boundaries(Boundaries Bounds){

}

/**********************************************************************/