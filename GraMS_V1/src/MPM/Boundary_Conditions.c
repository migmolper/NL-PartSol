#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../GRAMS/grams.h"

/**********************************************************************/

Curve BcDirichlet(char * Expresion)
/*
 BcDirichlet V 0 
*/
{
  Curve BcDir;

  BcDir = ReadCurve(Expresion);
  
  return BcDir;
}

/**********************************************************************/

void BCC_Nod_VALUE(Mesh FEM_Mesh, Matrix Nodal_VALUE, int TimeStep)
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
	if( (FEM_Mesh.Bounds.BCC_i[i].Dir[k] == 1) ||
	    (FEM_Mesh.Bounds.BCC_i[i].Dir[k] == -1)){
	  /* 8º Check if the curve it is on time */
	  if( (TimeStep < 0) ||
	      (TimeStep > FEM_Mesh.Bounds.BCC_i[i].Value[k].Num)){
	    printf("%s : %s \n",
		   "Error in BCC_Nodal_VALUE()",
		   "The time step is out of the curve !!");
	    exit(0);
	  }
	  /* 9º Assign the boundary condition */
	  Nodal_VALUE.nM[k][Id_BCC] =
	    FEM_Mesh.Bounds.BCC_i[i].Value[k].Fx[TimeStep]*
	    (double)FEM_Mesh.Bounds.BCC_i[i].Dir[k];
	}
      }
    }    
  } 

}

/*********************************************************************/

Matrix Eval_Body_Forces(LoadCase B, int NumGP, int TimeStep)
/*
  Evaluate body forces in a time step.
*/
{
  Matrix Body_Forces_t =
    MatAllocZ(NumberDimensions,NumGP);
  int GP_Force;
  
  if(B.NumLoads>0){
    for(int i = 0 ; i<B.NumLoads; i++){
      for(int j = 0 ; j<B.Load_i[i].NumNodes ; j++){
	GP_Force = B.Load_i[i].Nodes[j];
	for(int k = 0 ; k<B.Load_i[i].Dim ; k++){
	  if( (B.Load_i[i].Dir[k] == 1) ||
	      (B.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > B.Load_i[i].Value[k].Num)){
	      printf("%s : %s\n",
		     "Error in GetNodalForces()",
		     "The time step is out of the curve !!");
	      exit(0);
	    }
	    Body_Forces_t.nM[k][GP_Force] +=
	      B.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)B.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  return Body_Forces_t;

}


/*********************************************************************/

Matrix Eval_Contact_Forces(LoadCase F, int NumGP, int TimeStep)
/*
  Evaluate contact forces in a time step
 */
{
  int GP_Force;
  Matrix Contact_Forces_t = 
    MatAllocZ(NumberDimensions,NumGP);
  
  if(F.NumLoads>0){
    for(int i = 0 ; i<F.NumLoads; i++){
      for(int j = 0 ; j<F.Load_i[i].NumNodes ; j++){
	GP_Force = F.Load_i[i].Nodes[j];
	for(int k = 0 ; k<F.Load_i[i].Dim ; k++){
	  if( (F.Load_i[i].Dir[k] == 1) ||
	      (F.Load_i[i].Dir[k] == -1)){
	    if( (TimeStep < 0) ||
		(TimeStep > F.Load_i[i].Value[k].Num)){
	      printf("%s : %s \n",
		     "Error in GetNodalForces()",
		     "The time step is out of the curve !!");
	      exit(0);
	    }
	    Contact_Forces_t.nM[k][GP_Force] +=
	      F.Load_i[i].Value[k].Fx[TimeStep]*
	      (double)F.Load_i[i].Dir[k];
	  }
	}
      }
    }
  }

  return Contact_Forces_t;
}

/*********************************************************************/
