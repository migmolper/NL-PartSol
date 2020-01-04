#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/grams.h"

/**************************************************/
/************* Local Maximum-Entropy **************/
/**************************************************/

/*
  Shape functions based in :
  "" Local maximum-entropy approximation schemes : a seamless 
  bridge between finite elements and meshfree methods ""
  by M.Arroyo and M.Ortiz, 2006.

  Here the employed nomenclature for the code is the same as in the paper.
*/

/****************************************************************************/

Matrix Initialize_lambda(Matrix d, Matrix lambda, double Beta)
{
  int Ndim = NumberDimensions;;
  Matrix A = MatAlloc(Ndim,Ndim);
  Matrix B = MatAlloc(Ndim,1);

  /* Fill A */
  A.nM[0][0] = d.nM[1][0] - d.nM[0][0]; 
  A.nM[0][1] = d.nM[1][1] - d.nM[0][1]; 
  A.nM[1][0] = d.nM[2][0] - d.nM[0][0]; 
  A.nM[1][1] = d.nM[2][1] - d.nM[0][1];
    
  /* Fill B */
  B.nV[1] = -Beta*(d.nM[1][1]*d.nM[1][1] + d.nM[1][2]*d.nM[1][2] -
		   (d.nM[2][1]*d.nM[2][1] + d.nM[2][2]*d.nM[2][2]));
  B.nV[2] = -Beta*(d.nM[1][1]*d.nM[1][1] + d.nM[1][2]*d.nM[1][2] -
		   (d.nM[3][1]*d.nM[3][1] + d.nM[3][2]*d.nM[3][2]));

  /* Update the value of lambda */
  lambda = Conjugate_Gradient_Method(A,B,lambda);

  /* Free memory */
  FreeMat(A), FreeMat(B);

  /* Return lambda */
  return lambda;
}

/****************************************************************************/

Matrix LME_lambda_NR(Matrix l, Matrix lambda, double Beta)
/*
  Get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method for that is the Newton-Rapson.

  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> lambda : Initial value of the
  lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
  -> h : Grid spacing (scalar).
  -> TOL_zero : Tolerance for Newton-Rapson.
*/
{
  /* Definition of some parameters */
  int MaxIter = 300;
  int Ndim = NumberDimensions;
  int NumIter = 0; /* Iterator counter */
  double norm_r = 10; /* Value of the norm */
  Matrix p; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix D_lambda; /* Increment of lambda */
 
  /* Start with the Newton-Rapson method */
  while(norm_r > TOL_NR){
	
    /* Get vector with the shape functions evaluated in the nodes */
    p = LME_p(l,lambda,Beta);

    /* Get the gradient of log(Z) and its norm */
    r = LME_r(l,p), norm_r = Norm_Mat(r,2);

    /* Get the Hessian of log(Z) */    
    J = LME_J(l,p,r);

    /* /\* Check the conditioning number of the Hessian *\/ */
    /* if (fabs(Cond_Mat(J)) > 10){ */
    /*   printf(" %s : %s \n", */
    /* 	     "Error in LME_lambda_NR", */
    /* 	     "The Hessian is near to singular matrix"); */
    /*   exit(0); */
    /* } */
    
    /* Get the increment of lambda */
    D_lambda = Solve_Linear_Sistem(J,r,MatAllocZ(Ndim,1));

    /* Update the value of lambda */
    for(int i = 0 ; i<NumberDimensions ; i++){
      if(fabs(D_lambda.nV[i])<TOL_NR) continue;
      lambda.nV[i] -= D_lambda.nV[i];
    }

    /* Free memory */
    FreeMat(p), FreeMat(r), FreeMat(J), FreeMat(D_lambda);
    
    /* Update the number of iterations */
    NumIter ++;
    if(NumIter >= MaxIter){
      printf("%s : %s %i/%i %s \n",
	     "Error in LME_lambda",
	     "No convergence in",NumIter,MaxIter,"iterations");
      break;
    }
  
  }
  if(NumIter >= 300){
    printf("%s %i %s : %f \n",
	   "Error after",NumIter,
	   "iterations",norm_r);
  }
  
  /* Once the stopping criteria is reached, 
     return the lagrange multipliers value */
  return lambda;
}

/****************************************************************************/

double LME_fa(Matrix la, Matrix lambda, double Beta)
/*
  Output :
  -> fa : the function fa that appear in [1] (scalar).
  Input parameters :
  -> la : Matrix with the distance to the neighborhood node ''a'' (1 x dim).
  -> lambda : Initial value of the lagrange multipliers (dim x 1).
  -> Gamma : Tunning parameter (scalar).
*/
{  
  /* Get the scalar product the distance and the lagrange multipliers */
  Matrix Aux = Scalar_prod(lambda,la);
  double norm_dist = Norm_Mat(la,2);

  /* Return the value of fa */
  return -Beta*norm_dist*norm_dist + Aux.n;
}

Matrix LME_p(Matrix l, Matrix lambda, double Beta)
/*
  Get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.

  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neiborghood x dim).
  -> lambda : Initial value of the lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
*/
{
  
  /* Definition of some parameters */
  int N_a = l.N_rows, N_dim = NumberDimensions;
  double Z = 0, Z_m1 = 0;
  Matrix p = /* Vector with the values of the shape-function in the nodes */
    MatAlloc(1,N_a); 
  Matrix la = /* Distance to the neighbour (x-x_a) */
    MatAssign(1,N_dim,NAN,NULL,NULL);

  /* Get Z and the numerator */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a];
    p.nV[a] = exp(LME_fa(la,lambda,Beta));
    Z += p.nV[a];
  }

  /* Get the inverse of Z */
  Z_m1 = (double)1/Z;

  /* Divide by Z and get the final value */
  for(int a = 0 ; a<N_a ; a++){
    p.nV[a] *= Z_m1;
  }
  
  /* Return the value of the shape function */  
  return p;
}

Matrix LME_r(Matrix l, Matrix p)
/*
  Get the gradient "r" (dim x 1) of the function log(Z) = 0.

  Input parameters :
  -> l : Matrix with the distances to the 
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (1 x neighborhood).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows, N_dim = NumberDimensions;
  Matrix r /* Value of the gradient */
    = MatAllocZ(NumberDimensions,1);

  /* Fill ''r'' */
  for(int i = 0 ; i<N_dim ; i++){
    for(int a = 0 ; a<N_a ; a++){
      r.nV[i] += p.nV[a]*l.nM[a][i];
    }
  }

  /* Return the value of the gradient */
  return r;
}

Matrix LME_J(Matrix l, Matrix p, Matrix r)
/*
  Get the Hessian "J" (dim x dim) of the function log(Z) = 0.

  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (neighborhood x 1).
  -> r : Gradient of log(Z) (dim x 1).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int N_dim = NumberDimensions;
  Matrix J; /* Hessian definition */
  
  /* Allocate Hessian */
  J = MatAllocZ(N_dim,N_dim);

  /* Fill the Hessian */
  for(int i = 0 ; i<N_dim ; i++){
    for(int j = 0 ; j<N_dim ; j++){
      for(int a = 0 ; a<N_a ; a++){
	/* Get the first component of the Hessian looping
	   over the neighborhood nodes. */
	J.nM[i][j] += p.nV[a]*l.nM[a][i]*l.nM[a][j];
      }
      /* Get the second value of the Hessian */
      J.nM[i][j] -= r.nV[i]*r.nV[j];
    }
  }

  /* Return the value of the Hessian */
  return J;
}

Matrix LME_dp(Matrix l, Matrix p)
/*
  Value of the shape function gradient "dp" (dim x neighborhood) in 
  the neighborhood nodes.

  Input parameters :
  -> l : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (neighborhood x 1).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows, N_dim = NumberDimensions;
  Matrix dp = MatAllocZ(N_dim,N_a);
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Jm1_la; /* Auxiliar vector */
  Matrix la = /* Distance to the neighbour (x-x_a) */
    MatAssign(N_dim,1,NAN,NULL,NULL); 
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = LME_r(l,p);
  J = LME_J(l,p,r);
  
  /* Inverse of the Hessian */
  Jm1 = Get_Inverse(J);
  
  /* Free memory */
  FreeMat(r), FreeMat(J);

  /* Fill the gradient for each node */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a]; 
    Jm1_la = Scalar_prod(Jm1,la);    
    for(int i = 0 ; i<N_dim ; i++){
      dp.nM[i][a] = - p.nV[a]*Jm1_la.nV[i];
    }
    FreeMat(Jm1_la);
  }

  /* Free memory */
  FreeMat(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dp;  
}

ChainPtr LME_Tributary_Nodes(Matrix X_GP, int Elem_GP,
			     Mesh FEM_Mesh, double Gamma){

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  ChainPtr * Table_Elem = NULL;
  ChainPtr Triburary_Nodes = NULL;
  ChainPtr List_Nodes = NULL;
  ChainPtr * Table_ElemNodes = NULL;
  ChainPtr Triburary_Elements = NULL;
  ChainPtr iPtr = NULL;
  int * List_Elements;
  int * NodesElem;
  int Num_Elem;
  int NumNodesElem;
  double Ra;

  /* Get the search radius */
  Ra = FEM_Mesh.DeltaX*sqrt(-log(TOL_lambda)/Gamma);

  /* Number of nodes of the initial element and
     list of nodes */
  NumNodesElem = FEM_Mesh.NumNodesElem[Elem_GP];
  NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],
			   NumNodesElem);

  /* Chain with the tributary elements, this is the list of element near the
     gauss point, including where it is */

  /* Iterate in the list and select the union of the sets of nodes */
  Table_Elem = malloc(NumNodesElem*sizeof(ChainPtr));
  for(int i = 0 ; i<NumNodesElem ; i++){
    Table_Elem[i] = FEM_Mesh.NodeNeighbour[NodesElem[i]];
  }
  Triburary_Elements = ChainUnion(Table_Elem,NumNodesElem);
  /* Free memory */
  free(NodesElem);
  free(Table_Elem);
  Table_Elem = NULL;
  
  /* List with the tributary nodes */
  Num_Elem = LenghtChain(Triburary_Elements);
  List_Elements = ChainToArray(Triburary_Elements,Num_Elem);

  /* Free the chain with the tributary elements */
  FreeChain(Triburary_Elements);
  
  /* Fill the chain with the preliminary tributary nodes */
  Table_ElemNodes = malloc(Num_Elem*sizeof(ChainPtr));
  for(int i = 0 ; i<Num_Elem ; i++){
    Table_ElemNodes[i] = FEM_Mesh.Connectivity[List_Elements[i]];
  }

  List_Nodes = ChainUnion(Table_ElemNodes,Num_Elem);
  
  /* Free the array wit the list of tributary elements */
  free(List_Elements);
  free(Table_ElemNodes);
  Table_ElemNodes = NULL;
  
  /* Initialize the iterator to iterate over the list of tributary nodes */
  iPtr = List_Nodes;

  /* Loop over the chain with the tributary nodes */
  while(iPtr != NULL){

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[iPtr->I];

    /* Get a vector from the GP to the node */
    Distance = Sub_Mat(X_GP,X_I);

    /* If the node is near the GP push in the chain */
    if(Norm_Mat(Distance,2) <= Ra){
      PushNodeTop(&Triburary_Nodes,iPtr->I);
    }

    /* Free memory of the distrance vector */
    FreeMat(Distance);

    /* Update pointer index */
    iPtr = iPtr->next;
  }
  /* Free memory */
  FreeChain(List_Nodes);

  return Triburary_Nodes;
}

