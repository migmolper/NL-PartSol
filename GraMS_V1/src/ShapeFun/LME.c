#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/grams.h"

/**************************************************/
/************* Local Maximum-Entropy **************/
/**************************************************/

/*
  Shape functions based in :
  [1] : "" Local maximum-entropy approximation schemes : a seamless 
  bridge between finite elements and meshfree methods ""
  by M.Arroyo and M.Ortiz, 2006

  The employed nomenclature is the same 

*/

Matrix LME_lambda_NR(Matrix da, Matrix lambda, double Beta)
/*
  Output: 
  -> lambda : Lagrange multipliers lambda (1 x dim).
  Input parameters :
  -> da : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> lambda : Initial value of the
  lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
  -> h : Grid spacing (scalar).
  -> TOL_zero : Tolerance for Newton-Rapson.
*/
{  
  /* Definition of some parameters */
  Matrix pa; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Increment_lambda;
  double norm_r = 10; /* Value of the norm */
  int NumIter = 0; /* Iterator counter */

  /* Start with the Newton-Rapson method */
  while(norm_r > TOL_NR){
    
    /* Get vector with the shape functions evaluated in the nodes */
    pa = LME_pa(da,lambda,Beta);
    /* Get the gradient of log(Z) */
    r = LME_r(da,pa);
    /* Get the norm of r for the stopping criteria porpouse */
    norm_r = Norm_Mat(r,2);
    /* Get the Hessian of log(Z) */    
    J = LME_J(da,pa,r);

    /* Check the conditioning number of the Hessian */
    if (fabs(Cond_Mat(J)) > 10){
      printf(" %s : %s \n",
    	     "Error in LME_lambda",
    	     "The Hessian is near to singular matrix");
      exit(0);
    }

    /* Free the shape function nodal values */
    FreeMat(pa);
    
    /* Inverse of the Hessian */
    Jm1 = Get_Inverse(J);
    
    /* Get the increment for lambda */
    Increment_lambda = Scalar_prod(Jm1,r);

    /* Free r, J, and the inverse of J */
    FreeMat(r);
    FreeMat(J);
    FreeMat(Jm1);

    /* Update the value of lambda */
    for(int i = 0 ; i<NumberDimensions ; i++){
      lambda.nV[i] -= Increment_lambda.nV[i];
    }

    /* Free memory */
    FreeMat(Increment_lambda);

    /* Update the number of iterations */
    NumIter ++;
    if(NumIter > 100){
      printf(" %s : %s \n",
	     "Error in LME_lambda",
	     "No convergence in 100 iterations");
    }
  
  }
  
  /* Once the stopping criteria is reached, 
     return the lagrange multipliers value */
  return lambda;
}

double LME_fa(Matrix da, Matrix lambda, double Beta)
/*
  Output :
  -> fa : the function fa that appear in [1] (scalar).
  Input parameters :
  -> da : Matrix with the distance to the neighborhood node ''a'' (1 x dim).
  -> lambda : Initial value of the lagrange multipliers (dim x 1).
  -> Gamma : Tunning parameter (scalar).
*/
{  
  /* Get the scalar product the distance and the lagrange multipliers */
  Matrix Aux = Scalar_prod(lambda,da);
  double norm_dist = Norm_Mat(da,2);

  /* Return the value of f_a*/
  return -Beta*norm_dist*norm_dist + Aux.n;
}

Matrix LME_pa(Matrix da, Matrix lambda, double Beta)
/*
  Output :
  -> pa : Value of the shape function in the
  neibourhud nodes (1 x neighborhood).
  Input parameters :
  -> da : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> lambda : Initial value of the lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
*/
{
  
  /* Definition of some parameters */
  int N_dim = da.N_cols;
  int N_neibourg = da.N_rows;
  Matrix pa = MatAlloc(1,N_neibourg);
  Matrix da_i;
  double Z_a = 0;
  double Z_a_m1 = 0;

  /* Get Z and the numerator */
  for(int i = 0 ; i<N_neibourg ; i++){
    da_i = MatAssign(1,N_dim,NAN,da.nM[i],NULL);
    pa.nV[i] = exp(LME_fa(da_i,lambda,Beta));
    Z_a += pa.nV[i];
  }

  /* Get the inverse of Z */
  Z_a_m1 = (double)1/Z_a;

  /* Divide by Z and get the final value */
  for(int i = 0 ; i<N_neibourg ; i++){
    pa.nV[i] *= Z_a_m1;
  }
  
  /* Return the value of the shape function */  
  return pa;
}


Matrix LME_r(Matrix da, Matrix pa)
/*
  Output :
  -> r : Gradient of the log(Z) function (dim x 1).
  Input parameters :
  -> da : Matrix with the distances to the 
  neighborhood nodes (neighborhood x dim).
  -> pa : Shape function value in the
  neighborhood nodes (1 x neighborhood).
*/
{  
  /* Definition of some parameters */
  int N_dim = da.N_cols;
  int N_neibourg = da.N_rows;
  Matrix r = MatAllocZ(N_dim,1);

  /* Fill ''r'' */
  for(int i = 0 ; i<N_neibourg ; i++){
    for(int j = 0 ; j<N_dim ; j++){
      r.nV[j] += pa.nV[i]*da.nM[i][j];
    }
  }

  /* Return the value of the gradient */
  return r;
}

Matrix LME_J(Matrix da, Matrix pa, Matrix r)
/*
  Output :
  -> J : Hessian of the log(Z) function (dim x dim).
  Input parameters :
  -> da : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> pa : Shape function value in the
  neighborhood nodes (neighborhood x 1).
  -> r : Gradient of log(Z) (dim x 1).
*/
{  
  /* Definition of some parameters */
  int N_neibourg = da.N_rows;
  int N_dim = da.N_cols;
  Matrix J; /* Hessian definition */
  Matrix J_I; /* First component of the Hessian */
  Matrix J_II; /* Second component of the Hessian */
  Matrix da_da; /* Tensorial product of da */

  /* Get the first component of the Hessian (J_I) */
  J_I = MatAllocZ(N_dim,N_dim); 
  for(int i = 0 ; i<N_neibourg ; i++){
    /* Get the tensorial product for each neighborhood. */
    da_da = Tensorial_prod(MatAssign(N_dim,1,NAN,da.nM[i],NULL),
			   MatAssign(1,N_dim,NAN,da.nM[i],NULL));
    /* Fill the first component of the Hessian (J_I) */
    for(int j = 0 ; j<N_dim ; j++){
      for(int k = 0 ; k<N_dim ; k++){
	J_I.nM[j][k] += pa.nV[i]*da_da.nM[j][k];
      }
    }
    /* Free tensorial product */
    FreeMat(da_da);
  }

  /* Get the second component of the Hessian (J_II) */
  J_II = Tensorial_prod(r,MatAssign(1,N_dim,NAN,r.nV,NULL));

  /* Get the Hessian */
  J = Sub_Mat(J_I,J_II);

  /* Free the auxiliar components of the Hessian */
  FreeMat(J_I);
  FreeMat(J_II);

  /* Return the value of the Hessian */
  return J;
}

Matrix LME_dpa(Matrix da, Matrix pa)
/*
  Output :
  -> dpa : Value of the shape function gradient in 
  the neighborhood nodes (dim x neighborhood).
  Input parameters :
  -> da : Matrix with the distances to the
  neighborhood nodes (neighborhood x dim).
  -> pa : Shape function value in the
  neighborhood nodes (neighborhood x 1).
*/
{  
  /* Definition of some parameters */
  int N_neibourg = da.N_rows;
  int N_dim = da.N_cols;
  Matrix dpa = MatAllocZ(N_dim,N_neibourg);
  Matrix r;
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix da_i; /* Distance to the node */
  Matrix Jm1_da; /* Auxiliar vector */
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = LME_r(da,pa);
  J = LME_J(da,pa,r);

  /* Check the conditioning number of the Hessian */
  if (fabs(Cond_Mat(J)) > 10){
    printf(" %s : %s \n",
  	   "Error in LME_lambda",
  	   "The Hessian is near to singular matrix");
    exit(0);
  }
    
  /* Inverse of the Hessian */
  Jm1 = Get_Inverse(J);

  /* Free memory */
  FreeMat(r);
  FreeMat(J);

  /* Fill the gradient for each node */
  for(int i = 0 ; i<N_neibourg ; i++){
    da_i = MatAssign(N_dim,1,NAN,da.nM[i],NULL); 
    Jm1_da = Scalar_prod(Jm1,da_i);    
    for(int j = 0 ; j<N_dim ; j++){
      dpa.nM[j][i] = pa.nV[i]*Jm1_da.nV[j];
    }
    FreeMat(Jm1_da);
  }

  /* Free memory */
  FreeMat(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dpa;
  
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

    /* If the node is not near the GP pop out of the chain */
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

