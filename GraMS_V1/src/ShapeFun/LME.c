#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../GRAMS/grams.h"
#include "../GRAMS/Utils.h"

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

Matrix LME_lambda(Matrix da, Matrix lambda,
		  double DeltaX, double Gamma)
/*
  Output: 
  -> lambda : Lagrange multipliers lambda for
  a material point (1 x dim).
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
  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();
  
  /* Definition of some parameters */
  Matrix pa; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Increment_lambda;
  double norm_r = 10; /* Initial value of the norm */
  int NumIter = 0; /* Iterator counter */

  /* Start with the Newton-Rapson method */
  while(norm_r > TOL_zero){
  
    /* Get vector with the shape functions evaluated in the nodes */
    pa = LME_pa(da,lambda,DeltaX,Gamma);

    /* Get the gradient of log(Z) */
    r = LME_r(da,pa);

    /* Get the norm of r for the stopping criteria porpouse */
    norm_r = MO.Norm(r,2);

    /* Get the Hessian of log(Z) */    
    J = LME_J(da,pa,r);

    /* Check the conditioning number of the Hessian */
    if (fabs(MO.Cond(J)) < 1e-8){
      printf(" %s : %s \n",
	     "Error in LME_lambda",
	     "The Hessian is near to singular matrix");      
      exit(0);
    }

    /* Free the distance matrix */
    MO.FreeMat(da);
    /* Free the shape function nodal values */
    MO.FreeMat(pa);
    
    /* Inverse of the Hessian */
    Jm1 = MO.Inv(J);

    /* Get the increment for lambda */
    Increment_lambda = MO.Sprod(Jm1,r);

    /* Free r, J, and the inverse of J */
    MO.FreeMat(r);
    MO.FreeMat(J);
    MO.FreeMat(Jm1);   

    /* Update the value of lambda with the use of the increment */
    lambda = MO.Incr(lambda,Increment_lambda);

    /* Free memory */
    MO.FreeMat(Increment_lambda);

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

double LME_fa(Matrix da, Matrix lambda, double DeltaX, double Gamma)
/*
  Output :
  -> fa : the function fa that appear in [1] (scalar).
  Input parameters :
  -> da : Matrix with the distance to the neighborhood node ''a'' (dim x 1).
  -> lambda : Initial value of the lagrange multipliers (1 x dim).
  -> Gamma : Tunning parameter (scalar).
*/
{
  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();

  /* Get Beta */
  double Beta = Gamma/(DeltaX*DeltaX);
  
  /* Get the scalar product the distance and the lagrange multipliers */
  Matrix Aux = MO.Sprod(lambda,da);
  double norm_dist = MO.Norm(da,2);

  /* Return the value of f_a*/
  return -Beta*norm_dist*norm_dist + Aux.n;
}

Matrix LME_pa(Matrix da, Matrix lambda, double DeltaX, double Gamma)
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
  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();
  
  /* Definition of some parameters */
  int N_dim = da.N_cols;
  int N_neibourg = da.N_rows;
  Matrix pa = MO.Alloc(1,N_neibourg);
  Matrix da_i;
  double Z_a = 0;
  double Z_a_m1 = 0;

  /* Get Z and the numerator */
  for(int i = 0 ; i<N_neibourg ; i++){
    da_i = MO.Assign(N_dim,1,NAN,da.nM[i],NULL);
    pa.nV[i] = exp(LME_fa(da_i,lambda,DeltaX,Gamma));
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
  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();
  
  /* Definition of some parameters */
  int N_dim = da.N_cols;
  int N_neibourg = da.N_rows;
  Matrix r = MO.AllocZ(N_dim,1);

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
  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();
  
  /* Definition of some parameters */
  int N_neibourg = da.N_rows;
  int N_dim = da.N_cols;
  Matrix J;
  Matrix J_I;
  Matrix J_II;
  Matrix da_i;
  Matrix da_iT;
  Matrix r_T;
  Matrix da_da; /* Tensorial product of da */

  /* Get the first component of the Hessian (J_I) */
  J_I = MO.AllocZ(N_dim,N_dim); 
  for(int i = 0 ; i<N_neibourg ; i++){
    /* Get the tensorial product for each neighborhood. */
    da_i = MO.Assign(N_dim,1,NAN,da.nM[i],NULL);
    da_iT = MO.Assign(1,N_dim,NAN,da.nM[i],NULL);
    da_da = MO.Tprod(da_i,da_iT);
    /* Fill the first component of the Hessian (J_I) */
    for(int j = 0 ; j<N_dim ; j++){
      for(int k = 0 ; k<N_dim ; k++){
	J_I.nM[j][k] += pa.nV[i]*da_da.nM[j][k];
      }
    }
    MO.FreeMat(da_da);
  }

  /* Get the second component of the Hessian (J_II) */
  r_T = MO.Assign(1,N_dim,NAN,r.nV,NULL);
  J_II = MO.Tprod(r,r_T);

  /* Get the Hessian */
  J = MO.Sub(J_I,J_II);

  /* Free the auxiliar components of the Hessian */
  MO.FreeMat(J_I);
  MO.FreeMat(J_II);

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

  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();
  
  /* Definition of some parameters */
  int N_neibourg = da.N_rows;
  int N_dim = da.N_cols;
  Matrix dpa = MO.AllocZ(N_dim,N_neibourg);
  Matrix r;
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix da_i; /* Distance to the node */
  Matrix Jm1_da; /* Auxiliar vector */
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = LME_r(da,pa);
  J = LME_J(da,pa,r);

  /* Check the conditioning number of the Hessian */
  if (fabs(MO.Cond(J)) < 1e-8){
    printf(" %s : %s \n",
	   "Error in LME_lambda",
	   "The Hessian is near to singular matrix");      
    exit(0);
  }
    
  /* Inverse of the Hessian */
  Jm1 = MO.Inv(J);

  /* Free memory */
  MO.FreeMat(r);
  MO.FreeMat(J);

  /* Fill the gradient for each node */
  for(int i = 0 ; i<N_neibourg ; i++){
    da_i = MO.Assign(N_dim,1,NAN,da.nM[i],NULL); 
    Jm1_da = MO.Sprod(Jm1,da_i);    
    for(int j = 0 ; j<N_dim ; j++){
      dpa.nM[j][i] = pa.nV[i]*Jm1_da.nV[j];
    }
    MO.FreeMat(Jm1_da);
  }

  /* Free memory */
  MO.FreeMat(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dpa;
  
}

ChainPtr LME_Tributary_Nodes(Matrix X_GP, int Elem_GP,
			     Mesh FEM_Mesh, double Gamma){

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = MatAssign(NumberDimensions,1,NAN,NULL,NULL);
  ChainPtr Triburary_Nodes = NULL;
  ChainPtr Triburary_Elements = NULL;
  ChainPtr iPtr = NULL;
  ChainPtr PrevPtr = NULL;
  ChainPtr AuxPtr;
  int * List_Elements;
  int * NodesElem;
  int Num_Elem;
  int NumNodesElem;
  double Ra;

  /* CALL LIBRARIES */
  MatLib MO = MatrixOperators();

  /* Get the search radius */
  Ra = FEM_Mesh.DeltaX*sqrt(-log(TOL_lambda)/Gamma);

  /* Number of nodes of the initial element and
     list of nodes */
  NumNodesElem = FEM_Mesh.NumNodesElem[Elem_GP];
  NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],
			   NumNodesElem);

  /* Chain with the tributary elements, this is the list of element near the
     gauss point, including where it is */
  for(int i = 0 ; i<NumNodesElem ; i++){
    Triburary_Elements =
      ChainUnion(Triburary_Elements,FEM_Mesh.NodeNeighbour[NodesElem[i]]);
  }

  /* Free the array with the nodes of the initial element */
  free(NodesElem);

  /* List with the tributary nodes */
  Num_Elem = LenghtChain(Triburary_Elements);
  List_Elements = ChainToArray(Triburary_Elements,Num_Elem);

  /* Free the chain wit the tributary elements */
  FreeChain(Triburary_Elements);
  
  /* Fill the chain with the preliminary tributary nodes */
  for(int i = 0 ; i<Num_Elem ; i++){
    Triburary_Nodes =
      ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[List_Elements[i]]);
  }

  /* Free the array wit the list of tributary elements */
  free(List_Elements);

  /* Initialize the iterator to iterate over the list of tributary nodes */
  iPtr = Triburary_Nodes;

  /* Loop over the chain with the tributary nodes */
  while(iPtr != NULL){

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[iPtr->I];

    /* Get a vector from the GP to the node */
    Distance = MO.Sub(X_GP,X_I);

    /* If the node is not near the GP pop out of the chain */
    if(MO.Norm(Distance,2) > Ra){
      /* If the node is the first in the chain */
      if(PrevPtr == NULL){
	AuxPtr = iPtr->next;
	free(iPtr);
	Triburary_Nodes->next = AuxPtr;
      }
      /* If the node is in the middle or at the end */
      else{
	PrevPtr->next = iPtr->next;
	free(iPtr);
      }
      /* Once the node is located, breack the loop */
      break;
    }

    /* Free memory of the distrance vector */
    MO.FreeMat(Distance);

    /* The previous is the index */
    PrevPtr = iPtr;
    /* Update pointer index */
    iPtr = iPtr->next;
  }

  return Triburary_Nodes;
}

