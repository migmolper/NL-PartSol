#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ToolsLib/TypeDefinitions.h"
#include "../MathTools/MathTools.h"
#include "MeshTools.h"


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
		  double Beta, double TOL_lambda)
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
  -> TOL_lambda : Tolerance for the lambda calculations.
*/
{
  /* Definition of some parameters */
  Matrix pa; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Increment_lambda;
  double norm_r = 10; /* Initial value of the norm */
  int NumIter = 0; /* Iterator counter */

  /* Start with the Newton-Rapson method */
  while(norm_r > TOL_lambda){
  
    /* Get vector with the shape functions evaluated in the nodes */
    pa = LME_pa(da,lambda,Beta);

    /* Get the gradient of log(Z) */
    r = LME_r(da,pa);

    /* Get the norm of r for the stopping criteria porpouse */
    norm_r = Norm_Mat(r,2);

    /* Get the Hessian of log(Z) */    
    J = LME_J(da,pa,r);

    /* Check the conditioning number of the Hessian */
    if (fabs(Cond_Mat(J)) < 1e-8){
      printf(" %s : %s \n",
	     "Error in LME_lambda",
	     "The Hessian is near to singular matrix");      
      exit(0);
    }

    /* Free the distance matrix */
    FreeMat(da);
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

    /* Update the value of lambda with the use of the increment */
    lambda = Incr_Mat(lambda,Increment_lambda);

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
  -> da : Matrix with the distance to the neighborhood node ''a'' (dim x 1).
  -> lambda : Initial value of the lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
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
  da_i.N_rows = N_dim;
  da_i.N_cols = 1;
  double Z_a = 0;
  double Z_a_m1 = 0;

  /* Get Z and the numerator */
  for(int i = 0 ; i<N_neibourg ; i++){
    da_i.nV = da.nM[i];
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
  Matrix J;
  Matrix J_I;
  Matrix J_II;
  Matrix da_i;
  da_i.N_rows = N_dim;
  da_i.N_cols = 1;
  Matrix da_iT;
  da_i.N_rows = 1;
  da_i.N_cols = N_dim;
  Matrix r_T;
  r_T.N_rows = 1;
  r_T.N_cols = N_dim;
  r_T.nV = r.nV;
  Matrix da_da; /* Tensorial product of da */

  /* Get the first component of the Hessian (J_I) */
  J_I = MatAllocZ(N_dim,N_dim); 
  for(int i = 0 ; i<N_neibourg ; i++){
    /* Get the tensorial product for each neighborhood. */
    da_i.nV = da.nM[i];
    da_iT.nV = da.nM[i];
    da_da = Tensorial_prod(da_i,da_iT);
    /* Fill the first component of the Hessian (J_I) */
    for(int j = 0 ; j<N_dim ; j++){
      for(int k = 0 ; k<N_dim ; k++){
	J_I.nM[j][k] += pa.nV[i]*da_da.nM[j][k];
      }
    }
    FreeMat(da_da);
  }

  /* Get the second component of the Hessian (J_II) */
  J_II = Tensorial_prod(r,r_T);

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
  da_i.N_rows = N_dim;
  da_i.N_cols = 1;
  Matrix Jm1_da; /* Auxiliar vector */
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = LME_r(da,pa);
  J = LME_J(da,pa,r);

  /* Check the conditioning number of the Hessian */
  if (fabs(Cond_Mat(J)) < 1e-8){
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
    da_i.nV = da.nM[i];
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

ChainPtr LME_Tributary_Nodes(Matrix X_EC_GP, int Elem_GP,
			     Mesh FEM_Mesh, double h, double Tol0){

  ChainPtr Triburary_Nodes = NULL;
  ChainPtr ChainElements = NULL;
  int * Tributary_Elements;
  int * NodesElem;
  int Elem_i;
  int Num_Elem;

  /* Get the search radius */
  double Ra = h*sqrt(-log(Tol0));
  
  /* Check if I am in the central area */
  if ((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
      (fabs(X_EC_GP.nV[1]) <= Dist[1])){
    Triburary_Nodes = CopyChain(FEM_Mesh.Connectivity[Elem_GP]);
  }    
  /* Check if I am in the 1º Quadrant */
  else if((X_EC_GP.nV[0]>=0) &&
	  (X_EC_GP.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);
    if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */       
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      
    }
    else if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[2]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if I am in the 2º Quadrant */
  else if((X_EC_GP.nV[0]<=0) &&
	  (X_EC_GP.nV[1]>=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) >= Dist[1])){

      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[3]]);
      Num_Elem = LenghtChain(ChainElements);      
      Tributary_Elements =  ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[3]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements =  ChainToArray(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);      
    }

    /* Free memory */
    free(NodesElem);
    
  }  
  /* Check if I am in the 3º Quadrant */
  else if((X_EC_GP.nV[0]<=0) &&
	  (X_EC_GP.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);   
    if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[3]],
					FEM_Mesh.NodeNeighbour[NodesElem[0]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);
      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */      
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[1]) >= Dist[1]) &&
	    (fabs(X_EC_GP.nV[0]) >= Dist[0])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[0]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);      
    }
    /* Free memory */
    free(NodesElem);    
  }  
  /* Check if it I am the 4º Quadrant */
  else if((X_EC_GP.nV[0]>=0) &&
	  (X_EC_GP.nV[1]<=0)){
    /* Create an array with the nodes of the element */
    NodesElem = ChainToArray(FEM_Mesh.Connectivity[Elem_GP],4);

    if((fabs(X_EC_GP.nV[0]) <= Dist[0]) &&
       (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[0]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements =  ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */      
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) <= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */ 
      ChainElements = ChainIntersection(FEM_Mesh.NodeNeighbour[NodesElem[2]],
					FEM_Mesh.NodeNeighbour[NodesElem[1]]);
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      FreeChain(ChainElements);

      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */      
      free(Tributary_Elements);
    }
    else if((fabs(X_EC_GP.nV[0]) >= Dist[0]) &&
	    (fabs(X_EC_GP.nV[1]) >= Dist[1])){
      /* Generate the list of Elements whose nodes contributes to the GP */
      ChainElements = FEM_Mesh.NodeNeighbour[NodesElem[1]];
      Num_Elem = LenghtChain(ChainElements);
      Tributary_Elements = ChainToArray(ChainElements,Num_Elem);
      /* Iterate in the list and select the union of the sets of nodes */
      for(int i = 0 ; i<Num_Elem ; i++){
	Elem_i = Tributary_Elements[i];
	Triburary_Nodes =
	  ChainUnion(Triburary_Nodes,FEM_Mesh.Connectivity[Elem_i]);
      }
      /* Free memory */
      free(Tributary_Elements);      
    }
    /* Free memory */
    free(NodesElem);
  }
  else{
    printf("%s : %s \n",
	   "Error in Tributary_Nodes_GIMP",
	   "Unlocated GP in the element");
    printf("%s : (%f;%f) \n",
	   "Natural coordinates of the GP",
	   X_EC_GP.nV[0],X_EC_GP.nV[1]);
    exit(0);
  }
 
  return Triburary_Nodes;
}
