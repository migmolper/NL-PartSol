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
  "" Local maximum-entropy approximation schemes : a seamless 
  bridge between finite elements and meshfree methods ""
  by M.Arroyo and M.Ortiz, 2006

  The employed nomenclature is the same 

*/

Matrix LME_lambda(Matrix da, Matrix lambda,
		  double Beta, double h, double TOL_lambda)
/*
  Especialized Newton-Rapson algorithm to get the lagrange multipliers lambda
  for a material point.
  Input parameters :
  - da : Matrix with the distances to the neighborhood nodes (neighborhood x dim).
  - lambda : Initial value of the lagrange multipliers (1 x dim).
  - Beta : Tunning parameter (scalar).
  - h : Grid spacing (scalar).
  - TOL_lambda : Tolerance for the lambda calculations.
*/
{
  /* Definition of some parameters */
  Matrix pa; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
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

*/
{

  /* Get the scalar product the distance and the lagrange multipliers */
  Matrix Aux = Scalar_prod(lambda,da);

  /* Return the value of f_a*/
  return -Beta*Norm_Mat(da,2) + Aux.n;
}

Matrix LME_pa(Matrix da, Matrix lambda, double Beta)
/*
  Value of the shape function evaluated in the neibourhud nodes.
  Input parameters :
  - da : Matrix with the distances to the neighborhood nodes (neighborhood x dim).
  - lambda : Initial value of the lagrange multipliers (1 x dim).
  - Beta : Tunning parameter (scalar).
*/
{

  Matrix pa = MatAlloc(da.N_rows,1);
  Matrix da_i;
  da_i.N_rows = da.N_cols;
  da_i.N_cols = 1;
  int Size_pa = da.N_rows;
  double Z_a = 0;
  double Z_a_m1 = 0;

  /* Get Z and the numerator */
  for(int i = 0 ; i<Size_pa ; i++){
    da_i.nV = da.nM[i];
    pa.nV[i] = exp(LME_fa(da_i,lambda,Beta));
    Z_a += pa.nV[i];
  }

  /* Get the inverse of Z */
  Z_a_m1 = (double)1/Z_a;

  /* Divide by Z and get the final value */
  for(int i = 0 ; i<Size_pa ; i++){
    pa.nV[i] *= Z_a_m1;
  }
  
  /* Return the value of the shape function */  
  return pa;
}

Matrix LME_r(Matrix da, Matrix pa)
/*
  Gradient of the log(Z) function.
  Input parameters :
  - da : Matrix with the distances to the neighborhood nodes (neighborhood x dim).
  - pa : Shape function value in the neighborhood nodes (neighborhood x 1).
*/
{
  int N_dim = da.N_cols;
  int N_neibourg = da.N_rows;
  Matrix r = MatAllocZ(N_dim,1);

  /* Fill the ''r'' vector */
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
  Hessian of the log(Z) function.
  Input parameters :
  - da : Matrix with the distances to the neighborhood nodes (neighborhood x dim).
  - pa : Shape function value in the neighborhood nodes (neighborhood x 1).
  - r : Gradient of log(Z) (dim x 1).
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

double LME_dpa(){
  
}

ChainPtr LME_Tributary_Nodes(Matrix, int,
			     Matrix, Mesh){
}
