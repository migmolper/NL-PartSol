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
  -
  -
  -
*/
{
  /* Definition of some parameters */
  Matrix pa;
  Matrix r;
  Matrix J;
  Matrix Jm1;
  double norm_r = 10;
  int NumIter = 0;

  /* Start with the Newton-Rapson method */
  while(norm_r > TOL_lambda){
  
    /* Get vector with the shape functions evaluated in the nodes */
    pa = LME_pa(da,lambda,Beta);

    /* Get the gradient and the hessian */
    r = LME_r(da,pa);
    J = LME_J(da,pa,r);

    /* Check the conditioning number */
    if (fabs(Cond_Mat(J)) < 1e-8){
      printf(" %s : %s \n",
	     "Error in LME_lambda",
	     "The Hessian is near to singular matrix");      
      exit(0);
    }

    /* Free the shape function nodal values */
    FreeMat(pa);

    /* Get the norm of the gradient */
    norm_r = Norm_Mat(r,2);

    /* Inverse of the hessian */
    Jm1 = Get_Inverse(J);

    /* Update value of lambda */
    Increment_lambda = Scalar_prod(Jm1,r);
    lambda = Incr_Mat(lambda,Increment_lambda);

    /* Update the number of iterations */
    NumIter ++;
    if(NumIter > 100){
      printf(" %s : %s \n",
	     "Error in LME_lambda",
	     "No convergence in 100 iterations");
    }
  
  }

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
  Value of the shape function evaluated in the neibourhud;
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
  Gradient of the log(Z) function
*/
{
  int N_dim = da.N_cols;
  int N_neibourg = da.N_rows;
  Matrix r = MatAllocZ(N_dim,1);
  Matrix r_a;

  /* Fill the ''r'' vector */
  for(int i = 0 ; i<N_neibourg ; i++){
    for(int j = 0 ; j<N_dim ; j++){
      r.nV[j] += pa.nV[i]*da.nM[i][j];
    }
  }
  
  return r;
}

Matrix LME_J(Matrix da, Matrix pa, Matrix r)
/*
  Hessian of the log(Z) function
*/
{
  
  Matrix J;
  Matrix J_I;
  Matrix J_II;

  /* Get the first component of the hessian */
  for(){
  }

  /* Get the second component of the hessian */
  J_II = Tensorial_prod(r,r);

  /* Get the Hessian */
  J = Sub_Mat(J_I,J_II);
  
  return J;
}

double LME_dpa(){
  
}

ChainPtr LME_Tributary_Nodes(Matrix, int,
			     Matrix, Mesh){
}
