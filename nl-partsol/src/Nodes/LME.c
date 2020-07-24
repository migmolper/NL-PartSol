#include "nl-partsol.h"

/*
  Auxiliar functions 
 */
static double fa__LME__(Matrix, Matrix, Matrix);
static Matrix r__LME__(Matrix, Matrix);
static Matrix J__LME__(Matrix, Matrix, Matrix);

/****************************************************************************/

void initialize__LME__(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nelem = FEM_Mesh.NumElemMesh;

  /* Particle coordinates */  
  Matrix X_p;
  /* Lagrange multipliers */
  Matrix lambda_p;
  /* Tunning parameter */
  Matrix Beta_p; 
  
  ChainPtr Elem_p;
  
  ChainPtr Nodes_p;
  int I_p;

  /* Distance from GP to the nodes */
  Matrix Delta_Xip; 

  /* Loop over the particle mesh */
  for(int p = 0 ; p<Np ; p++){

    /* Get some properties for each particle */ 
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    Beta_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Beta.nM[p]);
    lambda_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lambda.nM[p]);

    /* Loop over the element mesh */
    for(int i = 0 ; i<Nelem ; i++){

      /* Get the element properties */
      Elem_p = FEM_Mesh.Connectivity[i];
      
      /* 6º Check out if the GP is in the Element */
      if(inout_convex_set__MeshTools__(X_p, Elem_p, FEM_Mesh.Coordinates)){

	/* With the element connectivity get the node close to the particle */
	I_p = get_closest_node__MeshTools__(X_p,Elem_p,FEM_Mesh.Coordinates);

	/* Calculate distance from particle to each node in the neibourhood */
	MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p);
	Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,
						  FEM_Mesh.Coordinates);

	/* Initialize Beta */
	Beta_p = beta_isotropic__LME__(Beta_p, Delta_Xip, gamma_LME);

	/* Get the initial connectivity of the particle */
	MPM_Mesh.ListNodes[p] = isotropic_tributary__LME__(X_p,Beta_p,I_p,FEM_Mesh);

	/* Measure the size of the connectivity */
	MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

	/* Asign to each particle the closest node in the mesh
	   and to this node asign the particle */
	MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,MPM_Mesh.ListNodes[p],
						       FEM_Mesh.Coordinates);
	
	/* Active those nodes that interact with the particle */
	asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
       	
	/* Calculate distance from particle to each node in the neibourhood */
	Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,
						  FEM_Mesh.Coordinates);

	/* Update the value of beta */
	Beta_p = beta_isotropic__LME__(Beta_p, Delta_Xip, gamma_LME);
			
	/* Update lagrange multipliers with Newton-Rapson */
	lambda_p = lambda__LME__(Delta_Xip, lambda_p, Beta_p);

	/* Free memory */
	free__MatrixLib__(Delta_Xip);

	break;
      }      
    }

  }

}

/****************************************************************************/

Matrix beta_isotropic__LME__(Matrix Beta, Matrix l, double Gamma)
{
  int Ndim = NumberDimensions;
  int NumNodes_GP = l.N_rows;
  double h = 0;
  Matrix l_GP_I = memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  
  /* Get the mean distande */
  for(int i = 0 ; i<NumNodes_GP ; i++){
    l_GP_I.nV = l.nM[i];
    h += norm__MatrixLib__(l_GP_I,2);
  }
  h = h/NumNodes_GP;

  /* Fill Beta */
  for(int j = 0 ; j<Ndim ; j++){
    Beta.nV[j] = Gamma/(h*h);
  }

  return Beta;
}

/****************************************************************************/

Matrix beta_anisotropic__LME__(Matrix Beta, Matrix f)
{

  int Ndim = NumberDimensions;

  Matrix f_m1 = inverse__MatrixLib__(f);

  Matrix Beta_x_f_m1 = matrix_product__MatrixLib__(Beta, f_m1);
  
  double f_m1T_x_Beta_x_f_m1;

  /*
    Compute and update :
    Beta = f_m1T_x_Beta_x_f_m1
   */
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{

	  f_m1T_x_Beta_x_f_m1 = 0;
	  
	  for(int k = 0 ; k<Ndim ; k++)
	    {
	      f_m1T_x_Beta_x_f_m1 += f_m1.nM[k][i]*Beta_x_f_m1.nM[k][j];
	    }
	  
	  /*
	    Update value of beta
	  */
	  Beta.nM[i][j] = f_m1T_x_Beta_x_f_m1;
	  
	}
    }

  /*
    Free memory
   */
  free__MatrixLib__(Beta_x_f_m1);
  free__MatrixLib__(f_m1);
    
  
  return Beta;
}

/****************************************************************************/

Matrix lambda__LME__(Matrix l, Matrix lambda, Matrix Beta)
/*!
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
  int MaxIter = 100;
  int Ndim = NumberDimensions;
  int NumIter = 0; /* Iterator counter */
  double norm_r = 10; /* Value of the norm */
  Matrix p; /* Shape function vector */
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix D_lambda; /* Increment of lambda */
 
  /* Start with the Newton-Rapson method */
  while(NumIter <= MaxIter){
	
    /* Get vector with the shape functions evaluated in the nodes */
    p = p__LME__(l,lambda,Beta);

    /* Get the gradient of log(Z) and its norm */
    r = r__LME__(l,p);
    norm_r = norm__MatrixLib__(r,2);

    /* Get the Hessian of log(Z) */    
    J = J__LME__(l,p,r);

    /* /\* Check the conditioning number of the Hessian *\/ */
    /* if (fabs(conditioning__MatrixLib__(J,TOL_lambda)) > 10){ */
    /*   printf(" %s (%s %i) : %s \n", */
    /* 	     "Error in lambda__LME__","Iter",NumIter, */
    /* 	     "The Hessian is near to singular matrix"); */
    /*   exit(0); */
    /* } */
    
    /* Get the increment of lambda */
    D_lambda = Solve_Linear_Sistem(J,r);

    /* Update the value of lambda */
    for(int i = 0 ; i<Ndim ; i++){
      lambda.nV[i] -= D_lambda.nV[i];
    }

    /* Free memory */
    free__MatrixLib__(p);
    free__MatrixLib__(r);
    free__MatrixLib__(J);
    free__MatrixLib__(D_lambda);
    
    if(norm_r > TOL_lambda){
      /* Update the iteration */
      NumIter ++;
    }
    else{
      /* Break the loop */
      break;
    }
    
  }

  if(NumIter == MaxIter){
    printf("%s : %s %i/%i %s \n",
	   "Warning in LME_lambda",
	   "No convergence in",NumIter,MaxIter,"iterations");
    printf("%s : %f \n",
	   "Error",norm_r);
  }
  
  /* Once the stopping criteria is reached, 
     return the lagrange multipliers value */
  return lambda;
}

/****************************************************************************/

static double fa__LME__(Matrix la, Matrix lambda, Matrix Beta)
/*!
  Output :
  -> fa : the function fa that appear in [1] (scalar).
  Input parameters :
  -> la : Matrix with the distance to the neighborhood node ''a'' (1 x dim).
  -> lambda : Initial value of the lagrange multipliers (dim x 1).
  -> Gamma : Tunning parameter (scalar).
*/
{  
  int Ndim = NumberDimensions;
  double fa = 0;

  for(int i = 0 ; i<Ndim ; i++){
    fa += - Beta.nV[i]*la.nV[i]*la.nV[i] + la.nV[i]*lambda.nV[i];
  }
    
  /* Return the value of fa */
  return fa;
}

/****************************************************************************/

Matrix p__LME__(Matrix l, Matrix lambda, Matrix Beta)
/*!
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
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0, Z_m1 = 0;
  Matrix p = /* Vector with the values of the shape-function in the nodes */
    alloc__MatrixLib__(1,N_a); 
  Matrix la = /* Distance to the neighbour (x-x_a) */
    memory_to_matrix__MatrixLib__(1,Ndim,NULL);

  /* Get Z and the numerator */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a];
    p.nV[a] = exp(fa__LME__(la,lambda,Beta));
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

/****************************************************************************/

static Matrix r__LME__(Matrix l, Matrix p)
/*!
  Get the gradient "r" (dim x 1) of the function log(Z) = 0.
  Input parameters :
  -> l : Matrix with the distances to the 
  neighborhood nodes (neighborhood x dim).
  -> p : Shape function value in the
  neighborhood nodes (1 x neighborhood).
*/
{  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix r /* Value of the gradient */
    = allocZ__MatrixLib__(Ndim,1);

  /* Fill ''r'' */
  for(int i = 0 ; i<Ndim ; i++){
    for(int a = 0 ; a<N_a ; a++){
      r.nV[i] += p.nV[a]*l.nM[a][i];
    }
  }

  /* Return the value of the gradient */
  return r;
}

/****************************************************************************/

static Matrix J__LME__(Matrix l, Matrix p, Matrix r)
/*!
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
  int Ndim = NumberDimensions;
  Matrix J; /* Hessian definition */
  
  /* Allocate Hessian */
  J = allocZ__MatrixLib__(Ndim,Ndim);

  /* Fill the Hessian */
  for(int i = 0 ; i<Ndim ; i++){
    for(int j = 0 ; j<Ndim ; j++){
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

/****************************************************************************/

Matrix dp__LME__(Matrix l, Matrix p)
/*!
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
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix dp = allocZ__MatrixLib__(N_a,Ndim);
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Jm1_la; /* Auxiliar vector */
  Matrix la = /* Distance to the neighbour (x-x_a) */
    memory_to_matrix__MatrixLib__(Ndim,1,NULL); 
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = r__LME__(l,p);
  J = J__LME__(l,p,r);
  
  /* Inverse of the Hessian */
  Jm1 = inverse__MatrixLib__(J);
  
  /* Free memory */
  free__MatrixLib__(r);
  free__MatrixLib__(J);

  /* Fill the gradient for each node */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a]; 
    Jm1_la = matrix_product__MatrixLib__(Jm1,la);    
    for(int i = 0 ; i<Ndim ; i++){
      dp.nM[a][i] = - p.nV[a]*Jm1_la.nV[i];
    }
    free__MatrixLib__(Jm1_la);
  }

  /* Free memory */
  free__MatrixLib__(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dp;  
}

/****************************************************************************/

ChainPtr isotropic_tributary__LME__(Matrix X_GP, Matrix Beta, int I0, Mesh FEM_Mesh)
{

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;
  /* Number of dimensionws of the problem */
  int Ndim = NumberDimensions;

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  
  ChainPtr Set_Nodes0 = NULL;
  int * Array_Nodes0;
  int NumNodes0;
  int Node0;

  /* Counter */
  int NumTributaryNodes = 0;

  /* Get the search radius */
  double Ra = sqrt(-log(TOL_lambda)/Beta.nV[0]);

  /* Get nodes close to the particle */
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0,NumNodes0);
     
  /* Loop over the chain with the tributary nodes */
  for(int i = 0 ; i<NumNodes0 ; i++){

    Node0 = Array_Nodes0[i];

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

    /* Get a vector from the GP to the node */
    Distance = substraction__MatrixLib__(X_GP,X_I);

    /* If the node is near the GP push in the chain */
    if(norm__MatrixLib__(Distance,2) <= Ra){
      push__SetLib__(&Triburary_Nodes,Node0);
      NumTributaryNodes++;
    }

    /* Free memory of the distrance vector */
    free__MatrixLib__(Distance);

  }

  /* If the Triburary_Nodes chain lenght is less than 3 assign al the node */
  if(NumTributaryNodes < Ndim + 1)
    {
      for(int i = 0 ; i<NumNodes0 ; i++)
	{
	  Node0 = Array_Nodes0[i];
	  if(!inout__SetLib__(Triburary_Nodes,Node0))
	    {
	      push__SetLib__(&Triburary_Nodes,Node0);
	    }
	}
    }
  
  /* Free memory */
  free(Array_Nodes0);

  return Triburary_Nodes;
}

/****************************************************************************/

ChainPtr anisotropic_tributary__LME__(Matrix X_p, Matrix M_p, Matrix f_p,
				      int I0, Mesh FEM_Mesh)
{

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;
  /* Number of dimensionws of the problem */
  int Ndim = NumberDimensions;

  Matrix Distance; /* Distance between node and GP */
  Matrix M_Distance;
  Matrix X_I = memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  
  ChainPtr Set_Nodes0 = NULL;
  int * Array_Nodes0;
  int NumNodes0;
  int Node0;

  /* Counter */
  int NumTributaryNodes = 0;

  /*
    Update the cut-off ellipsoid  
  */
  M_p = anisotropic_cut_off(M_p, f_p);
  double R_p;

  /* Get nodes close to the particle */
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0,NumNodes0);
     
  /* Loop over the chain with the tributary nodes */
  for(int i = 0 ; i<NumNodes0 ; i++){

    Node0 = Array_Nodes0[i];

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

    /*
      Check if the node is in/out of the ellipsoid 
    */
    Distance = substraction__MatrixLib__(X_p,X_I);
    M_Distance = matrix_product__MatrixLib__(M_p, Distance);
    R_p = scalar_product__MatrixLib__(Distance,M_Distance);

    /* If the node is near the GP push in the chain */
    if(R_p <= 1)
      {
	push__SetLib__(&Triburary_Nodes,Node0);
	NumTributaryNodes++;
      }

    /* Free memory of the distrance vector */
    free__MatrixLib__(Distance);
    free__MatrixLib__(M_Distance);

  }

  /* If the Triburary_Nodes chain lenght is less than 3 assign al the node */
  if(NumTributaryNodes < Ndim + 1)
    {
      for(int i = 0 ; i<NumNodes0 ; i++)
	{
	  Node0 = Array_Nodes0[i];
	  if(!inout__SetLib__(Triburary_Nodes,Node0))
	    {
	      push__SetLib__(&Triburary_Nodes,Node0);
	    }
	}
    }
  
  /* Free memory */
  free(Array_Nodes0);

  return Triburary_Nodes;
}

/****************************************************************************/

Matrix   anisotropic_cut_off(Matrix M, Matrix f)
{
  int Ndim = NumberDimensions;

  Matrix f_m1 = inverse__MatrixLib__(f);

  Matrix M_x_f_m1 = matrix_product__MatrixLib__(M, f_m1);
  
  double f_m1T_x_M_x_f_m1;

  /*
    Compute and update :
    M = f_m1T_x_M_x_f_m1
  */
  for(int i = 0 ; i<Ndim ; i++)
    {
      for(int j = 0 ; j<Ndim ; j++)
	{

	  f_m1T_x_M_x_f_m1 = 0;
	  
	  for(int k = 0 ; k<Ndim ; k++)
	    {
	      f_m1T_x_M_x_f_m1 += f_m1.nM[k][i]*M_x_f_m1.nM[k][j];
	    }
	  
	  /*
	    Update value of beta
	  */
	  M.nM[i][j] = f_m1T_x_M_x_f_m1;
	  
	}
    }

  /*
    Free memory
  */
  free__MatrixLib__(M_x_f_m1);
  free__MatrixLib__(f_m1);
    
  
  return M;
}
