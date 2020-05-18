#include "nl-partsol.h"

/***********************************************/
/************* Local Maximum-Entropy **************/
/**************************************************/

/*!
  Shape functions based in :
  "" Local maximum-entropy approximation schemes : a seamless 
  bridge between finite elements and meshfree methods ""
  by M.Arroyo and M.Ortiz, 2006.

  Here we employ the same nomenclature as in the paper. With the single
  different of the "l" variable wich represents the distances between the
  evaluation point and the neighborhood nodes.

  List of functions :
  - LME_Init_lambda
  - LME_lambda_NR
  - LME_fa
  - LME_p
  - LME_r
  - LME_J
  - LME_dp
  - LME_Tributary_Nodes
*/


/****************************************************************************/

void LME_Initialize_Beta(Matrix Beta, double DeltaX, int Np)
/*!
  Function to get a initial value of Beta
*/
{
  int Ndim = NumberDimensions;
  
  for(int p = 0 ; p<Np ; p++){
    for(int i = 0 ; i<Ndim ; i++){
      Beta.nM[p][i] = gamma_LME/(DeltaX*DeltaX);
    }
  }
  
}

/****************************************************************************/

void LME_Initialize(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
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

  /* Initialize Beta */
  LME_Initialize_Beta(MPM_Mesh.Beta, FEM_Mesh.DeltaX, Np);

  /* Loop over the particle mesh */
  for(int p = 0 ; p<Np ; p++){

    /* Get some properties for each particle */ 
    X_p = get_RowFrom(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    Beta_p = get_RowFrom(Ndim,1,MPM_Mesh.Beta.nM[p]);
    lambda_p = get_RowFrom(Ndim,1,MPM_Mesh.lambda.nM[p]);

    /* Loop over the element mesh */
    for(int i = 0 ; i<Nelem ; i++){

      /* Get the element properties */
      Elem_p = FEM_Mesh.Connectivity[i];
      
      /* 6º Check out if the GP is in the Element */
      if(InOut_Element(X_p, Elem_p, FEM_Mesh.Coordinates)){

	/* With the element connectivity get the node close to the particle */
	I_p = get_closest_node_to(X_p,Elem_p,FEM_Mesh.Coordinates);

	/* Get the initial connectivity of the particle */
	MPM_Mesh.ListNodes[p] = LME_Tributary_Nodes(X_p,Beta_p,I_p,FEM_Mesh);

	/* Measure the size of the connectivity */
	MPM_Mesh.NumberNodes[p] = get_Lenght_Set(MPM_Mesh.ListNodes[p]);

	/* Asign to each particle the closest node in the mesh
	   and to this node asign the particle */
	MPM_Mesh.I0[p] = get_closest_node_to(X_p,MPM_Mesh.ListNodes[p],
					     FEM_Mesh.Coordinates);
	
	/* Active those nodes that interact with the particle */
	asign_particle_to_nodes(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
       	
	/* Calculate distance from particle to each node in the neibourhood */
	Delta_Xip = get_set_Coordinates(MPM_Mesh.ListNodes[p],X_p,
					FEM_Mesh.Coordinates);

	/* Update the value of beta */
	Beta_p = LME_Beta(Beta_p, Delta_Xip, gamma_LME);
			
	/* Update lagrange multipliers with Newton-Rapson */
	lambda_p = LME_lambda_NR(Delta_Xip, lambda_p, Beta_p);

	/* Free memory */
	FreeMat(Delta_Xip);

	break;
      }      
    }

  }

}

/****************************************************************************/

Matrix LME_Beta(Matrix Beta, Matrix l, double Gamma)
/*!
  Function to update the value of beta
*/
{

  int Ndim = NumberDimensions;
  int NumNodes_GP = l.N_rows;
  double h = 0;
  Matrix l_GP_I = get_RowFrom(Ndim,1,NULL);
  
  /* Get the mean distande */
  for(int i = 0 ; i<NumNodes_GP ; i++){
    l_GP_I.nV = l.nM[i];
    h += Norm_Mat(l_GP_I,2);
  }
  h = h/NumNodes_GP;

  /* Fill Beta */
  for(int j = 0 ; j<Ndim ; j++){
    Beta.nV[j] = Gamma/(h*h);
  }

  return Beta;
}

/****************************************************************************/

Matrix LME_lambda_NR(Matrix l, Matrix lambda, Matrix Beta)
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
    p = LME_p(l,lambda,Beta);

    /* Get the gradient of log(Z) and its norm */
    r = LME_r(l,p);
    norm_r = Norm_Mat(r,2);

    /* Get the Hessian of log(Z) and update it with +||r||*I 
       according with Dennis M.Kochmann et al. 2019 (CMAME) */    
    J = LME_J(l,p,r);
    for(int i = 0 ; i<Ndim ; i++){
      J.nM[i][i] += norm_r;
    }

    /* /\* Check the conditioning number of the Hessian *\/ */
    /* if (fabs(Cond_Mat(J,TOL_NR)) > 10){ */
    /*   printf(" %s (%s %i) : %s \n", */
    /* 	     "Error in LME_lambda_NR","Iter",NumIter, */
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
    FreeMat(p);
    FreeMat(r);
    FreeMat(J);
    FreeMat(D_lambda);
    
    if(norm_r > TOL_NR){
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

double LME_fa(Matrix la, Matrix lambda, Matrix Beta)
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

Matrix LME_p(Matrix l, Matrix lambda, Matrix Beta)
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
    MatAlloc(1,N_a); 
  Matrix la = /* Distance to the neighbour (x-x_a) */
    get_RowFrom(1,Ndim,NULL);

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

/****************************************************************************/

Matrix LME_r(Matrix l, Matrix p)
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
    = MatAllocZ(Ndim,1);

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

Matrix LME_J(Matrix l, Matrix p, Matrix r)
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
  J = MatAllocZ(Ndim,Ndim);

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

Matrix LME_dp(Matrix l, Matrix p)
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
  Matrix dp = MatAllocZ(N_a,Ndim);
  Matrix r; /* Gradient of log(Z) */
  Matrix J; /* Hessian of log(Z) */
  Matrix Jm1; /* Inverse of J */
  Matrix Jm1_la; /* Auxiliar vector */
  Matrix la = /* Distance to the neighbour (x-x_a) */
    get_RowFrom(Ndim,1,NULL); 
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = LME_r(l,p);
  J = LME_J(l,p,r);
  
  /* Inverse of the Hessian */
  Jm1 = Get_Inverse(J);
  
  /* Free memory */
  FreeMat(r);
  FreeMat(J);

  /* Fill the gradient for each node */
  for(int a = 0 ; a<N_a ; a++){
    la.nV = l.nM[a]; 
    Jm1_la = Scalar_prod(Jm1,la);    
    for(int i = 0 ; i<Ndim ; i++){
      dp.nM[a][i] = - p.nV[a]*Jm1_la.nV[i];
    }
    FreeMat(Jm1_la);
  }

  /* Free memory */
  FreeMat(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dp;  
}

/****************************************************************************/

ChainPtr LME_Tributary_Nodes(Matrix X_GP, Matrix Beta, int I0, Mesh FEM_Mesh){

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;
  /* Number of dimensionws of the problem */
  int Ndim = NumberDimensions;

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = get_RowFrom(Ndim,1,NULL);
  
  ChainPtr Set_Nodes0 = NULL;
  int * Array_Nodes0;
  int NumNodes0;
  int Node0;

  /* Get the search radius */
  double Ra = sqrt(-log(TOL_lambda)/Beta.nV[0]);

  /* Get nodes close to the particle */
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = Set_to_Pointer(Set_Nodes0,NumNodes0);
     
  /* Loop over the chain with the tributary nodes */
  for(int i = 0 ; i<NumNodes0 ; i++){

    Node0 = Array_Nodes0[i];

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

    /* Get a vector from the GP to the node */
    Distance = Sub_Mat(X_GP,X_I);

    /* If the node is near the GP push in the chain */
    if(Norm_Mat(Distance,2) <= Ra){
      push_to_Set(&Triburary_Nodes,Node0);
    }

    /* Free memory of the distrance vector */
    FreeMat(Distance);

  }
  
  /* Free memory */
  free(Array_Nodes0);

  return Triburary_Nodes;
}

