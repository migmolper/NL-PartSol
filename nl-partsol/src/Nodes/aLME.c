#include "nl-partsol.h"

/*
  Call global variables
*/
double TOL_lambda;
double gamma_LME;

/*
  Define local global variable for the relative error
*/
double Error0;

/*
  Auxiliar functions 
 */
static Tensor initilise_beta__aLME__(Tensor, Matrix, double);
static Tensor initialise_cutoff__aLME__(Tensor, Matrix, double);
static double fa__aLME__(Tensor, Tensor, Tensor);
static Tensor r__aLME__(Matrix, Matrix);
static Tensor J__aLME__(Matrix, Matrix, Tensor);
static bool   check_convergence(Tensor,double,int,int);
static void   standard_error(char *);

/****************************************************************************/

void initialize__aLME__(GaussPoint MPM_Mesh, Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP;
  int Nelem = FEM_Mesh.NumElemMesh;

  /* Particle coordinates */  
  Matrix X_p;
  /* Lagrange multipliers */
  Tensor lambda_p;
  /* Tunning parameter */
  Tensor Beta_p; 
  /* Cutoff ellipsoid */
  Tensor M_p;
  /* Increment of deformation gradient */
  Tensor DF_p;
  
  ChainPtr Elem_p;
  
  ChainPtr Nodes_p;
  int I_p;

  /* Distance from GP to the nodes */
  Matrix Delta_Xip; 

  /* Loop over the particle mesh */
  for(int p = 0 ; p<Np ; p++)
  {

    /* Get some properties for each particle */ 
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);

    lambda_p = memory_to_tensor__TensorLib__(MPM_Mesh.lambda.nM[p],1);
    Beta_p   = memory_to_tensor__TensorLib__(MPM_Mesh.Beta.nM[p],2); 
    M_p      = memory_to_tensor__TensorLib__(MPM_Mesh.Cut_Off_Ellipsoid.nM[p],2); 
    DF_p     = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.DF.nM[p],2); 

    /* Loop over the element mesh */
    for(int i = 0 ; i<Nelem ; i++)
    {

      /* Get the element properties */
      Elem_p = FEM_Mesh.Connectivity[i];
      
      /* 6º Check out if the GP is in the Element */
      if(inout_convex_set__MeshTools__(X_p, Elem_p, FEM_Mesh.Coordinates))
      {

    	   /* With the element connectivity get the node close to the particle */
        I_p = get_closest_node__MeshTools__(X_p,Elem_p,FEM_Mesh.Coordinates);

    	   /* Calculate distance from particle to each node in the neibourhood */
        MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p);
        Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

        /* Initialize Beta */
        Beta_p = initilise_beta__aLME__(Beta_p, Delta_Xip, gamma_LME);

        /* Initilise cutoff ellipsoid */
        M_p = initialise_cutoff__aLME__(M_p, Delta_Xip, gamma_LME);

    	  /* Get the initial connectivity of the particle */
        MPM_Mesh.ListNodes[p] = tributary__aLME__(X_p,M_p,I_p,FEM_Mesh);

    	  /* Measure the size of the connectivity */
        MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

        /* Asign to each particle the closest node in the mesh
	       and to this node asign the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);

        /* Active those nodes that interact with the particle */
        asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);

        /* Calculate distance from particle to each node in the neibourhood */
        Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

	       /* Update the value of beta */
        Beta_p = beta__aLME__(Beta_p, DF_p);

      	/* Update lagrange multipliers with Newton-Rapson */
        lambda_p = lambda__aLME__(Delta_Xip, lambda_p, Beta_p);

      	/* Free memory */
        free__MatrixLib__(Delta_Xip);

       break;
     }      
   }

 }

}

/****************************************************************************/

Tensor beta__aLME__(Tensor Beta_n, Tensor f)
/*! 
  Using a covariant push forward, the Beta matrix tensor is updated
*/
{
  Tensor Beta_n1 = Beta_n;

  covariant_push_forward_tensor__TensorLib__(Beta_n1, Beta_n, f);

  return Beta_n1;
}

/****************************************************************************/

Tensor cut_off__aLME__(Tensor M_n, Tensor f)
{
  Tensor M_n1 = M_n;

  covariant_push_forward_tensor__TensorLib__(M_n1, M_n, f);

  return M_n1;
}

/****************************************************************************/

Tensor lambda__aLME__(Matrix l, Tensor lambda, Tensor Beta)
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
  
  int Ndim = NumberDimensions;  

  /* Definition of some parameters */
  int MaxIter = 100;
  int Iter = 0;
  bool Convergence = false;
  double TOL = TOL_lambda;


  double norm_r; /* Value of the norm */
  Matrix p; /* Shape function vector */
  Tensor r; /* Gradient of log(Z) */
  Tensor J; /* Hessian of log(Z) */
  Tensor D_lambda; /* Increment of lambda */
 
  /* Start with the Newton-Rapson method */
  while(Convergence == false)
  {
	
    /* Get vector with the shape functions evaluated in the nodes */
    p = p__aLME__(l,lambda,Beta);

    /* Get the gradient of log(Z) and check the convergence */
    r = r__aLME__(l,p);
    Convergence = check_convergence(r, TOL, Iter, MaxIter);

    if(Convergence == false)
    {
      /* Get the Hessian of log(Z) */    
      J = J__aLME__(l,p,r);
    
      /* Get the increment of lambda */
      D_lambda = Solve_system__TensorLib__(J,r);

      /* Update the value of lambda */
      for(int i = 0 ; i<Ndim ; i++)
      {
        lambda.n[i] -= D_lambda.n[i];
      }

      /* Free memory */
      free__MatrixLib__(p);
      free__TensorLib__(J);
      free__TensorLib__(D_lambda);
      free__TensorLib__(r);

      /* Update iterations */
      Iter++;
    }
    
  }

  /* Free memory */
  free__MatrixLib__(p);
  free__TensorLib__(r);
  
  /* Return the lagrange multipliers value */
  return lambda;
}

/****************************************************************************/

Matrix p__aLME__(Matrix l, Tensor lambda, Tensor Beta)
/*!
  Get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.

  Input parameters :
  -> l : Matrix with the distances to the neighborhood nodes (neiborghood x dim).
  -> lambda : Initial value of the lagrange multipliers (1 x dim).
  -> Beta : Tunning parameter (scalar).
*/
{
  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double Z_m1 = 0;
  Tensor la; /* Distance to the neighbour (x-x_a) */

  /* Vector with the values of the shape-function in the nodes */
  Matrix p = alloc__MatrixLib__(1,N_a); 

  /* Get Z and the numerator */
  for(int a = 0 ; a<N_a ; a++)
  {
    la = memory_to_tensor__TensorLib__(l.nM[a], 1);
    p.nV[a] = exp(fa__aLME__(la,lambda,Beta));
    Z += p.nV[a];
  }

  /* Get the inverse of Z */
  Z_m1 = (double)1/Z;

  /* Divide by Z and get the final value */
  for(int a = 0 ; a<N_a ; a++)
  {
    p.nV[a] *= Z_m1;
  }
  
  /* Return the value of the shape function */  
  return p;
}

/****************************************************************************/

Matrix dp__aLME__(Matrix l, Matrix p)
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

  /* Auxiliar variables */
  Tensor r;   /* Gradient of log(Z) */
  Tensor J;   /* Hessian of log(Z) */
  Tensor Jm1; /* Inverse of J */
  Tensor la; /* Vector with the nodal distance to the node a */ 
  Tensor Jm1_la; /* Auxiliar vector */
  
  /* Get the Gradient and the Hessian of log(Z) */
  r = r__aLME__(l,p);
  J = J__aLME__(l,p,r);
  
  /* Inverse of the Hessian */
  Jm1 = Inverse__TensorLib__(J);

  /* Fill the gradient for each node */
  for(int a = 0 ; a<N_a ; a++)
  {
    la = memory_to_tensor__TensorLib__(l.nM[a], 1);
    Jm1_la = vector_linear_mapping__TensorLib__(Jm1, la);
    
    for(int i = 0 ; i<Ndim ; i++)
    {
      dp.nM[a][i] = - p.nV[a]*Jm1_la.n[i];
    }

    free__TensorLib__(Jm1_la);
  }

  /* Free memory */
  free__TensorLib__(r);
  free__TensorLib__(J);
  free__TensorLib__(Jm1);
  
  /* Return the value of the shape function gradient */  
  return dp;  
}

/****************************************************************************/

ChainPtr tributary__aLME__(Matrix Coordinate_p, Tensor M_p, int I0, Mesh FEM_Mesh)
{

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;
  /* Number of dimensionws of the problem */
  int Ndim = NumberDimensions;

  /* Distance between node and GP */
  Tensor X_p;
  Tensor X_I;
  Tensor x_Ip;
  
  ChainPtr Set_Nodes0 = NULL;
  int * Array_Nodes0;
  int NumNodes0;
  int Node0;

  /* Counter */
  int NumTributaryNodes = 0;

  double R_p;

  X_p = memory_to_tensor__TensorLib__(Coordinate_p.nV, 1);

  /* Get nodes close to the particle */
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0,NumNodes0);
     
  /* Loop over the chain with the tributary nodes */
  for(int i = 0 ; i<NumNodes0 ; i++)
  {

    Node0 = Array_Nodes0[i];

    /* Assign to a pointer the coordinates of the nodes */
    X_I = memory_to_tensor__TensorLib__(FEM_Mesh.Coordinates.nM[Node0], 1);

    /* Check if the node is in/out of the ellipsoid */
    x_Ip = subtraction__TensorLib__(X_p,X_I);

    R_p  = Generalised_norm__TensorLib__(x_Ip, M_p);
    
    /* If the node is near the GP push in the chain */
    if(R_p <= 1.0)
      {
	       push__SetLib__(&Triburary_Nodes,Node0);
	       NumTributaryNodes++;
      }

    /* Free memory */
    free__TensorLib__(x_Ip);

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

static Tensor initilise_beta__aLME__(Tensor Beta_p, Matrix l, double Gamma)
{
  int Ndim = NumberDimensions;
  int NumNodes_p = l.N_rows;
  double h = 0;
  Tensor l_pI;
  
  /* Get the mean distande */
  for(int i = 0 ; i<NumNodes_p ; i++)
  {
    l_pI = memory_to_tensor__TensorLib__(l.nM[i],1);
    h += EuclideanNorm__TensorLib__(l_pI);
  }

  h = h/NumNodes_p;

  /* Fill Beta */
  for(int i = 0 ; i<Ndim ; i++)
  {
    Beta_p.N[i][i] = Gamma/(h*h);
  }

  return Beta_p;
}

/****************************************************************************/

static Tensor initialise_cutoff__aLME__(Tensor M_p, Matrix l, double Gamma)
{
  int Ndim = NumberDimensions;
  int NumNodes_p = l.N_rows;
  double h = 0;
  Matrix l_p_I = memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  double TOL = TOL_lambda;
  
  /* Get the mean distande */
  for(int i = 0 ; i<NumNodes_p ; i++)
  {
    l_p_I.nV = l.nM[i];
    h += norm__MatrixLib__(l_p_I,2);
  }

  h = h/NumNodes_p;

  double Ra2 = -log(TOL)*DSQR(h)/Gamma;

  /* Fill the cutoff ellipsoid */
  for(int i = 0 ; i<Ndim ; i++)
  {
    M_p.N[i][i] = 1/Ra2;
  }

  return M_p;
}

/**************************************************************/

static double fa__aLME__(Tensor la, Tensor lambda, Tensor Beta)
/*!
  Output :
  -> fa : the function fa that appear in [1] (scalar).
  Input parameters :
  -> la : Vector with the distance to the neighborhood node ''a'' (1 x dim).
  -> lambda : VEctor with the initial value of the lagrange multipliers (dim x 1).
  -> Beta : 
*/
{  
  int Ndim = NumberDimensions;
  double fa = 0;

  Tensor Beta_x_la = vector_linear_mapping__TensorLib__(Beta,la);

  for(int i = 0 ; i<Ndim ; i++)
  {
    fa += - Beta_x_la.n[i]*la.n[i] + la.n[i]*lambda.n[i];
  }

  free__TensorLib__(Beta_x_la);
    
  /* Return the value of fa */
  return fa;
}

/****************************************************************************/

static Tensor r__aLME__(Matrix l, Matrix p)
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
  Tensor r = alloc__TensorLib__(1);

  /* Fill ''r'' */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int a = 0 ; a<N_a ; a++)
    {
      r.n[i] += p.nV[a]*l.nM[a][i];
    }
  }

  /* Return the value of the gradient */
  return r;
}

/****************************************************************************/

static Tensor J__aLME__(Matrix l, Matrix p, Tensor r)
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
  
  /* Hessian */
  Tensor J = alloc__TensorLib__(2);

  /* Fill the Hessian */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      for(int a = 0 ; a<N_a ; a++)
      {
      /* Get the first component of the Hessian looping
        over the neighborhood nodes. */
         J.N[i][j] += p.nV[a]*l.nM[a][i]*l.nM[a][j];
      }
      /* Get the second value of the Hessian */
      J.N[i][j] -= r.n[i]*r.n[j];
    }
  }

  /* Return the value of the Hessian */
  return J;
}

/****************************************************************************/

static bool check_convergence(Tensor r, double TOL, int Iter, int MaxIter)
{
  bool convergence = false;
  double Error_relative;
  char Error_message[MAXW];

  double Error = EuclideanNorm__TensorLib__(r);

  if(Iter > MaxIter)
    {
      sprintf(Error_message,"%s","Convergence not reached in the maximum number of iterations");
      standard_error(Error_message); 
    }

  /*
    Compute relative error
  */
  if(Iter == 0)
    {
      Error0 = Error;
      Error_relative = Error/Error0;
//        printf("Error iter %i : %1.4e ; %1.4e \n",Iter,Error,Error_relative);
    }
    else
    {
      Error_relative = Error/Error0;
//        printf("Error iter %i : %1.4e ; %1.4e \n",Iter,Error,Error_relative);
    }
      
    /*
      Check convergence using the relative error
    */
  if(Error_relative < TOL)
    {
      convergence = true;
    }

  return convergence;
}


/****************************************************************************/

static void standard_error(char * Error_message)
{
  fprintf(stderr,"%s : %s !!! \n", "Error in lambda__aLME__",Error_message);
  exit(EXIT_FAILURE);
}

/****************************************************************************/

