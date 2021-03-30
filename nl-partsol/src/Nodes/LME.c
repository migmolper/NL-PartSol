#include "nl-partsol.h"

/*
  Auxiliar functions 
 */
static double fa__LME__(Matrix, Matrix, Matrix, double);
static Matrix r__LME__(Matrix, Matrix);
static Matrix J__LME__(Matrix, Matrix, Matrix);

/*
  Call global var¡ables
*/
double gamma_LME;
double TOL_LME;
char * Metric_LME;

/****************************************************************************/

void initialize__LME__(
  GaussPoint MPM_Mesh,
  Mesh FEM_Mesh)
{

  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP; // Number of gauss-points in the simulation
  int Nelem = FEM_Mesh.NumElemMesh; // Number of elements
  int I0; // Closest node to the particle
  ChainPtr Elem_p; // Surrounding elements
  ChainPtr Nodes_p; // Surrounding particles

  Matrix Metric_p; // Define a metric tensor
  Matrix X_p; // Particle coordinates  
  Matrix Delta_Xip; // Distance from GP to the nodes
  Matrix lambda_p; // Lagrange multiplier
  Matrix F_p; // Particle deformation gradient, only for anysotropic
  double Beta_p; // Thermalization or regularization parameter


  for(int p = 0 ; p<Np ; p++)
  {
    /* 
      Get some properties for each particle
    */ 
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    lambda_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lambda.nM[p]);
    Beta_p = MPM_Mesh.Beta.nV[p];

    /*
      Get the metric tensor
    */
    if(strcmp(Metric_LME,"Identity") == 0)
    {
      Metric_p = metric_I__LME__();
    }
    else if(strcmp(Metric_LME,"bm1") == 0)
    {
      F_p = memory_to_matrix__MatrixLib__(Ndim,Ndim,MPM_Mesh.Phi.F_n.nM[p]);
      Metric_p = metric_bm1__LME__(F_p);
    }

    /*
      Loop over the element mesh
    */
    for(int i = 0 ; i<Nelem ; i++)
    {

      /* Get the element properties */
      Elem_p = FEM_Mesh.Connectivity[i];
      
      /* Check out if the GP is in the Element */
      if(inout_convex_set__MeshTools__(X_p, Elem_p, FEM_Mesh.Coordinates))
      {

        /* With the element connectivity get the node close to the particle */
        I0 = get_closest_node__MeshTools__(X_p,Elem_p,FEM_Mesh.Coordinates);

        /* Calculate distance from particle to each node in the neibourhood */
        MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p);
        Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

        /* Initialize Beta */
        Beta_p = beta__LME__(Delta_Xip, gamma_LME, FEM_Mesh.DeltaX);

        /* Get the initial connectivity of the particle */
        MPM_Mesh.ListNodes[p] = tributary__LME__(X_p,Metric_p,Beta_p,I0,FEM_Mesh);

        /* Measure the size of the connectivity */
        MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

        /* Asign to each particle the closest node in the mesh
          and to this node asign the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,MPM_Mesh.ListNodes[p],FEM_Mesh.Coordinates);

        /* 
          Active those nodes that interact with the particle
        */
        asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
       	
        /*
          Calculate distance from particle to each node in the neibourhood
        */
        Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

        /*
          Update the value of the thermalization parameter
        */
        Beta_p = beta__LME__(Delta_Xip, gamma_LME, FEM_Mesh.DeltaX);
        MPM_Mesh.Beta.nV[p] = Beta_p;

        /* 
          Update lagrange multipliers with Newton-Rapson
        */
        lambda_p = lambda__LME__(Delta_Xip, lambda_p, Metric_p, Beta_p);

        /* 
          Free memory
        */
        free__MatrixLib__(Delta_Xip);

        break;
      }      
    }

    /*
      Free memory
    */
    free__MatrixLib__(Metric_p);
  }

}

/****************************************************************************/

double beta__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  double Gamma, // User define parameter to control the value of the thermalization parameter.
  double DeltaX) // Average mesh size
/*!
  Get the thermalization parameter beta using the global variable gamma_LME.
*/
{
  int Ndim = NumberDimensions;
  int NumNodes_GP = l.N_rows;
  double Beta = 0; // Intialise the thermalization parameter
  double avg_l = 0; // Initalise the average nodal distance
  double h = 0; // Distance parameter
  Matrix l_pI = memory_to_matrix__MatrixLib__(Ndim,1,NULL);
  
  /* 
    Get the mean distande
  */
  for(int i = 0 ; i<NumNodes_GP ; i++){
    l_pI.nV = l.nM[i];
    avg_l += norm__MatrixLib__(l_pI,2);
  }
  avg_l = avg_l/NumNodes_GP;


  h = DMIN(avg_l,DeltaX);

  /*
    Compute beta
  */
  Beta = Gamma/(h*h);

  return Beta;
}

/****************************************************************************/

Matrix metric_I__LME__()
/*!
  Return a metric tensor to compute the locality parameter
  in the LME shape functions
*/
{
  int Ndim = NumberDimensions;
  Matrix Metric = allocZ__MatrixLib__(Ndim,Ndim);

  for(int i = 0 ; i<Ndim ; i++)
  {
    Metric.nM[i][i] = 1.0;
  }

  return Metric;
}

/****************************************************************************/

 Matrix metric_bm1__LME__(Matrix F)
 /*!
   Return the metric tensor proposed by Molinos (b^{-1} = F^{-T}F^{-1})
 */
 {
    int Ndim = NumberDimensions;
    Matrix Metric = allocZ__MatrixLib__(Ndim,Ndim);

    Matrix Fm1 = inverse__MatrixLib__(F);

    for(int i = 0 ; i < Ndim  ; i++)
    {
      for(int j = 0 ; j < Ndim  ; j++)
      {
        for(int k = 0 ; k < Ndim  ; k++)
        {
          Metric.nM[i][j] += Fm1.nM[k][i]*Fm1.nM[k][j];
        }
      }
    } 

    free__MatrixLib__(Fm1);

   return Metric;
 }

/****************************************************************************/

Matrix lambda__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix lambda, // Lagrange multiplier.
  Matrix Metric, // Measure for the norm definition.
  double Beta) // Thermalization parameter.
/*!
  Get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method for that is the Newton-Rapson.
*/
{
  /*
    Definition of some parameters
  */
  int MaxIter = 100;
  int Ndim = NumberDimensions;
  int NumIter = 0; // Iterator counter
  double norm_r = 10; // Value of the norm
  Matrix p; // Shape function vector
  Matrix r; // Gradient of log(Z)
  Matrix J; // Hessian of log(Z)
  Matrix D_lambda; // Increment of lambda
 

  while(NumIter <= MaxIter)
  {
	
    /* 
      Get vector with the shape functions evaluated in the nodes
    */
    p = p__LME__(l,lambda,Metric,Beta);

    /*
      Get the gradient of log(Z) and its norm
    */
    r = r__LME__(l,p);
    norm_r = norm__MatrixLib__(r,2);

    /* 
      Check convergence
    */
    if(norm_r > TOL_LME)
    {
      /* 
        Get the Hessian of log(Z)
      */    
      J = J__LME__(l,p,r);
    
      /*
        Get the increment of lambda
      */
      D_lambda = Solve_Linear_Sistem(J,r);

      /*
        Update the value of lambda
      */
      for(int i = 0 ; i<Ndim ; i++)
      {
        lambda.nV[i] -= D_lambda.nV[i];
      }

      /*
        Free memory
      */
      free__MatrixLib__(p);
      free__MatrixLib__(r);
      free__MatrixLib__(J);
      free__MatrixLib__(D_lambda);

      NumIter ++;
    }
    else
    {
      free__MatrixLib__(r);
      break;
    }
    
  }

  if(NumIter >= MaxIter)
  {
    printf("%s : %s \n",
      "Warning in LME_lambda",
      "No convergence reached in the maximum number of interations");
    printf("%s : %f\n", "Total Error",norm_r);
  }
  

  return lambda;
}

/****************************************************************************/

static double fa__LME__(
  Matrix la, // Vector form node "a" to particle.
  Matrix lambda, // Lagrange multiplier.
  Matrix Metric, // Measure for the norm definition.
  double Beta) // Thermalization parameter.
/*!
  fa (scalar): the function fa that appear in [1].
*/
{  
  int Ndim = NumberDimensions;
  double la_x_lambda = 0;
  double Metric_x_la = 0;
  double norm_la = 0;
  double fa = 0;

  for(int i = 0 ; i<Ndim ; i++)
  {

    Metric_x_la = 0;

    for(int j = 0 ; j<Ndim ; j++)
    {
      Metric_x_la += Metric.nM[i][j]*la.nV[j];
    }

    norm_la += la.nV[i]*Metric_x_la;

    la_x_lambda += la.nV[i]*lambda.nV[i];
    
  }
  
  fa = - Beta*norm_la + la_x_lambda;

  return fa;
}

/****************************************************************************/

Matrix p__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix lambda, // Lagrange multiplier.
  Matrix Metric, // Measure for the norm definition.
  double Beta) // Thermalization parameter.
/*!
  Get the value of the shape function "pa" (1 x neighborhood) in the
  neighborhood nodes.
*/
{
  
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double Z_m1 = 0;
  Matrix p = alloc__MatrixLib__(1,N_a); // Shape function in the nodes
  Matrix la = memory_to_matrix__MatrixLib__(1,Ndim,NULL); // Vector form node "a" to particle.

  /*
    Get Z and the numerator
  */
  for(int a = 0 ; a<N_a ; a++)
  {
    la.nV = l.nM[a];
    p.nV[a] = exp(fa__LME__(la,lambda,Metric,Beta));
    Z += p.nV[a];
  }

  /*
    Divide by Z and get the final value of the shape function
  */
  Z_m1 = (double)1/Z;
  for(int a = 0 ; a<N_a ; a++)
  {
    p.nV[a] *= Z_m1;
  }
  

  return p;
}

/****************************************************************************/

static Matrix r__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix p) // Set with the evaluation of the shape function in the neighborhood nodes.
/*!
  Gradient dlogZ_dLambda "r"
*/
{  
  /* 
    Definition of some parameters
  */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix r = allocZ__MatrixLib__(Ndim,1); // Gradient definition

  /* 
    Fill the gradient
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int a = 0 ; a<N_a ; a++)
    {
      r.nV[i] += p.nV[a]*l.nM[a][i];
    }
  }

  return r;
}

/****************************************************************************/

static Matrix J__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix p, // Set with the evaluation of the shape function in the neighborhood nodes.
  Matrix r) // Gradient dlogZ_dLambda "r"
/*!
  Hessian d2logZ_dLambdadLambda "J"
*/
{  
  /* 
    Definition of some parameters
  */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix J = allocZ__MatrixLib__(Ndim,Ndim); // Hessian definition
  
  /*
    Fill the Hessian
  */
  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      /* 
        Get the first component of the Hessian looping
        over the neighborhood nodes.
      */
      for(int a = 0 ; a<N_a ; a++)
      {
        J.nM[i][j] += p.nV[a]*l.nM[a][i]*l.nM[a][j];
      }

      /*
        Get the second value of the Hessian
      */
      J.nM[i][j] -= r.nV[i]*r.nV[j];
    }
  }

  return J;
}

/****************************************************************************/

Matrix dp__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix p) // Set with the evaluation of the shape function in the neighborhood nodes.
/*!
  Value of the shape function gradient "dp" (dim x neighborhood) in the neighborhood nodes
*/
{  
  /* 
    Definition of some parameters
  */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix dp = allocZ__MatrixLib__(N_a,Ndim);
  Matrix r; // Gradient of log(Z)
  Matrix J; // Hessian of log(Z)
  Matrix Jm1; // Inverse of J
  Matrix Jm1_la; // Auxiliar vector
  Matrix la = memory_to_matrix__MatrixLib__(Ndim,1,NULL); // Distance to the neighbour (x-x_a)
  
  /*
    Get the Gradient and the Hessian of log(Z)
  */
  r = r__LME__(l,p);
  J = J__LME__(l,p,r);
  
  /*
    Inverse of the Hessian
  */
  Jm1 = inverse__MatrixLib__(J);
  
  /* 
    Fill the gradient for each node
  */
  for(int a = 0 ; a<N_a ; a++)
  {
    la.nV = l.nM[a]; 
    Jm1_la = matrix_product__MatrixLib__(Jm1,la);    
  
    for(int i = 0 ; i<Ndim ; i++)
    {
      dp.nM[a][i] = - p.nV[a]*Jm1_la.nV[i];
    }

    free__MatrixLib__(Jm1_la);
  }

  /* 
    Free memory
  */
  free__MatrixLib__(r);
  free__MatrixLib__(J);
  free__MatrixLib__(Jm1);

  return dp;  
}

/****************************************************************************/

ChainPtr tributary__LME__(
  Matrix X_p,
  Matrix Metric,
  double Beta_p,
  int I0,
  Mesh FEM_Mesh)
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
  double Ra = sqrt(-log(TOL_LME)/Beta_p);

  /* Get nodes close to the particle */
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0,NumNodes0);
     
  /* Loop over the chain with the tributary nodes */
  for(int i = 0 ; i<NumNodes0 ; i++)
  {

    Node0 = Array_Nodes0[i];

    /* Assign to a pointer the coordinates of the nodes */
    X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

    /* Get a vector from the GP to the node */
    Distance = substraction__MatrixLib__(X_p,X_I);

    /* If the node is near the GP push in the chain */
    if(generalised_Euclidean_distance__MatrixLib__(Distance, Metric) <= Ra)
    {
      push__SetLib__(&Triburary_Nodes,Node0);
      NumTributaryNodes++;
    }

    /* Free memory of the distrance vector */
    free__MatrixLib__(Distance);

  }


  /* 
    If the Triburary_Nodes chain lenght is less than 3 assign al the node
  */
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
