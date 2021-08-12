#include "nl-partsol.h"

/*!
  \fn double fa__aLME__(Matrix la,Matrix lambda,Matrix Beta)
  \brief Function with computes the scalar values which contains the restrictions related 
  with the widh and information entropy of the shape function.
  \param la : Vector form node "a" to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.
*/
static double fa__aLME__(Matrix, Matrix, Matrix);

/****************************************************************************/

/*!
  \fn double logZ__aLME__(Matrix l,Matrix lambda,Matrix Beta)
  \brief Compute the value of the objetive function to compute the lagrange multiplier.
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.
*/
static double logZ__aLME__(Matrix, Matrix, Matrix);

/****************************************************************************/

/*!
  \fn Matrix r__aLME__(Matrix l,Matrix p);
  \brief Compute the value of the gradient of the objetive function to compute the lagrange multiplier.
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param p : Value of the shape function in the nodes.
*/
static Matrix r__aLME__(Matrix, Matrix);

/****************************************************************************/

/*!
  \fn Matrix J__aLME__(Matrix l,Matrix p,Matrix r)
  \brief Compute the hessian matrix for the function logZ
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param p : Set with the evaluation of the shape function in the neighborhood nodes.
  \param r : Gradient of the function logZ_dLambda "r".
*/
static Matrix J__aLME__(Matrix, Matrix, Matrix);

/****************************************************************************/

/*!  
 \fn void initialize_beta__aLME__(Matrix Beta, double Gamma, double DeltaX); 
 \brief Compute the value of the thermalization parameter using a circular support. 
 \param Beta : Termalization tensor 
 \param Gamma : Adimensional paramter to control the regularization parameter. 
 \param DeltaX : Minimum size in the all nodal set.
*/
static void initialize_beta__aLME__(Matrix,double,double);

/****************************************************************************/

/*!  
 \fn void update_beta__aLME__(Matrix Beta, Matrix Delta_F); 
 \brief Update the termalization matrix employing the increment of the deformation gradient. 
 \param Beta : Termalization matrix.
 \param Delta_F : Increment of the deformation gradient.
*/
static void update_beta__aLME__(Matrix,Matrix);

/****************************************************************************/

/*!
  \fn void initialize_Cut_off_Ellipsoid__aLME__(Matrix Cut_off_Ellipsoid,double Gamma,double h_avg);
  \brief Initializae the definition of the cutt-off matrix
  \param Gamma : User define parameter to control the value of the thermalization parameter.
  \param h_avg : Average mesh size
*/
static void initialize_Cut_off_Ellipsoid__aLME__(Matrix,double,double);

/****************************************************************************/

/*!  
 \fn void update_cut_off_ellipsoid__aLME__(Matrix Cut_off_Ellipsoid, Matrix Delta_F); 
 \brief  Update the shape of the cut-off region employing the increment of the deformation gradient. 
 \param Cut_off_Ellipsoid : Cut-off matrix.
 \param Delta_F : Increment of the deformation gradient.
*/
static void update_cut_off_ellipsoid__aLME__(Matrix,Matrix);

/****************************************************************************/

/*!
\fn void initialise_lambda__aLME__(int Idx_particle,Matrix X_p,Matrix Elem_p_Coordinates,Matrix lambda,Matrix Beta)
  \param Idx_particle : Index of the particle.
  \param X_p : Particle position.
  \param Elem_p_Coordinates : Coordinates of the nodes close to the particle.
  \param lambda: Lagrange multiplier of the particle.
  \param Beta: Thermalization matrix.
*/
static void initialise_lambda__aLME__(int,Matrix,Matrix,Matrix,Matrix);

/****************************************************************************/

/*!
  \fn void update_lambda_Newton_Rapson__aLME__(int p,Matrix l,Matrix lambda,Matrix Beta)
  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical methodis the Newton-Rapson.
  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.
*/
static void update_lambda_Newton_Rapson__aLME__(int,Matrix,Matrix,Matrix);

/****************************************************************************/

/*!
  \fn void lambda_Nelder_Mead__aLME__(int p,Matrix l,Matrix lambda,Matrix Beta)
  \brief Function ot get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method is the Nelder-Mead.
  \param p : Current particle
  \param l : Set than contanins vector form neighborhood nodes to particle.
  \param lambda : Lagrange multiplier.
  \param Beta : Thermalization tensor.
*/
void update_lambda_Nelder_Mead__aLME__(int, Matrix, Matrix, Matrix);

// Auxiliar functions for the Neldel Mead in the LME
static Matrix gravity_center_Nelder_Mead__aLME__(Matrix);
static void order_logZ_simplex_Nelder_Mead__aLME__(Matrix, Matrix);
static void expansion_Nelder_Mead__aLME__(Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,double);
static void contraction_Nelder_Mead__aLME__(Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,double);
static void shrinkage_Nelder_Mead__aLME__(Matrix,Matrix,Matrix,Matrix);

// Nelder-Mead parameters
double NM_rho_aLME = 1.0;
double NM_chi_aLME = 2.0;
double NM_gamma_aLME = 0.5;
double NM_sigma_aLME = 0.5;
double NM_tau_aLME = 1E-3;

/****************************************************************************/

/*!
  \fn ChainPtr tributary__aLME__(int Indx_p,Matrix X_p,Matrix Cut_off_Ellipsoid,int I0,Mesh FEM_Mesh);
  \brief Update the list of tributary nodes for each particle
  \param Indx_p : Index of the particle.
  \param X_p : Position of the particle.
  \param Cut_off_Ellipsoid : Cut-off matrix.
  \param I0 : Index of the closest node to the particle.
  \param FEM_Mesh : Information of the background set of nodes.
*/
static ChainPtr tributary__aLME__(int,Matrix, Matrix, int, Mesh);

/****************************************************************************/

// Call global var¡ables
char wrapper_LME[MAXC];
double gamma_LME;
double TOL_zero_LME;
double TOL_wrapper_LME;
int max_iter_LME;

/****************************************************************************/

void initialize__aLME__(
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{
  int Ndim = NumberDimensions;
  int Np = MPM_Mesh.NumGP; // Number of gauss-points in the simulation
  int Nelem = FEM_Mesh.NumElemMesh; // Number of elements
  int I0; // Closest node to the particle
  ChainPtr Elem_p_Connectivity; // Surrounding elements
  Matrix Elem_p_Coordinates;
  ChainPtr Nodes_p; // Surrounding particles
  bool Init_p;
  Matrix X_p; // Particle coordinates  
  Matrix Delta_Xip; // Distance from GP to the nodes
  Matrix lambda_p; // Lagrange multiplier
  Matrix Cut_off_Ellipsoid; // Ellipsoid
  Matrix Beta_p; // Thermalization or regularization tensor
  ChainPtr Locality_I0; // List of nodes close to the node I0_p

  for(int p = 0 ; p<Np ; p++)
  {

    // Supose that the particle was not initilise
    Init_p = false;

    // Get some properties for each particle
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    Beta_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Beta.nM[p]);

    // Loop over the element mesh
    for(int i = 0 ; i<Nelem ; i++)
    {

      // Get the element properties
      Elem_p_Connectivity = FEM_Mesh.Connectivity[i];
      Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(Elem_p_Connectivity, FEM_Mesh.Coordinates);

      // Check out if the GP is in the Element
      if(FEM_Mesh.In_Out_Element(X_p,Elem_p_Coordinates) == true)
      {

        // Assign the index of the element
        MPM_Mesh.Element_p[p] = i;

        // Particle will be initilise
        Init_p = true;

        // Asign to each particle the closest node in the mesh and to this node asign the particle
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p_Connectivity,FEM_Mesh.Coordinates);

        // Initialize Beta
        initialize_beta__aLME__(Beta_p, gamma_LME, FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);

        // Initialise lambda for the Nelder-Mead using Bo-Li approach
        if(strcmp(wrapper_LME,"Nelder-Mead") == 0)
        {
          initialise_lambda__aLME__(p, X_p, Elem_p_Coordinates, lambda_p, Beta_p);
        }

        // Select the closest nodes to the particle and activate them
        Locality_I0 = FEM_Mesh.NodalLocality_0[MPM_Mesh.I0[p]];

        while(Locality_I0 != NULL)
        {
          if(FEM_Mesh.ActiveNode[Locality_I0->I] == false)
          {
            FEM_Mesh.ActiveNode[Locality_I0->I] = true;
          }

          Locality_I0 = Locality_I0->next; 

        }

        free__MatrixLib__(Elem_p_Coordinates);
       	
        break;
      }      

      free__MatrixLib__(Elem_p_Coordinates);

    }

    if(!Init_p)
    {
      fprintf(stderr,"%s : %s %i\n",
        "Error in initialize__aLME__()",
        "The search algorithm was unable to find particle",p);
      exit(EXIT_FAILURE);
    }

  } 

  for(int p = 0 ; p<Np ; p++)
  {

    // Get some properties for each particle
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    lambda_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lambda.nM[p]);
    
    // Get the metric tensor and initialize it
    Cut_off_Ellipsoid = memory_to_matrix__MatrixLib__(Ndim,Ndim,MPM_Mesh.Cut_off_Ellipsoid.nM[p]);
    initialize_Cut_off_Ellipsoid__aLME__(Cut_off_Ellipsoid,gamma_LME,FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);

    // Get the initial connectivity of the particle
    MPM_Mesh.ListNodes[p] = tributary__aLME__(p,X_p,Cut_off_Ellipsoid,MPM_Mesh.I0[p],FEM_Mesh);

    // Calculate number of nodes
    MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

    // Generate nodal distance list
    Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

    // Update lagrange multiplier with Newton-Rapson or with Nelder-Mead
    if(strcmp(wrapper_LME,"Newton-Raphson") == 0)
    {
      update_lambda_Newton_Rapson__aLME__(p, Delta_Xip, lambda_p, Beta_p);
    }
    else if(strcmp(wrapper_LME,"Nelder-Mead") == 0)
    {
      update_lambda_Nelder_Mead__aLME__(p, Delta_Xip, lambda_p, Beta_p);
    }
    else
    {
      fprintf(stderr,"%s : %s \n","Error in local_search__aLME__","Unrecognaised wrapper");
      exit(EXIT_FAILURE);      
    }

    // Active those nodes that interact with the particle
    asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);

    free__MatrixLib__(Delta_Xip);
  }


}

/****************************************************************************/

static void initialize_beta__aLME__(
  Matrix Beta,
  double Gamma,
  double h_avg)
{
  int Ndim = NumberDimensions;
  double aux = Gamma/(h_avg*h_avg);

  for(int i = 0 ; i<Ndim ; i++)
  {
    Beta.nM[i][i] = aux;
  }
}

/****************************************************************************/

static void initialize_Cut_off_Ellipsoid__aLME__(
  Matrix Cut_off_Ellipsoid,
  double Gamma,
  double h_avg)
{
  int Ndim = NumberDimensions;
  double aux = Gamma/(-log(TOL_zero_LME)*h_avg*h_avg);

  for(int i = 0 ; i<Ndim ; i++)
  {
    Cut_off_Ellipsoid.nM[i][i] = aux;
  }
}

/****************************************************************************/

static void initialise_lambda__aLME__(
  int Idx_particle,
  Matrix X_p,
  Matrix Elem_p_Coordinates,
  Matrix lambda,
  Matrix Beta)
{

  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;
  int Size_element = Elem_p_Coordinates.N_rows;
  double sqr_dist_i;

  int * simplex;

  Matrix Norm_l = allocZ__MatrixLib__(Size_element,1);
  Matrix l = allocZ__MatrixLib__(Size_element,Ndim);

  Matrix A = allocZ__MatrixLib__(Ndim,Ndim);
  Matrix b = allocZ__MatrixLib__(Ndim,1);
  Matrix x;

  // Initialise a list with distances and order
  for(int i = 0 ; i<Size_element ; i++)
  {

    sqr_dist_i = 0.0;

    for(int j = 0 ; j<Ndim ; j++)
    {
      l.nM[i][j] = X_p.nV[i] - Elem_p_Coordinates.nM[i][j];
      sqr_dist_i += DSQR(l.nM[i][j]);
    }

    Norm_l.nV[i] = sqr_dist_i;
  
  }

  if(Size_element == 3)
  {
    simplex = (int *)Allocate_ArrayZ(Nnodes_simplex,sizeof(int));
    simplex[0] = 0;
    simplex[1] = 1;
    simplex[2] = 2;
  }
  else if(Size_element == 4)
  {
    simplex = (int *)Allocate_ArrayZ(Nnodes_simplex,sizeof(int));
    simplex[0] = 0;
    simplex[1] = 1;
    simplex[2] = 2;
  }
  else
  {
    exit(0);
  }

  // Assemble matrix to solve the system Ax = b
  for(int i = 1 ; i<Nnodes_simplex ; i++)
  {

    b.nV[i-1] = - Beta.nM[i][i]*(Norm_l.nV[simplex[0]] - Norm_l.nV[simplex[i]]);

    for(int j = 0 ; j<Ndim ; j++)
    {
      A.nM[i-1][j] = l.nM[simplex[i]][j] - l.nM[simplex[0]][j];
    }
  }

  // Solve the system
  if(rcond__MatrixLib__(A) < 1E-8)
  {
    fprintf(stderr,"%s %i : %s \n",
      "Error in initialise_lambda__aLME__ for particle",
      Idx_particle,"The Hessian near to singular matrix!");
    exit(EXIT_FAILURE);
  }

  x = solve__MatrixLib__(A,b);

  // Update the value of lambda
  for(int i = 0 ; i<Ndim ; i++)
  {
    lambda.nV[i] = x.nV[i];
  }

  free(simplex);
  free__MatrixLib__(Norm_l);
  free__MatrixLib__(l);
  free__MatrixLib__(A);
  free__MatrixLib__(b);
  free__MatrixLib__(x);
}

/****************************************************************************/

static void update_lambda_Newton_Rapson__aLME__(
  int Idx_particle,
  Matrix l,
  Matrix lambda,
  Matrix Beta)
{
  int MaxIter = max_iter_LME;
  int Ndim = NumberDimensions;
  int NumIter = 0; // Iterator counter.
  double norm_r = 10; // Value of the norm.
  Matrix p; // Shape function vector.
  Matrix r; // Gradient of log(Z).
  Matrix J; // Hessian of log(Z).
  Matrix D_lambda; // Increment of lambda.
 
  while(NumIter <= MaxIter)
  {
	
    // Get vector with the shape functions evaluated in the nodes.
    p = p__aLME__(l,lambda,Beta);

    // Get the gradient of log(Z) and its norm.
    r = r__aLME__(l,p);
    norm_r = norm__MatrixLib__(r,2);

    // Check convergence.
    if(norm_r > TOL_wrapper_LME)
    {
      // Get the Hessian of log(Z).
      J = J__aLME__(l,p,r);

      if(rcond__MatrixLib__(J) < 1E-8)
      {
        fprintf(stderr,"%s %i : %s \n",
          "Error in lambda_Newton_Rapson__aLME__ for particle",
          Idx_particle,"The Hessian near to singular matrix!");
        exit(EXIT_FAILURE);
      }
    
      // Get the increment of lambda.
      D_lambda = solve__MatrixLib__(J,r);

      // Update the value of lambda.
      for(int i = 0 ; i<Ndim ; i++)
      {
        lambda.nV[i] -= D_lambda.nV[i];
      }

      free__MatrixLib__(p);
      free__MatrixLib__(r);
      free__MatrixLib__(J);
      free__MatrixLib__(D_lambda);

      NumIter ++;
    }
    else
    {
      free__MatrixLib__(r);
      free__MatrixLib__(p);
      break;
    }
    
  }

  if(NumIter >= MaxIter)
  {
    fprintf(stderr,"%s %i : %s (%i)\n",
      "Warning in lambda_Newton_Rapson__aLME__ for particle",Idx_particle,
      "No convergence reached in the maximum number of interations",MaxIter);
    fprintf(stderr,"%s : %e\n", "Total Error",norm_r);
    exit(EXIT_FAILURE);
  }
  
}

/****************************************************************************/

void update_lambda_Nelder_Mead__aLME__(
  int Idx_particle,
  Matrix l,
  Matrix lambda,
  Matrix Beta)
{
  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;
  int MaxIter = max_iter_LME;
  int NumIter = 0;

  // Simplex generated with lagrange multipliers
  Matrix simplex = allocZ__MatrixLib__(Nnodes_simplex, Ndim); 
  // Vector with the evaluation of the objective function in each vertex of the symplex
  Matrix logZ = allocZ__MatrixLib__(Nnodes_simplex, 1); 
  // Auxiliar variables
  Matrix simplex_a = memory_to_matrix__MatrixLib__(1,Ndim,NULL);
  Matrix gravity_center;
  Matrix reflected_point;
  double logZ_reflected_point;
  double logZ_0;
  double logZ_n;
  double logZ_n1;

  // Compute the initial positions of the nodes in the simplex (P.Navas)
  for(int a = 0 ; a<Nnodes_simplex ; a++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      if(i == a)
      {
        simplex.nM[a][i] = lambda.nV[i]/10;
      }
      else
      {
        simplex.nM[a][i] = lambda.nV[i];
      }
    }
  }

  // Compute the initial values of logZ in each vertex of the simplex.
  for(int a = 0 ; a<Nnodes_simplex ; a++)
  {    
    simplex_a.nV = simplex.nM[a];
    logZ.nV[a] = logZ__aLME__(l, simplex_a, Beta);
  }

  // Nelder-Mead main loop.
  while(NumIter <= MaxIter)
  {

    order_logZ_simplex_Nelder_Mead__aLME__(logZ, simplex);

    logZ_0 = logZ.nV[0];
    logZ_n = logZ.nV[Nnodes_simplex-2];
    logZ_n1 = logZ.nV[Nnodes_simplex-1];

    // Check convergence.
    if(fabs(logZ_0 - logZ_n1) > TOL_wrapper_LME)
    {

      // Spin the simplex to get the simplex with the smallest normalized volume.
      // spin_Nelder_Mead__aLME__(simplex);

      // Compute the gravity center of the simplex.
      gravity_center = gravity_center_Nelder_Mead__aLME__(simplex);

      // Compute the reflected point and evaluate the objetive function in this point.
      reflected_point = allocZ__MatrixLib__(1, Ndim);

      for(int i = 0 ; i<Ndim ; i++)
      {
        reflected_point.nV[i] = gravity_center.nV[i] + NM_rho_aLME*(gravity_center.nV[i] - simplex.nM[Nnodes_simplex-1][i]);
      }

      logZ_reflected_point = logZ__aLME__(l,reflected_point,Beta);

      // Do an expansion using the reflected point.
      if(logZ_reflected_point < logZ_0)
      {
        expansion_Nelder_Mead__aLME__(simplex,logZ,reflected_point,gravity_center,l,Beta,logZ_reflected_point);
      }
      // Take the reflected point.
      else if((logZ_reflected_point > logZ_0) && (logZ_reflected_point < logZ_n))
      {
        for(int i = 0 ; i<Ndim ; i++)
        {
          simplex.nM[Nnodes_simplex-1][i] = reflected_point.nV[i];
        }

        logZ.nV[Nnodes_simplex-1] = logZ_reflected_point; 

      }
      // Do a contraction using the reflected point (or a shrinkage).
      else if(logZ_reflected_point >= logZ_n)
      {
        contraction_Nelder_Mead__aLME__(simplex,logZ,reflected_point,gravity_center,l,Beta,logZ_reflected_point);
      }

      free__MatrixLib__(reflected_point);  
      NumIter ++;

    }
    else
    {
      break;
    }

  }

  if(NumIter >= MaxIter)
  {
    fprintf(stderr,"%s %i : %s (%i) \n",
      "Warning in lambda_Nelder_Mead__aLME__ for particle",Idx_particle,
      "No convergence reached in the maximum number of interations",MaxIter);
    fprintf(stderr,"%s : %e\n", "Total Error",fabs(logZ_0 - logZ_n1));
    exit(EXIT_FAILURE);
  }

  // Update the value of lambda.
  for(int i = 0 ; i<Ndim ; i++)
  {
    lambda.nV[i] = simplex.nM[0][i];
  }

  free__MatrixLib__(simplex);
  free__MatrixLib__(logZ);

}


/****************************************************************************/

static void order_logZ_simplex_Nelder_Mead__aLME__(
  Matrix logZ, 
  Matrix simplex)
{

  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;
  bool swapped = false;
  double aux;

  // Ordenate the list from lowest to higher (bubble sort)
  for(int i = 1 ; i<Nnodes_simplex ; i++)
  {
    swapped = false;

    for(int j = 0 ; j<(Nnodes_simplex - i) ; j++)
    {

      if(logZ.nV[j] > logZ.nV[j+1]) 
      {
      
        // swap values of logZ
        aux = logZ.nV[j];
        logZ.nV[j] = logZ.nV[j+1];
        logZ.nV[j+1] = aux;

        // swap values of logZ
        for(int k = 0 ; k<Ndim ; k++)
        {
          aux = simplex.nM[j][k];
          simplex.nM[j][k] = simplex.nM[j+1][k];
          simplex.nM[j+1][k] = aux;
        }

        swapped = true;
      }

    }

    if(!swapped) 
    {
      break;
    }

  }

}

/****************************************************************************/

static Matrix gravity_center_Nelder_Mead__aLME__(Matrix simplex)
{
  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;

  Matrix gravity_center = allocZ__MatrixLib__(1, Ndim);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int a = 0 ; a<Nnodes_simplex ; a++)
    {
      gravity_center.nV[i] += simplex.nM[a][i]/Nnodes_simplex;
    }
  }

  return gravity_center;
}

/****************************************************************************/

static void expansion_Nelder_Mead__aLME__(
  Matrix simplex,
  Matrix logZ,
  Matrix reflected_point,
  Matrix gravity_center,
  Matrix l,
  Matrix Beta, 
  double logZ_reflected_point)
{
  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;
  double logZ_expanded_point;
  Matrix expanded_point;

  // Compute the expanded point and evaluate the objetive function in this point
  expanded_point = allocZ__MatrixLib__(1, Ndim);
  
  for(int i = 0 ; i<Ndim ; i++)
  {
    expanded_point.nV[i] = gravity_center.nV[i] + NM_chi_aLME*(reflected_point.nV[i] - gravity_center.nV[i]);
  }

  logZ_expanded_point = logZ__aLME__(l, expanded_point, Beta); 

  // Take the expanded point
  if(logZ_expanded_point < logZ_reflected_point)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      simplex.nM[Nnodes_simplex-1][i] = expanded_point.nV[i];
    }

    logZ.nV[Nnodes_simplex-1] = logZ_expanded_point;
  }

  // Take the reflected point
  else
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      simplex.nM[Nnodes_simplex-1][i] = reflected_point.nV[i];
    }

    logZ.nV[Nnodes_simplex-1] = logZ_reflected_point;
  }

  free__MatrixLib__(expanded_point);
}

/****************************************************************************/

static void contraction_Nelder_Mead__aLME__(
  Matrix simplex,
  Matrix logZ,
  Matrix reflected_point,
  Matrix gravity_center,
  Matrix l,
  Matrix Beta, 
  double logZ_reflected_point)
{
  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;
  double logZ_n1 = logZ.nV[Nnodes_simplex-1];
  double logZ_contracted_point;
  Matrix contracted_point;

  contracted_point = allocZ__MatrixLib__(1, Ndim);

  // External contraction
  if(logZ_reflected_point < logZ_n1)
  {

    for(int i = 0 ; i<Ndim ; i++)
    {
      contracted_point.nV[i] = gravity_center.nV[i] + NM_gamma_aLME*(reflected_point.nV[i] - gravity_center.nV[i]);
    }

    logZ_contracted_point = logZ__aLME__(l, contracted_point, Beta); 

    // Take the contracted point
    if(logZ_contracted_point < logZ_reflected_point)
    {
      for(int i = 0 ; i<Ndim ; i++)
      {
        simplex.nM[Nnodes_simplex-1][i] = contracted_point.nV[i];
      }

      logZ.nV[Nnodes_simplex-1] = logZ_contracted_point;
    }
    // Do a shrinkage
    else
    {
      shrinkage_Nelder_Mead__aLME__(simplex,logZ,l,Beta);
    }
  }
  // Internal contraction
  else if(logZ_reflected_point > logZ_n1)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      contracted_point.nV[i] = gravity_center.nV[i] - NM_gamma_aLME*(gravity_center.nV[i] - simplex.nM[Nnodes_simplex-1][i]);
    }

    logZ_contracted_point = logZ__aLME__(l, contracted_point, Beta); 

    // Take the contracted point
    if(logZ_contracted_point < logZ_n1)
    {
      for(int i = 0 ; i<Ndim ; i++)
      {
        simplex.nM[Nnodes_simplex-1][i] = contracted_point.nV[i];
      }

      logZ.nV[Nnodes_simplex-1] = logZ_contracted_point;
    }
    // Do a shrinkage
    else
    {
      shrinkage_Nelder_Mead__aLME__(simplex,logZ,l,Beta);
    }
  }

  free__MatrixLib__(contracted_point);

}

/****************************************************************************/

static void shrinkage_Nelder_Mead__aLME__(
  Matrix simplex,
  Matrix logZ,
  Matrix l,
  Matrix Beta)
{
  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;

  // Axiliar function to get the coordinates of the simplex in the node a
  Matrix simplex_a = memory_to_matrix__MatrixLib__(1,Ndim,NULL);

  for (int a = 0 ; a<Nnodes_simplex ; a++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      simplex.nM[a][i] = simplex.nM[0][i] + NM_sigma_aLME*(simplex.nM[a][i] - simplex.nM[0][i]);
    }
    simplex_a.nV = simplex.nM[a];
    logZ.nV[a] = logZ__aLME__(l, simplex_a, Beta);
  }
}

/****************************************************************************/

static double fa__aLME__(
  Matrix la,
  Matrix lambda,
  Matrix Beta)
{  
  int Ndim = NumberDimensions;
  double la_x_lambda = 0;
  double Beta_x_la = 0;
  double Beta_norm_la = 0;
  double fa = 0;

  for(int i = 0 ; i<Ndim ; i++)
  {

    Beta_x_la = 0;

    for(int j = 0 ; j<Ndim ; j++)
    {
      Beta_x_la += Beta.nM[i][j]*la.nV[j];
    }

    Beta_norm_la += la.nV[i]*Beta_x_la;

    la_x_lambda += la.nV[i]*lambda.nV[i];
    
  }
  
  fa = - Beta_norm_la + la_x_lambda;

  return fa;
}

/****************************************************************************/

Matrix p__aLME__(
  Matrix l,
  Matrix lambda,
  Matrix Beta)
{
  
  // Definition of some parameters
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double Z_m1 = 0;
  Matrix p = allocZ__MatrixLib__(1,N_a); // Shape function in the nodes
  Matrix la = memory_to_matrix__MatrixLib__(1,Ndim,NULL); // Vector form node "a" to particle.

  // Get Z and the numerator.
  for(int a = 0 ; a<N_a ; a++)
  {
    la.nV = l.nM[a];
    p.nV[a] = exp(fa__aLME__(la,lambda,Beta));
    Z += p.nV[a];
  }

  // Divide by Z and get the final value of the shape function.
  Z_m1 = (double)1/Z;
  for(int a = 0 ; a<N_a ; a++)
  {
    p.nV[a] *= Z_m1;
  }
  

  return p;
}

/****************************************************************************/

static double logZ__aLME__(
  Matrix l,
  Matrix lambda,
  Matrix Beta)
{
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double logZ = 0;
  Matrix la = memory_to_matrix__MatrixLib__(1,Ndim,NULL); 

  for(int a = 0 ; a<N_a ; a++)
  {
    la.nV = l.nM[a];
    Z += exp(fa__aLME__(la,lambda,Beta));
  }

  logZ = log(Z);

  return logZ;
}

/****************************************************************************/

static Matrix r__aLME__(
  Matrix l,
  Matrix p)
{  
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix r = allocZ__MatrixLib__(Ndim,1); // Gradient definition

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

static Matrix J__aLME__(
  Matrix l,
  Matrix p,
  Matrix r)
{  
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix J = allocZ__MatrixLib__(Ndim,Ndim);

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {

      for(int a = 0 ; a<N_a ; a++)
      {
        J.nM[i][j] += p.nV[a]*l.nM[a][i]*l.nM[a][j];
      }

      J.nM[i][j] -= r.nV[i]*r.nV[j];
    }
  }

  return J;
}

/****************************************************************************/

Matrix dp__aLME__(
  Matrix l,
  Matrix p)
{  
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix dp = allocZ__MatrixLib__(N_a,Ndim);
  Matrix r; // Gradient of log(Z)
  Matrix J; // Hessian of log(Z)
  Matrix Jm1; // Inverse of J
  Matrix Jm1_la; // Auxiliar vector
  Matrix la = memory_to_matrix__MatrixLib__(Ndim,1,NULL); // Distance to the neighbour (x-x_a)
  
  // Get the Gradient and the Hessian of log(Z).
  r = r__aLME__(l,p);
  J = J__aLME__(l,p,r);
  
  // Inverse of the Hessian.
  Jm1 = inverse__MatrixLib__(J);
  
  // Fill the gradient for each node.
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

  free__MatrixLib__(r);
  free__MatrixLib__(J);
  free__MatrixLib__(Jm1);

  return dp;  
}

/****************************************************************************/

void local_search__aLME__(
  Particle MPM_Mesh,
  Mesh FEM_Mesh)
{
  int Ndim = NumberDimensions;
  int Num_Particles_Node_i;
  int Num_Particles_Element_i;
  Matrix X_p; // Velocity of the particle.
  Matrix V_p; // Velocity of the particle.
  ChainPtr Locality_I0; // List of nodes close to the node I0_p

  // Set to zero the active/non-active node, and the GPs in each element
  for(int i = 0 ; i<FEM_Mesh.NumNodesMesh ; i++)
  {
    Num_Particles_Node_i = FEM_Mesh.Num_Particles_Node[i];

    if(Num_Particles_Node_i != 0)
    {
      FEM_Mesh.Num_Particles_Node[i] = 0;
      free__SetLib__(&FEM_Mesh.List_Particles_Node[i]);
    }

    FEM_Mesh.ActiveNode[i] = false;

  }

  for(int i = 0 ; i<FEM_Mesh.Num_Patch_Mesh ; i++)
  {
    FEM_Mesh.Vol_Patch_n[i] = 0.0;
    FEM_Mesh.Vol_Patch_n1[i] = 0.0;
  }


  // Loop over the particles to create the list with active nodes.
  for(int p = 0 ; p<MPM_Mesh.NumGP ; p++)
  {

    // Get the global coordinates and velocity of the particle.
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);
    V_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.vel.nM[p]);


    // Update the index of the node close to the particle if there is advection.
    if(norm__MatrixLib__(V_p,2) > 0)
    {
      // Update the index of the closest node to the particle.
      Locality_I0 = FEM_Mesh.NodalLocality_0[MPM_Mesh.I0[p]];
      MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Locality_I0,FEM_Mesh.Coordinates);

      // Search particle in the sourrounding elements to this node.
      MPM_Mesh.Element_p[p] = search_particle_in_surrounding_elements__Particles__(p,X_p,FEM_Mesh.NodeNeighbour[MPM_Mesh.I0[p]],FEM_Mesh);
    }
  
    // Select the closest nodes to the particle.
     Locality_I0 = FEM_Mesh.NodalLocality_0[MPM_Mesh.I0[p]];

    // Activate the nodes near the particle.
    while(Locality_I0 != NULL)
    {
      if(FEM_Mesh.ActiveNode[Locality_I0->I] == false)
      {
        FEM_Mesh.ActiveNode[Locality_I0->I] = true;
      }
      Locality_I0 = Locality_I0->next; 
    }

  }

  // Loop over the particles to compute the tributary nodes.
  for(int p = 0 ; p<MPM_Mesh.NumGP ; p++)
  {
    Matrix Metric_p; // Define a metric tensor
    Matrix Delta_Xip; // Distance from particles to the nodes
    Matrix lambda_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.lambda.nM[p]);
    Matrix DF_p = memory_to_matrix__MatrixLib__(Ndim,Ndim,MPM_Mesh.Phi.DF.nM[p]); // Particle deformation gradient
    Matrix Beta_p = memory_to_matrix__MatrixLib__(Ndim,Ndim,MPM_Mesh.Beta.nM[p]);
    Matrix Cut_off_Ellipsoid_p = memory_to_matrix__MatrixLib__(Ndim,Ndim,MPM_Mesh.Cut_off_Ellipsoid.nM[p]);
    
    // Get the global coordinates of the particle.
    X_p = memory_to_matrix__MatrixLib__(Ndim,1,MPM_Mesh.Phi.x_GC.nM[p]);

    // Update the termalization matrix and cut-off ellipsoid.
    update_beta__aLME__(Beta_p, DF_p);
    update_cut_off_ellipsoid__aLME__(Cut_off_Ellipsoid_p,DF_p);

    // Free previous list of tributary nodes to the particle.
    free__SetLib__(&MPM_Mesh.ListNodes[p]);
    MPM_Mesh.ListNodes[p] = NULL;

    // Calculate the new connectivity with the previous value of beta.
    MPM_Mesh.ListNodes[p] = tributary__aLME__(p,X_p,Cut_off_Ellipsoid_p,MPM_Mesh.I0[p],FEM_Mesh);

    // Calculate number of nodes.
    MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

    // Generate nodal distance list.
    Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p], X_p, FEM_Mesh.Coordinates);            

    // Update lagrange multiplier with Newton-Rapson or with Nelder-Mead.
    if(strcmp(wrapper_LME,"Newton-Raphson") == 0)
    {
      update_lambda_Newton_Rapson__aLME__(p, Delta_Xip, lambda_p, Beta_p);
    }
    else if(strcmp(wrapper_LME,"Nelder-Mead") == 0)
    {
      update_lambda_Nelder_Mead__aLME__(p, Delta_Xip, lambda_p, Beta_p);
    }
    else
    {
      fprintf(stderr,"%s : %s \n","Error in local_search__aLME__","Unrecognaised wrapper");
      exit(EXIT_FAILURE);      
    }
    
    free__MatrixLib__(Metric_p);
    free__MatrixLib__(Delta_Xip);

    // Active those nodes that interact with the particle.
    asign_to_nodes__Particles__(p, MPM_Mesh.Element_p[p], MPM_Mesh.I0[p], MPM_Mesh.ListNodes[p], FEM_Mesh);
  }


}

/****************************************************************************/

static void update_beta__aLME__(
   Matrix Beta,
   Matrix Delta_F)
 {
    int Ndim = NumberDimensions;
    Matrix Delta_F_m1 = inverse__MatrixLib__(Delta_F);
    Matrix updated_Beta = allocZ__MatrixLib__(Ndim,Ndim);

#if NumberDimensions == 2
    double aux_00 = Beta.nM[0][0]*Delta_F_m1.nM[0][0] + Beta.nM[0][1]*Delta_F_m1.nM[1][0];
    double aux_01 = Beta.nM[0][0]*Delta_F_m1.nM[0][1] + Beta.nM[0][1]*Delta_F_m1.nM[1][1];
    double aux_10 = Beta.nM[1][0]*Delta_F_m1.nM[0][0] + Beta.nM[1][1]*Delta_F_m1.nM[1][0];
    double aux_11 = Beta.nM[1][0]*Delta_F_m1.nM[0][1] + Beta.nM[1][1]*Delta_F_m1.nM[1][1];

    updated_Beta.nM[0][0] = 
      Delta_F_m1.nM[0][0]*Beta.nM[0][0]*Delta_F_m1.nM[0][0] + 
      Delta_F_m1.nM[0][0]*Beta.nM[0][1]*Delta_F_m1.nM[1][0] +
      Delta_F_m1.nM[1][0]*Beta.nM[1][0]*Delta_F_m1.nM[0][0] + 
      Delta_F_m1.nM[1][0]*Beta.nM[1][1]*Delta_F_m1.nM[1][0];

    updated_Beta.nM[0][1] = Delta_F_m1.nM[0][0]*aux_01 + Delta_F_m1.nM[1][0]*aux_11;
    updated_Beta.nM[1][0] = Delta_F_m1.nM[0][1]*(Beta.nM[0][0]*Delta_F_m1.nM[0][0] + Beta.nM[0][1]*Delta_F_m1.nM[1][0]) + Delta_F_m1.nM[1][1]*aux_10;
    updated_Beta.nM[1][1] = Delta_F_m1.nM[0][1]*aux_01 + Delta_F_m1.nM[1][1]*aux_11;
#endif

#if NumberDimensions == 3 
    fprintf(stderr,"%s : %s !!! \n",
        "Error in update_beta__aLME__()",
        "This operation it is not implemented for 3D");
    exit(EXIT_FAILURE);  
#endif


  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Beta.nM[i][j] = updated_Beta.nM[i][j];
    }
  }

  free__MatrixLib__(Delta_F_m1);
  free__MatrixLib__(updated_Beta);

 }

/****************************************************************************/

static void update_cut_off_ellipsoid__aLME__(
   Matrix Cut_off_Ellipsoid,
   Matrix Delta_F)
 {
    int Ndim = NumberDimensions;
    Matrix Delta_F_m1 = inverse__MatrixLib__(Delta_F);
    Matrix updated_Cut_off_Ellipsoid = allocZ__MatrixLib__(Ndim,Ndim);

#if NumberDimensions == 2
    double aux_00 = Cut_off_Ellipsoid.nM[0][0]*Delta_F_m1.nM[0][0] + Cut_off_Ellipsoid.nM[0][1]*Delta_F_m1.nM[1][0];
    double aux_01 = Cut_off_Ellipsoid.nM[0][0]*Delta_F_m1.nM[0][1] + Cut_off_Ellipsoid.nM[0][1]*Delta_F_m1.nM[1][1];
    double aux_10 = Cut_off_Ellipsoid.nM[1][0]*Delta_F_m1.nM[0][0] + Cut_off_Ellipsoid.nM[1][1]*Delta_F_m1.nM[1][0];
    double aux_11 = Cut_off_Ellipsoid.nM[1][0]*Delta_F_m1.nM[0][1] + Cut_off_Ellipsoid.nM[1][1]*Delta_F_m1.nM[1][1];

    updated_Cut_off_Ellipsoid.nM[0][0] = Delta_F_m1.nM[0][0]*aux_00 + Delta_F_m1.nM[1][0]*aux_10;
    updated_Cut_off_Ellipsoid.nM[0][1] = Delta_F_m1.nM[0][0]*aux_01 + Delta_F_m1.nM[1][0]*aux_11;
    updated_Cut_off_Ellipsoid.nM[1][0] = Delta_F_m1.nM[0][1]*aux_00 + Delta_F_m1.nM[1][1]*aux_10;
    updated_Cut_off_Ellipsoid.nM[1][1] = Delta_F_m1.nM[0][1]*aux_01 + Delta_F_m1.nM[1][1]*aux_11;
#endif

#if NumberDimensions == 3 
    fprintf(stderr,"%s : %s !!! \n",
        "Error in update_cut_off_ellipsoid__aLME__()",
        "This operation it is not implemented for 3D");
    exit(EXIT_FAILURE);  
#endif

  for(int i = 0 ; i<Ndim ; i++)
  {
    for(int j = 0 ; j<Ndim ; j++)
    {
      Cut_off_Ellipsoid.nM[i][j] = updated_Cut_off_Ellipsoid.nM[i][j];
    }
  }

  free__MatrixLib__(Delta_F_m1);
  free__MatrixLib__(updated_Cut_off_Ellipsoid);

 }

/****************************************************************************/

static ChainPtr tributary__aLME__(
  int Indx_p,
  Matrix X_p,
  Matrix Cut_off_Ellipsoid,
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

  // Counter
  int NumTributaryNodes = 0;

  // Get nodes close to the particle
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0,NumNodes0);

  // Loop over the chain with the tributary nodes
  for(int i = 0 ; i<NumNodes0 ; i++)
  {

    Node0 = Array_Nodes0[i];

    if(FEM_Mesh.ActiveNode[Node0] == true)
    {
      // Assign to a pointer the coordinates of the nodes
      X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

      // Get a vector from the GP to the node
      Distance = substraction__MatrixLib__(X_p,X_I);

      // If the node is near the GP push in the chain 
      if(generalised_Euclidean_distance__MatrixLib__(Distance, Cut_off_Ellipsoid) <= 1.0)
      {
        push__SetLib__(&Triburary_Nodes,Node0);
        NumTributaryNodes++;
      }

      free__MatrixLib__(Distance);
    }

  }

  // If the Triburary_Nodes chain lenght is less than 3 assign al the node
  if(NumTributaryNodes < Ndim + 1)
  {
    fprintf(stderr,"%s %i : %s -> %i\n",
      "Warning in tributary__aLME__ for particle",Indx_p,
      "Insufficient nodal connectivity",NumTributaryNodes);
    exit(EXIT_FAILURE);
  }

  free(Array_Nodes0);

  return Triburary_Nodes;
}

/****************************************************************************/
