#include "nl-partsol.h"


// Auxiliar functions to compute the shape functions
static double fa__LME__(Matrix, Matrix, Matrix, double);
static double logZ__LME__(Matrix, Matrix, Matrix, double);
static Matrix r__LME__(Matrix, Matrix);
static Matrix J__LME__(Matrix, Matrix, Matrix);

// Auxiliar functions for the Neldel Mead in the LME
static void initialise_lambda__LME__(int,Matrix,Matrix,Matrix,double, double);
static Matrix gravity_center_Nelder_Mead__LME__(Matrix);
static void order_logZ_simplex_Nelder_Mead__LME__(Matrix, Matrix);
static void expansion_Nelder_Mead__LME__(Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,double,double);
static void contraction_Nelder_Mead__LME__(Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,double,double);
static void shrinkage_Nelder_Mead__LME__(Matrix,Matrix,Matrix,Matrix,double);

// Call global var¡ables
double gamma_LME;
double curvature_LME;
double TOL_LME;
int max_iter_LME;

// Nelder-Mead parameters
double NM_rho_LME = 1.0;
double NM_chi_LME = 2.0;
double NM_gamma_LME = 0.5;
double NM_sigma_LME = 0.5;
double tau_NM = 1E-3;

/****************************************************************************/

void initialize__LME__(
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

  Matrix Metric_p; // Define a metric tensor
  Matrix X_p; // Particle coordinates  
  Matrix Delta_Xip; // Distance from GP to the nodes
  Matrix lambda_p; // Lagrange multiplier
  Tensor F_p; // Particle deformation gradient, only for anysotropic
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
    F_p = memory_to_tensor__TensorLib__(MPM_Mesh.Phi.F_n.nM[p],2);
    Metric_p = metric__LME__(F_p);

    /*
      Loop over the element mesh
    */
    for(int i = 0 ; i<Nelem ; i++)
    {

      /* Get the element properties */
      Elem_p_Connectivity = FEM_Mesh.Connectivity[i];
      Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(Elem_p_Connectivity, FEM_Mesh.Coordinates);

      /* Check out if the GP is in the Element */
      if(FEM_Mesh.In_Out_Element(X_p,Elem_p_Coordinates))
      {

        /* Asign to each particle the closest node in the mesh
          and to this node asign the particle */
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p,Elem_p_Connectivity,FEM_Mesh.Coordinates);

        /* Calculate distance from particle to each node in the neibourhood */
        MPM_Mesh.ListNodes[p] = copy__SetLib__(Elem_p_Connectivity);
        Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

        /* Initialize Beta */
        Beta_p = beta__LME__(Delta_Xip, gamma_LME, FEM_Mesh.DeltaX);

        /* Initialise lambda using Bo-Li approach */
        initialise_lambda__LME__(p, X_p, Elem_p_Coordinates, lambda_p, Beta_p, FEM_Mesh.DeltaX);

        /* Free memory */ 
        free__MatrixLib__(Elem_p_Coordinates);
        free__MatrixLib__(Delta_Xip);

        /* Get the initial connectivity of the particle */
        MPM_Mesh.ListNodes[p] = tributary__LME__(X_p,Metric_p,Beta_p,MPM_Mesh.I0[p],FEM_Mesh);

        /* Measure the size of the connectivity */
        MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

        /* Active those nodes that interact with the particle */
        asign_to_nodes__Particles__(p, MPM_Mesh.ListNodes[p], FEM_Mesh);
       	
        /* Calculate distance from particle to each node in the neibourhood */
        Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p],X_p,FEM_Mesh.Coordinates);

        /* Update the value of the thermalization parameter */
        Beta_p = beta__LME__(Delta_Xip, gamma_LME, FEM_Mesh.DeltaX);
        MPM_Mesh.Beta.nV[p] = Beta_p;

        /* Update lagrange multiplier with Newton-Rapson or with Nelder-Mead */
        MPM_Mesh.update_lambda(p, Delta_Xip, lambda_p, Metric_p, Beta_p);

        /* 
          Free memory
        */
        free__MatrixLib__(Delta_Xip);

        break;
      }      

      /* 
        Free memory
      */
      free__MatrixLib__(Elem_p_Coordinates);

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

 Matrix metric__LME__(Tensor F)
 /*!
   Return the metric tensor intrucing curvature as a convex combination of the 
   right Cauch-Green tensor (C = F^{T}F) and the identiy (Euclidean norm).
 */
 {
    int Ndim = NumberDimensions;
    double C_ij;
    Matrix Metric = allocZ__MatrixLib__(Ndim,Ndim);

    for(int i = 0 ; i < Ndim ; i++)
    {

      // Introduce Euclidean metric contribution
      Metric.nM[i][i] += (1.0 - curvature_LME);

      // Include include non-Euclidean metric
      if(curvature_LME > 0.0)
      {

        for(int j = 0 ; j < Ndim ; j++)
        {
          C_ij = 0.0;

          for(int k = 0 ; k < Ndim ; k++)
          {
            C_ij += F.N[k][i]*F.N[k][j];
          }

          Metric.nM[i][j] += curvature_LME*C_ij;
        
        }
      }
    } 

   return Metric;
 }


/****************************************************************************/

static void initialise_lambda__LME__(
  int Idx_particle,
  Matrix X_p,
  Matrix Elem_p_Coordinates, //
  Matrix lambda, // Lagrange multiplier.
  double Beta,
  double DeltaX) // Thermalization parameter.
{

  int Ndim = NumberDimensions;
/*  int Nnodes_simplex = Ndim + 1;
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

  if(Size_element == 8)
  {
    simplex = (int *)Allocate_ArrayZ(Nnodes_simplex,sizeof(int));
    simplex[0] = 0;
    simplex[1] = 1;
    simplex[2] = 3;
    simplex[3] = 4;
  }

  // Assemble matrix to solve the system Ax = b
  for(int i = 1 ; i<Nnodes_simplex ; i++)
  {

    b.nV[i-1] = - Beta*(Norm_l.nV[simplex[0]] - Norm_l.nV[simplex[i]]);

    for(int j = 0 ; j<Ndim ; j++)
    {
      A.nM[i-1][j] = l.nM[simplex[i]][j] - l.nM[simplex[0]][j];
    }
  }


  if(fabs(I3__MatrixLib__(A)) < TOL_zero)
  {
    strcpy(A.Info,"A");
    print__MatrixLib__(A,Ndim,Ndim);
    fprintf(stderr,"%s %i : %s \n",
      "Error in initialise_lambda__LME__ for particle",
      Idx_particle,"Determinant of A is null !");
    exit(EXIT_FAILURE);
  }

  // Solve the system
  x = Solve_Linear_Sistem(A,b);
*/
  // Update the value of lambda
  for(int i = 0 ; i<Ndim ; i++)
  {
//    lambda.nV[i] = x.nV[i];
    lambda.nV[i] = 1.0;
  }

/*  // Free memory
  free(simplex);
  free__MatrixLib__(Norm_l);
  free__MatrixLib__(l);
  free__MatrixLib__(A);
  free__MatrixLib__(b);
  free__MatrixLib__(x);*/
}

/****************************************************************************/

void update_lambda_Newton_Rapson__LME__(
  int Idx_particle,
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix lambda, // Lagrange multiplier.
  Matrix Metric, // Measure for the norm definition.
  double Beta) // Thermalization parameter.
/*!
  Get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method is the Newton-Rapson.
*/
{
  /*
    Definition of some parameters
  */
  int MaxIter = max_iter_LME;
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

      if(fabs(I3__MatrixLib__(J)) < TOL_zero)
      {
        fprintf(stderr,"%s %i : %s \n",
          "Error in lambda_Newton_Rapson__LME__ for particle",
          Idx_particle,"The Hessian is singular !");
        exit(EXIT_FAILURE);
      }
    
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
    fprintf(stderr,"%s %i : %s (%i)\n",
      "Warning in lambda_Newton_Rapson__LME__ for particle",Idx_particle,
      "No convergence reached in the maximum number of interations",MaxIter);
    fprintf(stderr,"%s : %e\n", "Total Error",norm_r);
    exit(EXIT_FAILURE);
  }
  
}

/****************************************************************************/

void update_lambda_Nelder_Mead__LME__(
  int Idx_particle,
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix lambda, // Lagrange multiplier.
  Matrix Metric, // Measure for the norm definition.
  double Beta) // Thermalization parameter.
/*!
  Get the lagrange multipliers "lambda" (1 x dim) for the LME 
  shape function. The numerical method is the Nelder-Mead.
  In this method, each vertex is represented by a lagrange multiplier
*/
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

  
  for(int a = 0 ; a<Nnodes_simplex ; a++)
  {    
    // Compute the initial values of logZ in each vertex of the simplex
    simplex_a.nV = simplex.nM[a];
    logZ.nV[a] = logZ__LME__(l, simplex_a, Metric, Beta);
  }

  // Nelder-Mead main loop
  while(NumIter <= MaxIter)
  {

    order_logZ_simplex_Nelder_Mead__LME__(logZ, simplex);

    logZ_0 = logZ.nV[0];
    logZ_n = logZ.nV[Nnodes_simplex-2];
    logZ_n1 = logZ.nV[Nnodes_simplex-1];

    // Check convergence
    if(fabs(logZ_0 - logZ_n1) > TOL_LME)
    {

      // Spin the simplex to get the simplex with the smallest normalized volume
      // spin_Nelder_Mead__LME__(simplex);

      // Compute the gravity center of the simplex
      gravity_center = gravity_center_Nelder_Mead__LME__(simplex);

      // Compute the reflected point and evaluate the objetive function in this point
      reflected_point = allocZ__MatrixLib__(1, Ndim);

      for(int i = 0 ; i<Ndim ; i++)
      {
        reflected_point.nV[i] = gravity_center.nV[i] + NM_rho_LME*(gravity_center.nV[i] - simplex.nM[Nnodes_simplex-1][i]);
      }

      logZ_reflected_point = logZ__LME__(l,reflected_point,Metric,Beta);

      // Do an expansion using the reflected point
      if(logZ_reflected_point < logZ_0)
      {
        expansion_Nelder_Mead__LME__(simplex,logZ,reflected_point,gravity_center,l,Metric,Beta,logZ_reflected_point);
      }
      // Take the reflected point
      else if((logZ_reflected_point > logZ_0) && (logZ_reflected_point < logZ_n))
      {
        for(int i = 0 ; i<Ndim ; i++)
        {
          simplex.nM[Nnodes_simplex-1][i] = reflected_point.nV[i];
        }

        logZ.nV[Nnodes_simplex-1] = logZ_reflected_point; 

      }
      // Do a contraction using the reflected point (or a shrinkage)
      else if(logZ_reflected_point >= logZ_n)
      {
        contraction_Nelder_Mead__LME__(simplex,logZ,reflected_point,gravity_center,l,Metric,Beta,logZ_reflected_point);
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
      "Warning in lambda_Nelder_Mead__LME__ for particle",Idx_particle,
      "No convergence reached in the maximum number of interations",MaxIter);
    fprintf(stderr,"%s : %e\n", "Total Error",fabs(logZ_0 - logZ_n1));
    exit(EXIT_FAILURE);
  }

  // Update the value of lambda
  for(int i = 0 ; i<Ndim ; i++)
  {
    lambda.nV[i] = simplex.nM[0][i];
  }

  /*
    Free memory
  */
  free__MatrixLib__(simplex);
  free__MatrixLib__(logZ);

}


/****************************************************************************/

static void order_logZ_simplex_Nelder_Mead__LME__(
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

/*static double spin_Nelder_Mead__LME__(
  Matrix simplex)
{

    int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;

  Matrix L_simplex = allocZ__MatrixLib__();


[n,sp]=size(LL);
mat=ones(n);

mat(:,1:2)=LL;
DT=abs(det(mat)/factorial(sp));

dist(n,n)=0;
for i=1:n
  for k=1:n
    for j=1:sp
      dist(i,k)=dist(i,k)+(LL(i,j)-LL(k,j))^2;
    end
    dist(i,k)=sqrt(dist(i,k));
  end
end

diam=max(max(dist));

Vn=DT/diam^sp;


[~,sp]=size(LL);

DIS=max(diam/50,eps*2);


LL1(1,:)=LL(1,:);
LL2(1,:)=LL(1,:);

LL1(2,1)=LL(2,1)+DIS;
LL1(2,2)=LL(2,2);
LL1(3,1)=LL(3,1);
LL1(3,2)=LL(3,2)-DIS;

LL2(2,1)=LL(2,1)-DIS;
LL2(2,2)=LL(2,2);
LL2(3,2)=LL(3,2)+DIS;

f_1(1)=1e10;
f_2(1)=1e10;

[f_1(2),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL1(2,:),near);
[f_1(3),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL1(3,:),near);
[f_2(2),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL2(2,:),near);
[f_2(3),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL2(3,:),near);

f1=min(f_1);
f2=min(f_2);

if f1<f2
  LL2=LL1;
  f_2=f_1;
end

f_2(1),~,~,~]=LME.Gamma_(sp,x_a,x,beta,LL1(1,:),near);

[LL2,f_m,f_mm,f_w]=LME.order(LL2,f_2);

}*/

/****************************************************************************/

static Matrix gravity_center_Nelder_Mead__LME__(Matrix simplex)
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

static void expansion_Nelder_Mead__LME__(
  Matrix simplex,
  Matrix logZ,
  Matrix reflected_point,
  Matrix gravity_center,
  Matrix l,
  Matrix Metric,
  double Beta, 
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
    expanded_point.nV[i] = gravity_center.nV[i] + NM_chi_LME*(reflected_point.nV[i] - gravity_center.nV[i]);
  }

  logZ_expanded_point = logZ__LME__(l, expanded_point, Metric, Beta); 

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

static void contraction_Nelder_Mead__LME__(
  Matrix simplex,
  Matrix logZ,
  Matrix reflected_point,
  Matrix gravity_center,
  Matrix l,
  Matrix Metric,
  double Beta, 
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
      contracted_point.nV[i] = gravity_center.nV[i] + NM_gamma_LME*(reflected_point.nV[i] - gravity_center.nV[i]);
    }

    logZ_contracted_point = logZ__LME__(l, contracted_point, Metric, Beta); 

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
      shrinkage_Nelder_Mead__LME__(simplex,logZ,l,Metric,Beta);
    }
  }
  // Internal contraction
  else if(logZ_reflected_point > logZ_n1)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      contracted_point.nV[i] = gravity_center.nV[i] - NM_gamma_LME*(gravity_center.nV[i] - simplex.nM[Nnodes_simplex-1][i]);
    }

    logZ_contracted_point = logZ__LME__(l, contracted_point, Metric, Beta); 

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
      shrinkage_Nelder_Mead__LME__(simplex,logZ,l,Metric,Beta);
    }
  }

  free__MatrixLib__(contracted_point);

}

/****************************************************************************/

static void shrinkage_Nelder_Mead__LME__(
  Matrix simplex,
  Matrix logZ,
  Matrix l,
  Matrix Metric,
  double Beta)
{
  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;

  // Axiliar function to get the coordinates of the simplex in the node a
  Matrix simplex_a = memory_to_matrix__MatrixLib__(1,Ndim,NULL);

  for (int a = 0 ; a<Nnodes_simplex ; a++)
  {
    for(int i = 0 ; i<Ndim ; i++)
    {
      simplex.nM[a][i] = simplex.nM[0][i] + NM_sigma_LME*(simplex.nM[a][i] - simplex.nM[0][i]);
    }
    simplex_a.nV = simplex.nM[a];
    logZ.nV[a] = logZ__LME__(l, simplex_a, Metric, Beta);
  }
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
  Matrix p = allocZ__MatrixLib__(1,N_a); // Shape function in the nodes
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

static double logZ__LME__(
  Matrix l, // Set than contanins vector form neighborhood nodes to particle.
  Matrix lambda, // Lagrange multiplier.
  Matrix Metric, // Measure for the norm definition.
  double Beta) // Thermalization parameter.
{
  /* Definition of some parameters */
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double logZ = 0;
  Matrix la = memory_to_matrix__MatrixLib__(1,Ndim,NULL); // Vector form node "a" to particle.

  for(int a = 0 ; a<N_a ; a++)
  {
    la.nV = l.nM[a];
    Z += exp(fa__LME__(la,lambda,Metric,Beta));
  }

  logZ = log(Z);

  return logZ;
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
