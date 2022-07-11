// clang-format off
#include <petscsystypes.h>
#include <petscvec.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "Macros.h"
#include "Types.h"
#include "Globals.h"
#include "Matlib.h"
#include "Particles.h"
#include "Nodes/LME.h"
#include <petsctao.h>
// clang-format on

/****************************************************************************/

// Auxiliar functions to compute the shape functions

typedef struct LME_ctx *LME_ctx_ptr;
struct LME_ctx {
  PetscInt N_a;
  PetscScalar beta;
  Vec l_a;
  Vec p_a;
};

/* -------------- User-defined routines ---------- */
static PetscErrorCode __function_gradient(Tao tao, Vec lambda,
                                           PetscScalar *log_Z, Vec r,
                                           void *logZ_ctx);

static PetscErrorCode __hessian(Tao tao, Vec lambda, Mat H, Mat Hpre,
                                  void *logZ_ctx);

static PetscScalar __eval_f_a(const PetscScalar *la, const PetscScalar *lambda,
                              PetscScalar Beta);

static void initialise_lambda__LME__(int, Matrix, Matrix, Matrix, double);

static ChainPtr __tributary_nodes(int, Matrix, double, int, Mesh);

/****************************************************************************/

void initialize__LME__(Particle MPM_Mesh, Mesh FEM_Mesh) {

  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP; // Number of gauss-points in the simulation
  unsigned Nelem = FEM_Mesh.NumElemMesh; // Number of elements
  int I0;                                // Closest node to the particle
  unsigned p;
  bool Is_particle_reachable;

  ChainPtr Locality_I0; // List of nodes close to the node I0_p
  Matrix lambda_p;      // Lagrange multiplier
  double Beta_p;        // Thermalization or regularization parameter

#pragma omp parallel shared(Np, Nelem)
  {

#pragma omp for private(p, Is_particle_reachable, lambda_p, Beta_p)
    for (p = 0; p < Np; p++) {

      Is_particle_reachable = false;
      unsigned idx_element = 0;

      // Particle position
      Matrix X_p =
          memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);

      // Loop over the element mesh
      while ((Is_particle_reachable == false) && (idx_element < Nelem)) {

        /* Get the element properties */
        ChainPtr Elem_p_Connectivity = FEM_Mesh.Connectivity[idx_element];
        Matrix Elem_p_Coordinates = get_nodes_coordinates__MeshTools__(
            Elem_p_Connectivity, FEM_Mesh.Coordinates);

        /* Check out if the GP is in the Element */
        if (FEM_Mesh.In_Out_Element(X_p, Elem_p_Coordinates) == true) {

          // Particle will be initilise
          Is_particle_reachable = true;

          // Assign the index of the element
          MPM_Mesh.Element_p[p] = idx_element;

          // Asign to each particle the closest node in the mesh and to this
          // node asign the particle
          MPM_Mesh.I0[p] = get_closest_node__MeshTools__(
              X_p, Elem_p_Connectivity, FEM_Mesh.Coordinates);

          // Initialize Beta
          Beta_p = beta__LME__(gamma_LME, FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);

          // Initialise lambda for the Nelder-Mead using Bo-Li approach
          if (strcmp(wrapper_LME, "Nelder-Mead") == 0) {
            initialise_lambda__LME__(p, X_p, Elem_p_Coordinates, lambda_p,
                                     Beta_p);
          }
        }

        /* Free coordinates of the element */
        free__MatrixLib__(Elem_p_Coordinates);

        ++idx_element;
      }

      if (!Is_particle_reachable) {
        fprintf(stderr, "%s : %s %i\n", "Error in initialize__LME__()",
                "The search algorithm was unable to find particle", p);
        exit(EXIT_FAILURE);
      }
    }

#pragma omp barrier

    /*
      Activate the nodes near the particles
    */
    for (p = 0; p < Np; p++) {

      I0 = MPM_Mesh.I0[p];

      Locality_I0 = FEM_Mesh.NodalLocality_0[I0];

      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
        push__SetLib__(&FEM_Mesh.List_Particles_Node[I0], p);
      }

      //    FEM_Mesh.Num_Particles_Node[I0] += 1;

      while (Locality_I0 != NULL) {
        if (FEM_Mesh.ActiveNode[Locality_I0->Idx] == false) {
          FEM_Mesh.ActiveNode[Locality_I0->Idx] = true;
        }

        Locality_I0 = Locality_I0->next;
      }
    }

#pragma omp for private(p, lambda_p, Beta_p)
    for (p = 0; p < Np; p++) {

      // Get some properties for each particle
      Matrix X_p =
          memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);
      lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[p]);
      Beta_p = MPM_Mesh.Beta.nV[p];

      // Get the initial connectivity of the particle
      MPM_Mesh.ListNodes[p] =
          __tributary_nodes(p, X_p, Beta_p, MPM_Mesh.I0[p], FEM_Mesh);

      // Calculate number of nodes
      MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

      // Update the value of the thermalization parameter
      Beta_p = beta__LME__(gamma_LME, FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);
      MPM_Mesh.Beta.nV[p] = Beta_p;
    }
  }
}

/****************************************************************************/

/**
 * @brief Get the thermalization parameter beta using the global variable
 * gamma_LME.
 *
 * @param Gamma User define parameter to control the value of the thermalization
 * parameter.
 * @param h_avg Average mesh size
 * @return double
 */
double beta__LME__(double Gamma, double h_avg) {
  return Gamma / (h_avg * h_avg);
}

/****************************************************************************/

static void initialise_lambda__LME__(int Idx_particle, Matrix X_p,
                                     Matrix Elem_p_Coordinates, //
                                     Matrix lambda, // Lagrange multiplier.
                                     double Beta)   // Thermalization parameter.
{

  int Ndim = NumberDimensions;
  int Nnodes_simplex = Ndim + 1;
  int Size_element = Elem_p_Coordinates.N_rows;
  double sqr_dist_i;

  int *simplex;

  Matrix Norm_l = allocZ__MatrixLib__(Size_element, 1);
  Matrix l = allocZ__MatrixLib__(Size_element, Ndim);

  Matrix A = allocZ__MatrixLib__(Ndim, Ndim);
  Matrix b = allocZ__MatrixLib__(Ndim, 1);
  Matrix x;

  // Initialise a list with distances and order
  for (int i = 0; i < Size_element; i++) {

    sqr_dist_i = 0.0;

    for (int j = 0; j < Ndim; j++) {
      l.nM[i][j] = X_p.nV[i] - Elem_p_Coordinates.nM[i][j];
      sqr_dist_i += DSQR(l.nM[i][j]);
    }

    Norm_l.nV[i] = sqr_dist_i;
  }

  if (Size_element == 3) {
    simplex = (int *)Allocate_ArrayZ(Nnodes_simplex, sizeof(int));
    simplex[0] = 0;
    simplex[1] = 1;
    simplex[2] = 2;
  } else if (Size_element == 4) {
    simplex = (int *)Allocate_ArrayZ(Nnodes_simplex, sizeof(int));
    simplex[0] = 0;
    simplex[1] = 1;
    simplex[2] = 2;
  } else {
    exit(0);
  }

  // Assemble matrix to solve the system Ax = b
  for (int i = 1; i < Nnodes_simplex; i++) {

    b.nV[i - 1] = -Beta * (Norm_l.nV[simplex[0]] - Norm_l.nV[simplex[i]]);

    for (int j = 0; j < Ndim; j++) {
      A.nM[i - 1][j] = l.nM[simplex[i]][j] - l.nM[simplex[0]][j];
    }
  }

  // Solve the system
  if (rcond__TensorLib__(A.nV) < 1E-8) {
    fprintf(stderr, "%s %i : %s \n",
            "Error in initialise_lambda__LME__ for particle", Idx_particle,
            "The Hessian near to singular matrix!");
    exit(EXIT_FAILURE);
  }

  x = solve__MatrixLib__(A, b);

  // Update the value of lambda
  for (int i = 0; i < Ndim; i++) {
    lambda.nV[i] = x.nV[i];
  }

  // Free memory
  free(simplex);
  free__MatrixLib__(Norm_l);
  free__MatrixLib__(l);
  free__MatrixLib__(A);
  free__MatrixLib__(b);
  free__MatrixLib__(x);
}

/****************************************************************************/

static double * p__LME__(const double *l_a,
                               double *lambda_tr, double Beta, unsigned N_a) {

  double *p_a = (double *)calloc(N_a, sizeof(double));

  /* Definition of some parameters */
  LME_ctx_ptr user_ptr;
  PetscErrorCode STATUS = EXIT_SUCCESS;
  PetscInt MaxIter = max_iter_LME;
  PetscInt Ndim = NumberDimensions;
  PetscInt NumIter;
  PetscScalar zero = 0.0;
  Tao tao;
  Mat H;
  Vec l_aux;
  Vec p_aux;
  Vec lambda_aux;

  /* Create user-defined variable */
  user_ptr->beta = Beta;
  user_ptr->l_a = l_aux;
  user_ptr->N_a = N_a;
  user_ptr->p_a = p_aux;

  /* Create Hessian */
  PetscCall(MatCreate(PETSC_COMM_SELF, &H));
  PetscCall(MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, Ndim, Ndim));

  /* Create TAO solver with desired solution method */
  PetscCall(TaoCreate(PETSC_COMM_SELF, &tao));
  PetscCall(TaoSetType(tao, TAONTL));

  /* Set solution vec */
  PetscCall(TaoSetSolution(tao, lambda_aux));

  /* Set routines for function, gradient, hessian evaluation */
  PetscCall(
      TaoSetObjectiveAndGradient(tao, NULL, __function_gradient, user_ptr));
  PetscCall(TaoSetHessian(tao, H, H, __hessian, user_ptr));

  /* Solve the system */
  PetscCall(TaoSolve(tao));

  /* Destroy auxiliar variables */
  PetscCall(TaoDestroy(&tao));
  PetscCall(MatDestroy(&H));

  return p_a;
}

/****************************************************************************/

static PetscErrorCode dp__LME__(double *dp_a, const double *l_a,
                                double *lambda_tr, double Beta, int N_a) {

  /* Definition of some parameters */
  LME_ctx_ptr user_ptr;
  PetscErrorCode STATUS = EXIT_SUCCESS;
  PetscInt MaxIter = max_iter_LME;
  PetscInt Ndim = NumberDimensions;
  PetscInt NumIter;
  PetscScalar zero = 0.0;
  Tao tao;
  Mat H;
  Vec l_aux;
  Vec p_aux;
  Vec lambda_aux;

  /* Create user-defined variable */
  user_ptr->beta = Beta;
  user_ptr->l_a = l_aux;
  user_ptr->N_a = N_a;
  user_ptr->p_a = p_aux;

  /* Create Hessian */
  PetscCall(MatCreate(PETSC_COMM_WORLD, &H));
  PetscCall(MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, Ndim, Ndim));

  /* Create TAO solver with desired solution method */
  PetscCall(TaoCreate(PETSC_COMM_SELF, &tao));
  PetscCall(TaoSetType(tao, TAONTL));

  /* Set solution vec */
  PetscCall(TaoSetSolution(tao, lambda_aux));

  /* Set routines for function, gradient, hessian evaluation */
  PetscCall(
      TaoSetObjectiveAndGradient(tao, NULL, __function_gradient, user_ptr));
  PetscCall(TaoSetHessian(tao, H, H, __hessian, user_ptr));

  /* Solve the system */
  PetscCall(TaoSolve(tao));

  /* Inverse of the Hessian */
  Jm1 = inverse__MatrixLib__(J);

  /*
    Fill the gradient for each node
  */
  for (int a = 0; a < N_a; a++) {
    la.nV = l.nM[a];
    Jm1_la = matrix_product__MatrixLib__(Jm1, la);

    for (int i = 0; i < Ndim; i++) {
      dp_a[a * Ndim + i] = -p.nV[a] * Jm1_la.nV[i];
    }

    free__MatrixLib__(Jm1_la);
  }

  PetscCall(TaoDestroy(&tao));
  PetscCall(MatDestroy(&H));

  return STATUS;
}

/****************************************************************************/

/**
 * @brief Evaluates the function and corresponding gradient.
 *
 * @param tao the Tao context
 * @param lambda
 * @param log_Z
 * @param r
 * @param logZ_ctx
 * @return PetscErrorCode
 */
static PetscErrorCode __function_gradient(Tao tao, Vec lambda,
                                           PetscScalar *log_Z, Vec r,
                                           void *logZ_ctx) {
  /* Definition of some parameters */
  LME_ctx_ptr user_ptr = (LME_ctx_ptr)logZ_ctx;

  /* Get constants */
  PetscInt Ndim = NumberDimensions;
  PetscScalar beta = user_ptr->beta;
  PetscInt N_a = user_ptr->N_a;

  /* Read auxiliar variables */
  PetscScalar *p_a_ptr;
  VecGetArray(user_ptr->p_a, &p_a_ptr);
  PetscScalar *r_ptr;
  VecGetArray(r, &r_ptr);
  const PetscScalar *l_a_ptr;
  VecGetArrayRead(user_ptr->l_a, &l_a_ptr);

  /* Compute partition function Z and the shape function p_a */
  PetscScalar Z = 0;
  for (PetscInt a = 0; a < N_a; a++) {
    p_a_ptr[a] = exp(__eval_f_a(&l_a_ptr[a * Ndim], lambda, beta));
    Z += p.nV[a];
  }

  /* Divide by Z and get the final value of the shape function */
  PetscScalar Z_m1 = 1.0 / Z;
  for (PetscInt a = 0; a < N_a; a++) {
    p_a_ptr[a] *= Z_m1;
  }

  /* Evaluate objective function */
  *log_Z = log(Z);

  /* Compute the gradient */
  for (int i = 0; i < Ndim; i++) {
    r_ptr[i] = 0.0 for (int a = 0; a < N_a; a++) {
      r_ptr[i] += p_a_ptr[a] * l_a_ptr[a * Ndim + i];
    }
  }

  /* Restore auxiliar pointers */
  VecRestoreArray(user_ptr->p_a, &p_a_ptr);
  VecRestoreArray(r, &r_ptr);
  VecRestoreArrayRead(user_ptr->l_a, &l_a_ptr);

  return logZ;
}

/****************************************************************************/

/**
 * @brief
 *
 * @param tao
 * @param lambda
 * @param H
 * @param Hpre
 * @param logZ_ctx
 * @return PetscErrorCode
 */
static PetscErrorCode __hessian(Tao tao, Vec lambda, Mat H, Mat Hpre,
                                  void *logZ_ctx) {

  /* Definition of some parameters */
  LME_ctx_ptr user_ptr = (LME_ctx_ptr)logZ_ctx;

  /* Get constants */
  PetscInt Ndim = NumberDimensions;
  PetscInt N_a = user_ptr->N_a;
  PetscScalar H_a;

  /* Read auxiliar variables */
  PetscScalar *p_a_ptr;
  VecGetArray(user_ptr->p_a, &p_a_ptr);
  const PetscScalar *l_a_ptr;
  VecGetArrayRead(user_ptr->l_a, &l_a_ptr);

  Vec r;
  TaoGetGradient(tao, &r, NULL, NULL);
  const PetscScalar *r_ptr;
  VecGetArrayRead(r, &r_ptr);

  /* Set to zero the Hessian */
  PetscCall(MatZeroEntries(H));

  /* Fill the Hessian */
  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      /* First component */
      for (int a = 0; a < N_a; a++) {
        H_a = p_a_ptr[a] * l_a_ptr[a * Ndim + i] * l_a_ptr[a * Ndim + j];
        MatSetValues(H, 1, i, 1, j, H_a, ADD_VALUES);
      }
      /* Second component */
      H_a = -r_ptr[i] * r_ptr[j];
      MatSetValues(H, 1, i, 1, j, H_a, ADD_VALUES);
    }
  }

  /* Restore auxiliar pointers */
  VecRestoreArray(user_ptr->p_a, &p_a_ptr);
  VecRestoreArrayRead(user_ptr->l_a, &l_a_ptr);
  VecRestoreArrayRead(r, &r_ptr);

  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
  if (H != Hpre) {
    MatAssemblyBegin(Hpre, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Hpre, MAT_FINAL_ASSEMBLY);
  }
  return 0;
}

/****************************************************************************/

/**
 * @brief fa (scalar): the function fa that appear in [1].
 *
 * @param la Vector form node "a" to particle.
 * @param lambda Lagrange multiplier.
 * @param Beta Thermalization parameter.
 * @return PetscScalar
 */
static PetscScalar __eval_f_a(const PetscScalar *la, const PetscScalar *lambda,
                              PetscScalar Beta) {
  PetscInt Ndim = NumberDimensions;
  PetscScalar la_x_la = 0.0;
  PetscScalar la_x_lambda = 0.0;

  for (PetscInt i = 0; i < Ndim; i++) {
    la_x_la += la[i] * la[i];
    la_x_lambda += la[i] * lambda[i];
  }

  PetscScalar fa = -Beta * la_x_la + la_x_lambda;

  return fa;
}

/****************************************************************************/

int local_search__LME__(Particle MPM_Mesh, Mesh FEM_Mesh)
/*
  Search the closest node to the particle based in its previous position.
*/
{
  int STATUS = EXIT_SUCCESS;
  int STATUS_p = EXIT_SUCCESS;
  int Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned p;

  Matrix X_p;           // Particle position
  Matrix dis_p;         // Particle displacement
  Matrix lambda_p;      // Lagrange multiplier of the shape function
  double Beta_p;        // Themalization parameter
  ChainPtr Locality_I0; // List of nodes close to the node I0_p

#pragma omp parallel shared(Np)
  {

#pragma omp for private(p, X_p, dis_p, Locality_I0)
    for (p = 0; p < Np; p++) {

      // Get the global coordinates and displacement of the particle
      X_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);
      dis_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.dis.nM[p]);

      // Check if the particle is static or is in movement
      if (norm__MatrixLib__(dis_p, 2) > 0) {

        // Update the index of the closest node to the particle
        Locality_I0 = FEM_Mesh.NodalLocality_0[MPM_Mesh.I0[p]];
        MPM_Mesh.I0[p] = get_closest_node__MeshTools__(X_p, Locality_I0,
                                                       FEM_Mesh.Coordinates);

        // Search particle in the sourrounding elements to this node
        MPM_Mesh.Element_p[p] =
            search_particle_in_surrounding_elements__Particles__(
                p, X_p, FEM_Mesh.NodeNeighbour[MPM_Mesh.I0[p]], FEM_Mesh);
        if (MPM_Mesh.Element_p[p] == -999) {
          fprintf(stderr,
                  "" RED " Error in " RESET "" BOLDRED
                  "search_particle_in_surrounding_elements__Particles__(%i,,)"
                  " " RESET " \n",
                  p);
          STATUS = EXIT_FAILURE;
        }
      }
    }

#pragma omp barrier

    // Activate the nodes near the particle
    for (p = 0; p < Np; p++) {

      int I0 = MPM_Mesh.I0[p];

      Locality_I0 = FEM_Mesh.NodalLocality_0[MPM_Mesh.I0[p]];
      while (Locality_I0 != NULL) {
        if (FEM_Mesh.ActiveNode[Locality_I0->Idx] == false) {
          FEM_Mesh.ActiveNode[Locality_I0->Idx] = true;
        }

        Locality_I0 = Locality_I0->next;
      }

      if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
        push__SetLib__(&FEM_Mesh.List_Particles_Node[I0], p);
      }
    }

    // Update the shape function

#pragma omp for private(p, X_p, lambda_p, Beta_p)
    for (p = 0; p < Np; p++) {

      lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[p]);
      Beta_p = MPM_Mesh.Beta.nV[p]; // Thermalization parameter

      //  Get the global coordinates of the particle
      X_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);

      //  Free previous list of tributary nodes to the particle
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      MPM_Mesh.ListNodes[p] = NULL;

      //  Calculate the new connectivity with the previous value of beta
      MPM_Mesh.ListNodes[p] =
          __tributary_nodes(p, X_p, Beta_p, MPM_Mesh.I0[p], FEM_Mesh);

      //  Calculate number of nodes
      MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

      //  Compute the thermalization parameter for the new set of nodes
      //  and update it
      Beta_p = beta__LME__(gamma_LME, FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);
      MPM_Mesh.Beta.nV[p] = Beta_p;

    }
  }

  if (STATUS == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/****************************************************************************/

/**
 * @brief Compute a set with the sourrounding nodes of the particle
 * 
 * @param Indx_p 
 * @param X_p Coordinates of the particle
 * @param Beta_p Thermalization parameter of the particle
 * @param I0 Index of the closest node to the particle
 * @param FEM_Mesh Variable wih information of the background set of nodes
 * @return ChainPtr 
 */
static ChainPtr __tributary_nodes(int Indx_p, Matrix X_p, double Beta_p, int I0,
                                 Mesh FEM_Mesh)
{

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;
  /* Number of dimensionws of the problem */
  int Ndim = NumberDimensions;

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = memory_to_matrix__MatrixLib__(Ndim, 1, NULL);
  Matrix Metric = Identity__MatrixLib__(Ndim);
  ChainPtr Set_Nodes0 = NULL;
  int *Array_Nodes0;
  int NumNodes0;
  int Node0;

  /* Counter */
  int NumTributaryNodes = 0;

  /* Get the search radius */
  double Ra = sqrt(-log(TOL_zero_LME) / Beta_p);

  /*
    Get nodes close to the particle
  */
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0, NumNodes0);

  /* Loop over the chain with the tributary nodes */
  for (int i = 0; i < NumNodes0; i++) {

    Node0 = Array_Nodes0[i];

    if (FEM_Mesh.ActiveNode[Node0] == true) {
      /* Assign to a pointer the coordinates of the nodes */
      X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

      /* Get a vector from the GP to the node */
      Distance = substraction__MatrixLib__(X_p, X_I);

      /* If the node is near the GP push in the chain */
      if (generalised_Euclidean_distance__MatrixLib__(Distance, Metric) <= Ra) {
        push__SetLib__(&Triburary_Nodes, Node0);
        NumTributaryNodes++;
      }

      /* Free memory of the distrance vector */
      free__MatrixLib__(Distance);
    }
  }

  /*
    If the Triburary_Nodes chain lenght is less than 3 assign al the node
  */
  if (NumTributaryNodes < Ndim + 1) {
    fprintf(stderr, "%s %i : %s -> %i\n",
            "Warning in __tributary_nodes for particle", Indx_p,
            "Insufficient nodal connectivity", NumTributaryNodes);
    exit(EXIT_FAILURE);
  }

  /* Free memory */
  free(Array_Nodes0);
  free__MatrixLib__(Metric);

  return Triburary_Nodes;
}

/****************************************************************************/
