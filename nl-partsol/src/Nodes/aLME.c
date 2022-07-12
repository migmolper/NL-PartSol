// clang-format off
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
#include "Nodes/aLME.h"
// clang-format on

/****************************************************************************/

// Auxiliar functions to compute the shape functions
static double fa__aLME__(Matrix, Matrix, Matrix);
static double logZ__aLME__(Matrix, Matrix, Matrix);
static Matrix r__aLME__(Matrix, Matrix);
static Matrix J__aLME__(Matrix, Matrix, Matrix);
static void initialize_beta__aLME__(Matrix, double, double);
static void update_beta__aLME__(Matrix, Matrix);
static void initialize_Cut_off_Ellipsoid__aLME__(Matrix, double, double);
static void update_cut_off_ellipsoid__aLME__(Matrix, Matrix);
static void initialise_lambda__aLME__(int, Matrix, Matrix, Matrix, Matrix);
static void update_lambda_Newton_Rapson__aLME__(int, Matrix, Matrix, Matrix);
static ChainPtr tributary__aLME__(int, Matrix, Matrix, int, Mesh);

/****************************************************************************/

void initialize__aLME__(Particle MPM_Mesh, Mesh FEM_Mesh) {

  unsigned Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP; // Number of gauss-points in the simulation
  unsigned Nelem = FEM_Mesh.NumElemMesh; // Number of elements
  int I0;                                // Closest node to the particle
  unsigned p;
  bool Is_particle_reachable;

  ChainPtr Locality_I0; // List of nodes close to the node I0_p
  Matrix Delta_Xip;     // Distance from GP to the nodes
  Matrix lambda_p;      // Lagrange multiplier

#pragma omp parallel shared(Np, Nelem)
  {

#pragma omp for private(p, Is_particle_reachable, lambda_p)
    for (p = 0; p < Np; p++) {

      Is_particle_reachable = false;
      unsigned idx_element = 0;

      // Particle position
      Matrix X_p =
          memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);
      Matrix Beta_p =
          memory_to_matrix__MatrixLib__(Ndim, Ndim, MPM_Mesh.Beta.nM[p]);

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
          initialize_beta__aLME__(Beta_p, gamma_LME,
                                  FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);
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

#pragma omp for private(p, Delta_Xip, lambda_p)
    for (p = 0; p < Np; p++) {

      // Get some properties for each particle
      Matrix X_p =
          memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);
      lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[p]);
      Matrix Beta_p =
          memory_to_matrix__MatrixLib__(Ndim, Ndim, MPM_Mesh.Beta.nM[p]);

      // Get the metric tensor and initialize it
      Matrix Cut_off_Ellipsoid = memory_to_matrix__MatrixLib__(
          Ndim, Ndim, MPM_Mesh.Cut_off_Ellipsoid.nM[p]);

      initialize_Cut_off_Ellipsoid__aLME__(Cut_off_Ellipsoid, gamma_LME,
                                           FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);

      // Get the initial connectivity of the particle
      MPM_Mesh.ListNodes[p] = tributary__aLME__(p, X_p, Cut_off_Ellipsoid,
                                                MPM_Mesh.I0[p], FEM_Mesh);

      // Calculate number of nodes
      MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

      // Generate nodal distance list
      Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p], X_p,
                                                FEM_Mesh.Coordinates);

      // Update the value of the thermalization parameter
      initialize_beta__aLME__(Beta_p, gamma_LME,
                              FEM_Mesh.h_avg[MPM_Mesh.I0[p]]);

      update_lambda_Newton_Rapson__aLME__(p, Delta_Xip, lambda_p, Beta_p);

      // Free memory
      free__MatrixLib__(Delta_Xip);
      free(Beta_p.nM);
    }
  }
}

/****************************************************************************/

/**
 * @brief Compute the value of the thermalization parameter using a circular
 * support.
 *
 * @param Beta Termalization tensor
 * @param Gamma Adimensional paramter to control the regularization parameter
 * @param h_avg Minimum size in the all nodal set.
 */
static void initialize_beta__aLME__(Matrix Beta, double Gamma, double h_avg) {
  int Ndim = NumberDimensions;
  double aux = Gamma / (h_avg * h_avg);

  for (int i = 0; i < Ndim; i++) {
    Beta.nM[i][i] = aux;
  }
}

/****************************************************************************/

/**
 * @brief Initializae the definition of the cutt-off matrix
 *
 * @param Cut_off_Ellipsoid
 * @param Gamma User define parameter to control the value of the thermalization
 * parameter
 * @param h_avg Average mesh size
 */
static void initialize_Cut_off_Ellipsoid__aLME__(Matrix Cut_off_Ellipsoid,
                                                 double Gamma, double h_avg) {
  int Ndim = NumberDimensions;
  double aux = Gamma / (-log(TOL_zero_LME) * h_avg * h_avg);

  for (int i = 0; i < Ndim; i++) {
    Cut_off_Ellipsoid.nM[i][i] = aux;
  }
}

/****************************************************************************/

/**
 * @brief
 *
 * @param Idx_particle Index of the particle.
 * @param X_p Particle position
 * @param Elem_p_Coordinates Coordinates of the nodes close to the particle.
 * @param lambda Lagrange multiplier of the particle
 * @param Beta Thermalization matrix.
 */
static void initialise_lambda__aLME__(int Idx_particle, Matrix X_p,
                                      Matrix Elem_p_Coordinates, Matrix lambda,
                                      Matrix Beta) {

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

    b.nV[i - 1] =
        -Beta.nM[i][i] * (Norm_l.nV[simplex[0]] - Norm_l.nV[simplex[i]]);

    for (int j = 0; j < Ndim; j++) {
      A.nM[i - 1][j] = l.nM[simplex[i]][j] - l.nM[simplex[0]][j];
    }
  }

  // Solve the system
  if (rcond__TensorLib__(A.nV) < 1E-8) {
    fprintf(stderr, "%s %i : %s \n",
            "Error in initialise_lambda__aLME__ for particle", Idx_particle,
            "The Hessian near to singular matrix!");
    exit(EXIT_FAILURE);
  }

  x = solve__MatrixLib__(A, b);

  // Update the value of lambda
  for (int i = 0; i < Ndim; i++) {
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

/**
 * @brief unction ot get the lagrange multipliers "lambda" (1 x dim) for the LME
 * shape function. The numerical methodis the Newton-Rapson.
 *
 * @param Idx_particle Current particle
 * @param l Set than contanins vector form neighborhood nodes to particle.
 * @param lambda Lagrange multiplier.
 * @param Beta Thermalization tensor.
 */
static void update_lambda_Newton_Rapson__aLME__(int Idx_particle, Matrix l,
                                                Matrix lambda, Matrix Beta) {
  int MaxIter = max_iter_LME;
  int Ndim = NumberDimensions;
  int NumIter = 0;    // Iterator counter.
  double norm_r = 10; // Value of the norm.
  Matrix p;           // Shape function vector.
  Matrix r;           // Gradient of log(Z).
  Matrix J;           // Hessian of log(Z).
  Matrix D_lambda;    // Increment of lambda.

  while (NumIter <= MaxIter) {

    // Get vector with the shape functions evaluated in the nodes.
    p = p__aLME__(l, lambda, Beta);

    // Get the gradient of log(Z) and its norm.
    r = r__aLME__(l, p);
    norm_r = norm__MatrixLib__(r, 2);

    // Check convergence.
    if (norm_r > TOL_wrapper_LME) {
      // Get the Hessian of log(Z).
      J = J__aLME__(l, p, r);

      if (rcond__TensorLib__(J.nV) < 1E-8) {
        fprintf(stderr, "%s %i : %s \n",
                "Error in lambda_Newton_Rapson__aLME__ for particle",
                Idx_particle, "The Hessian near to singular matrix!");
        exit(EXIT_FAILURE);
      }

      // Get the increment of lambda.
      D_lambda = solve__MatrixLib__(J, r);

      // Update the value of lambda.
      for (int i = 0; i < Ndim; i++) {
        lambda.nV[i] -= D_lambda.nV[i];
      }

      free__MatrixLib__(p);
      free__MatrixLib__(r);
      free__MatrixLib__(J);
      free__MatrixLib__(D_lambda);

      NumIter++;
    } else {
      free__MatrixLib__(r);
      free__MatrixLib__(p);
      break;
    }
  }

  if (NumIter >= MaxIter) {
    fprintf(
        stderr, "%s %i : %s (%i)\n",
        "Warning in lambda_Newton_Rapson__aLME__ for particle", Idx_particle,
        "No convergence reached in the maximum number of interations", MaxIter);
    fprintf(stderr, "%s : %e\n", "Total Error", norm_r);
    exit(EXIT_FAILURE);
  }
}

/****************************************************************************/

/**
 * @brief Function with computes the scalar values which contains the
 * restrictions related with the widh and information entropy of the shape
 * function.
 *
 * @param la Vector form node "a" to particle
 * @param lambda Lagrange multiplier
 * @param Beta Thermalization tensor
 * @return double
 */
static double fa__aLME__(Matrix la, Matrix lambda, Matrix Beta) {
  int Ndim = NumberDimensions;
  double la_x_lambda = 0;
  double Beta_x_la = 0;
  double Beta_norm_la = 0;
  double fa = 0;

  for (int i = 0; i < Ndim; i++) {

    Beta_x_la = 0;

    for (int j = 0; j < Ndim; j++) {
      Beta_x_la += Beta.nM[i][j] * la.nV[j];
    }

    Beta_norm_la += la.nV[i] * Beta_x_la;

    la_x_lambda += la.nV[i] * lambda.nV[i];
  }

  fa = -Beta_norm_la + la_x_lambda;

  return fa;
}

/****************************************************************************/

Matrix p__aLME__(Matrix l, Matrix lambda, Matrix Beta) {

  // Definition of some parameters
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double Z_m1 = 0;
  Matrix p = allocZ__MatrixLib__(1, N_a); // Shape function in the nodes
  Matrix la = memory_to_matrix__MatrixLib__(
      1, Ndim, NULL); // Vector form node "a" to particle.

  // Get Z and the numerator.
  for (int a = 0; a < N_a; a++) {
    la.nV = l.nM[a];
    p.nV[a] = exp(fa__aLME__(la, lambda, Beta));
    Z += p.nV[a];
  }

  // Divide by Z and get the final value of the shape function.
  Z_m1 = (double)1 / Z;
  for (int a = 0; a < N_a; a++) {
    p.nV[a] *= Z_m1;
  }

  return p;
}

/****************************************************************************/

/**
 * @brief Compute the value of the objetive function to compute the lagrange
 * multiplier.
 *
 * @param l Set than contanins vector form neighborhood nodes to particle.
 * @param lambda Lagrange multiplier.
 * @param Beta  Thermalization tensor.
 * @return double
 */
static double logZ__aLME__(Matrix l, Matrix lambda, Matrix Beta) {
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  double Z = 0;
  double logZ = 0;
  Matrix la = memory_to_matrix__MatrixLib__(1, Ndim, NULL);

  for (int a = 0; a < N_a; a++) {
    la.nV = l.nM[a];
    Z += exp(fa__aLME__(la, lambda, Beta));
  }

  logZ = log(Z);

  return logZ;
}

/****************************************************************************/

/**
 * @brief Compute the value of the gradient of the objetive function to compute
 * the lagrange multiplier.
 *
 * @param l Set than contanins vector form neighborhood nodes to particle.
 * @param p Value of the shape function in the nodes.
 * @return Matrix
 */
static Matrix r__aLME__(Matrix l, Matrix p) {
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix r = allocZ__MatrixLib__(Ndim, 1); // Gradient definition

  for (int i = 0; i < Ndim; i++) {
    for (int a = 0; a < N_a; a++) {
      r.nV[i] += p.nV[a] * l.nM[a][i];
    }
  }

  return r;
}

/****************************************************************************/

/**
 * @brief Compute the hessian matrix for the function logZ
 *
 * @param l Set than contanins vector form neighborhood nodes to particle.
 * @param p Set with the evaluation of the shape function in the neighborhood
 * @param r  Gradient of the function logZ_dLambda "r".
 * @return Matrix
 */
static Matrix J__aLME__(Matrix l, Matrix p, Matrix r) {
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix J = allocZ__MatrixLib__(Ndim, Ndim);

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {

      for (int a = 0; a < N_a; a++) {
        J.nM[i][j] += p.nV[a] * l.nM[a][i] * l.nM[a][j];
      }

      J.nM[i][j] -= r.nV[i] * r.nV[j];
    }
  }

  return J;
}

/****************************************************************************/

Matrix dp__aLME__(Matrix l, Matrix p) {
  int N_a = l.N_rows;
  int Ndim = NumberDimensions;
  Matrix dp = allocZ__MatrixLib__(N_a, Ndim);
  Matrix r;      // Gradient of log(Z)
  Matrix J;      // Hessian of log(Z)
  Matrix Jm1;    // Inverse of J
  Matrix Jm1_la; // Auxiliar vector
  Matrix la = memory_to_matrix__MatrixLib__(
      Ndim, 1, NULL); // Distance to the neighbour (x-x_a)

  // Get the Gradient and the Hessian of log(Z).
  r = r__aLME__(l, p);
  J = J__aLME__(l, p, r);

  // Inverse of the Hessian.
  Jm1 = inverse__MatrixLib__(J);

  // Fill the gradient for each node.
  for (int a = 0; a < N_a; a++) {
    la.nV = l.nM[a];
    Jm1_la = matrix_product__MatrixLib__(Jm1, la);

    for (int i = 0; i < Ndim; i++) {
      dp.nM[a][i] = -p.nV[a] * Jm1_la.nV[i];
    }

    free__MatrixLib__(Jm1_la);
  }

  free__MatrixLib__(r);
  free__MatrixLib__(J);
  free__MatrixLib__(Jm1);

  return dp;
}

/****************************************************************************/

int local_search__aLME__(Particle MPM_Mesh, Mesh FEM_Mesh)
/*
  Search the closest node to the particle based in its previous position.
*/
{
  int STATUS = EXIT_SUCCESS;
  int STATUS_p = EXIT_SUCCESS;
  int Ndim = NumberDimensions;
  unsigned Np = MPM_Mesh.NumGP;
  unsigned p;

  Matrix X_p;       // Particle position
  Matrix dis_p;     // Particle displacement
  Matrix Delta_Xip; // Distance from particles to the nodes
  Matrix lambda_p;  // Lagrange multiplier of the shape function
  Matrix DF_p;      // Particle deformation gradient
  Matrix Beta_p;
  Matrix Cut_off_Ellipsoid_p;

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

#pragma omp for private(p, X_p, Delta_Xip, lambda_p, DF_p, Beta_p,             \
                        Cut_off_Ellipsoid_p)
    for (p = 0; p < Np; p++) {

      //  Get the global coordinates of the particle
      X_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_GC.nM[p]);
      lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[p]);
      DF_p = memory_to_matrix__MatrixLib__(Ndim, Ndim, MPM_Mesh.Phi.DF.nM[p]);
      Beta_p = memory_to_matrix__MatrixLib__(Ndim, Ndim, MPM_Mesh.Beta.nM[p]);
      Cut_off_Ellipsoid_p = memory_to_matrix__MatrixLib__(
          Ndim, Ndim, MPM_Mesh.Cut_off_Ellipsoid.nM[p]);

      // Update the termalization matrix and cut-off ellipsoid.
      update_beta__aLME__(Beta_p, DF_p);
      update_cut_off_ellipsoid__aLME__(Cut_off_Ellipsoid_p, DF_p);

      //  Free previous list of tributary nodes to the particle
      free__SetLib__(&MPM_Mesh.ListNodes[p]);
      MPM_Mesh.ListNodes[p] = NULL;

      // Calculate the new connectivity with the previous value of beta.
      MPM_Mesh.ListNodes[p] = tributary__aLME__(p, X_p, Cut_off_Ellipsoid_p,
                                                MPM_Mesh.I0[p], FEM_Mesh);

      //  Calculate number of nodes
      MPM_Mesh.NumberNodes[p] = lenght__SetLib__(MPM_Mesh.ListNodes[p]);

      //  Generate nodal distance list
      Delta_Xip = compute_distance__MeshTools__(MPM_Mesh.ListNodes[p], X_p,
                                                FEM_Mesh.Coordinates);

      // Update lagrange multiplier with Newton-Rapson or with Nelder-Mead.
      update_lambda_Newton_Rapson__aLME__(p, Delta_Xip, lambda_p, Beta_p);

      // Free memory
      free__MatrixLib__(Delta_Xip);
      free(Cut_off_Ellipsoid_p.nM);
      free(Beta_p.nM);
      free(DF_p.nM);
    }
  }

  if (STATUS == EXIT_FAILURE) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/****************************************************************************/

/**
 * @brief Update the termalization matrix employing the increment of the
 * deformation gradient.
 *
 * @param Beta Termalization matrix.
 * @param Delta_F Increment of the deformation gradient.
 */
static void update_beta__aLME__(Matrix Beta, Matrix Delta_F) {
  int Ndim = NumberDimensions;
  Matrix Delta_F_m1 = inverse__MatrixLib__(Delta_F);
  Matrix updated_Beta = allocZ__MatrixLib__(Ndim, Ndim);

#if NumberDimensions == 2

  updated_Beta.nM[0][0] =
      Delta_F_m1.nM[0][0] * Beta.nM[0][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[0][0] * Beta.nM[0][1] * Delta_F_m1.nM[1][0] +
      Delta_F_m1.nM[1][0] * Beta.nM[1][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[1][0] * Beta.nM[1][1] * Delta_F_m1.nM[1][0];

  updated_Beta.nM[0][1] =
      Delta_F_m1.nM[0][0] * Beta.nM[0][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[0][0] * Beta.nM[0][1] * Delta_F_m1.nM[1][1] +
      Delta_F_m1.nM[1][0] * Beta.nM[1][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[1][0] * Beta.nM[1][1] * Delta_F_m1.nM[1][1];

  updated_Beta.nM[1][0] =
      Delta_F_m1.nM[0][1] * Beta.nM[0][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[0][1] * Beta.nM[0][1] * Delta_F_m1.nM[1][0] +
      Delta_F_m1.nM[1][1] * Beta.nM[1][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[1][1] * Beta.nM[1][1] * Delta_F_m1.nM[1][0];

  updated_Beta.nM[1][1] =
      Delta_F_m1.nM[0][1] * Beta.nM[0][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[0][1] * Beta.nM[0][1] * Delta_F_m1.nM[1][1] +
      Delta_F_m1.nM[1][1] * Beta.nM[1][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[1][1] * Beta.nM[1][1] * Delta_F_m1.nM[1][1];
#endif

#if NumberDimensions == 3
  fprintf(stderr, "%s : %s !!! \n", "Error in update_beta__aLME__()",
          "This operation it is not implemented for 3D");
  exit(EXIT_FAILURE);
#endif

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      Beta.nM[i][j] = updated_Beta.nM[i][j];
    }
  }

  free__MatrixLib__(Delta_F_m1);
  free__MatrixLib__(updated_Beta);
}

/****************************************************************************/

/**
 * @brief Update the shape of the cut-off region employing the increment of the
 * deformation gradient
 *
 * @param Cut_off_Ellipsoid Cut-off matrix.
 * @param Delta_F Increment of the deformation gradient.
 */
static void update_cut_off_ellipsoid__aLME__(Matrix Cut_off_Ellipsoid,
                                             Matrix Delta_F) {
  int Ndim = NumberDimensions;
  Matrix Delta_F_m1 = inverse__MatrixLib__(Delta_F);
  Matrix updated_Cut_off_Ellipsoid = allocZ__MatrixLib__(Ndim, Ndim);

#if NumberDimensions == 2
  updated_Cut_off_Ellipsoid.nM[0][0] =
      Delta_F_m1.nM[0][0] * Cut_off_Ellipsoid.nM[0][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[0][0] * Cut_off_Ellipsoid.nM[0][1] * Delta_F_m1.nM[1][0] +
      Delta_F_m1.nM[1][0] * Cut_off_Ellipsoid.nM[1][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[1][0] * Cut_off_Ellipsoid.nM[1][1] * Delta_F_m1.nM[1][0];

  updated_Cut_off_Ellipsoid.nM[0][1] =
      Delta_F_m1.nM[0][0] * Cut_off_Ellipsoid.nM[0][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[0][0] * Cut_off_Ellipsoid.nM[0][1] * Delta_F_m1.nM[1][1] +
      Delta_F_m1.nM[1][0] * Cut_off_Ellipsoid.nM[1][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[1][0] * Cut_off_Ellipsoid.nM[1][1] * Delta_F_m1.nM[1][1];

  updated_Cut_off_Ellipsoid.nM[1][0] =
      Delta_F_m1.nM[0][1] * Cut_off_Ellipsoid.nM[0][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[0][1] * Cut_off_Ellipsoid.nM[0][1] * Delta_F_m1.nM[1][0] +
      Delta_F_m1.nM[1][1] * Cut_off_Ellipsoid.nM[1][0] * Delta_F_m1.nM[0][0] +
      Delta_F_m1.nM[1][1] * Cut_off_Ellipsoid.nM[1][1] * Delta_F_m1.nM[1][0];

  updated_Cut_off_Ellipsoid.nM[1][1] =
      Delta_F_m1.nM[0][1] * Cut_off_Ellipsoid.nM[0][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[0][1] * Cut_off_Ellipsoid.nM[0][1] * Delta_F_m1.nM[1][1] +
      Delta_F_m1.nM[1][1] * Cut_off_Ellipsoid.nM[1][0] * Delta_F_m1.nM[0][1] +
      Delta_F_m1.nM[1][1] * Cut_off_Ellipsoid.nM[1][1] * Delta_F_m1.nM[1][1];
#endif

#if NumberDimensions == 3
  fprintf(stderr, "%s : %s !!! \n",
          "Error in update_cut_off_ellipsoid__aLME__()",
          "This operation it is not implemented for 3D");
  exit(EXIT_FAILURE);
#endif

  for (int i = 0; i < Ndim; i++) {
    for (int j = 0; j < Ndim; j++) {
      Cut_off_Ellipsoid.nM[i][j] = updated_Cut_off_Ellipsoid.nM[i][j];
    }
  }

  free__MatrixLib__(Delta_F_m1);
  free__MatrixLib__(updated_Cut_off_Ellipsoid);
}

/****************************************************************************/

/**
 * @brief  Update the list of tributary nodes for each particle
 *
 * @param Indx_p  Index of the particle.
 * @param X_p Position of the particle.
 * @param Cut_off_Ellipsoid Cut-off matrix.
 * @param I0 Index of the closest node to the particle.
 * @param FEM_Mesh Information of the background set of nodes.
 * @return ChainPtr
 */
static ChainPtr tributary__aLME__(int Indx_p, Matrix X_p,
                                  Matrix Cut_off_Ellipsoid, int I0,
                                  Mesh FEM_Mesh) {

  /* Define output */
  ChainPtr Triburary_Nodes = NULL;
  /* Number of dimensionws of the problem */
  int Ndim = NumberDimensions;

  Matrix Distance; /* Distance between node and GP */
  Matrix X_I = memory_to_matrix__MatrixLib__(Ndim, 1, NULL);

  ChainPtr Set_Nodes0 = NULL;
  int *Array_Nodes0;
  int NumNodes0;
  int Node0;

  // Counter
  int NumTributaryNodes = 0;

  // Get nodes close to the particle
  Set_Nodes0 = FEM_Mesh.NodalLocality[I0];
  NumNodes0 = FEM_Mesh.SizeNodalLocality[I0];
  Array_Nodes0 = set_to_memory__SetLib__(Set_Nodes0, NumNodes0);

  // Loop over the chain with the tributary nodes
  for (int i = 0; i < NumNodes0; i++) {

    Node0 = Array_Nodes0[i];

    if (FEM_Mesh.ActiveNode[Node0] == true) {
      // Assign to a pointer the coordinates of the nodes
      X_I.nV = FEM_Mesh.Coordinates.nM[Node0];

      // Get a vector from the GP to the node
      Distance = substraction__MatrixLib__(X_p, X_I);

      // If the node is near the GP push in the chain
      if (generalised_Euclidean_distance__MatrixLib__(
              Distance, Cut_off_Ellipsoid) <= 1.0) {
        push__SetLib__(&Triburary_Nodes, Node0);
        NumTributaryNodes++;
      }

      free__MatrixLib__(Distance);
    }
  }

  // If the Triburary_Nodes chain lenght is less than 3 assign al the node
  if (NumTributaryNodes < Ndim + 1) {
    fprintf(stderr, "%s %i : %s -> %i\n",
            "Warning in tributary__aLME__ for particle", Indx_p,
            "Insufficient nodal connectivity", NumTributaryNodes);
    exit(EXIT_FAILURE);
  }

  free(Array_Nodes0);

  return Triburary_Nodes;
}

/****************************************************************************/
