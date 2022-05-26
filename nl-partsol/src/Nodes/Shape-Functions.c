// clang-format off
#include "Nodes/Shape-Functions.h"
#include "Globals.h"
// clang-format on

/*********************************************************************/

void initialise_shapefun__MeshTools__(Particle MPM_Mesh, Mesh FEM_Mesh) {

  if (strcmp(ShapeFunctionGP, "FEM") == 0) {
    if (strcmp(FEM_Mesh.TypeElem, "Triangle") == 0) {
      initialize__T3__(MPM_Mesh, FEM_Mesh);
    } else if (strcmp(FEM_Mesh.TypeElem, "Quadrilateral") == 0) {
      initialize__Q4__(MPM_Mesh, FEM_Mesh);
    } else if (strcmp(FEM_Mesh.TypeElem, "Tetrahedra") == 0) {
      initialize__T4__(MPM_Mesh, FEM_Mesh);
    } else if (strcmp(FEM_Mesh.TypeElem, "Hexahedra") == 0) {
      initialize__H8__(MPM_Mesh, FEM_Mesh);
    }
  } else if (strcmp(ShapeFunctionGP, "uGIMP") == 0) {
    initialize__GIMP__(MPM_Mesh, FEM_Mesh);
  } else if (strcmp(ShapeFunctionGP, "LME") == 0) {
    initialize__LME__(MPM_Mesh, FEM_Mesh);
  } else if (strcmp(ShapeFunctionGP, "aLME") == 0) {
    initialize__aLME__(MPM_Mesh, FEM_Mesh);
  }
}

/*********************************************************************/

void local_search__MeshTools__(Particle MPM_Mesh, Mesh FEM_Mesh)
/*
  Search the closest node to the particle based in its previous position.
*/
{

  //! Set to zero the active/non-active node, and the GPs in each element
  for (unsigned i = 0; i < FEM_Mesh.NumNodesMesh; i++) {

    FEM_Mesh.ActiveNode[i] = false;

    if ((Driver_EigenErosion == true) || (Driver_EigenSoftening == true)) {
      free__SetLib__(&FEM_Mesh.List_Particles_Node[i]);
    }
  }

  if (FEM_Mesh.Locking_Control_Fbar) {
    for (unsigned i = 0; i < FEM_Mesh.Num_Patch_Mesh; i++) {
      FEM_Mesh.Vol_Patch_n[i] = 0.0;
      FEM_Mesh.Vol_Patch_n1[i] = 0.0;
    }
  }

  //! Local search
  if (strcmp(ShapeFunctionGP, "FEM") == 0) {
    if ((strcmp(FEM_Mesh.TypeElem, "Triangle") == 0) &&
        (FEM_Mesh.NumNodesElem[0] == 3)) {
      local_search__T3__(MPM_Mesh, FEM_Mesh);
    } else if ((strcmp(FEM_Mesh.TypeElem, "Triangle") == 0) &&
               (FEM_Mesh.NumNodesElem[0] == 6)) {
      local_search__T3__(MPM_Mesh, FEM_Mesh);
    } else if (strcmp(FEM_Mesh.TypeElem, "Quadrilateral") == 0) {
      local_search__Q4__(MPM_Mesh, FEM_Mesh);
    } else if (strcmp(FEM_Mesh.TypeElem, "Tetrahedra") == 0) {
      local_search__T4__(MPM_Mesh, FEM_Mesh);
    } else if (strcmp(FEM_Mesh.TypeElem, "Hexahedra") == 0) {
      local_search__H8__(MPM_Mesh, FEM_Mesh);
    }
  } else if (strcmp(ShapeFunctionGP, "uGIMP") == 0) {
    exit(1);
  } else if (strcmp(ShapeFunctionGP, "LME") == 0) {
    local_search__LME__(MPM_Mesh, FEM_Mesh);
  } else if (strcmp(ShapeFunctionGP, "aLME") == 0) {
    local_search__aLME__(MPM_Mesh, FEM_Mesh);
  }
}

/*********************************************************************/

Matrix compute_N__MeshTools__(Element GP_Element, Particle MPM_Mesh,
                              Mesh FEM_Mesh) {

  /*
    Auxiliar variables
  */
  int Ndim = NumberDimensions;
  int i_GP = GP_Element.i_GP; // Current intex of the particle
  int I_p; // Index of a node in the neiborhood of the particle
  int NumNodes_p =
      GP_Element
          .NumberNodes; // Number of nodes in the neiborhood of the particle
  int *Connectivity_p =
      GP_Element
          .Connectivity; // List of nodes in the neiborhood of the particle

  /*
    Matrix with the nodal shape functions
  */
  Matrix ShapeFunction_p;

  if (strcmp(ShapeFunctionGP, "FEM") == 0) {
    Matrix xi_p; // Element coordinates of the Gauss-Point

    /*
      Get the element coordinates of the GP
    */
    xi_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_EC.nM[i_GP]);

    /*
      Evaluate the shape function
    */
    ShapeFunction_p = FEM_Mesh.N_ref(xi_p);

  }

  else if (strcmp(ShapeFunctionGP, "uGIMP") == 0) {

    Matrix lp;   // Particle voxel
    Matrix l_Ip; // Distance from GP to the nodes

    /*
      Get the distance of the particle to the nodes
    */
    l_Ip = alloc__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        l_Ip.nM[k][l] =
            MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get the GP voxel
    */
    lp.nV = MPM_Mesh.lp.nM[i_GP];

    /*
      Evaluate the shape function
    */
    ShapeFunction_p = N__GIMP__(l_Ip, lp, FEM_Mesh.DeltaX);

    /*
      Free memory
    */
    free__MatrixLib__(l_Ip);
  }

  else if (strcmp(ShapeFunctionGP, "LME") == 0) {
    Matrix l_Ip;     // Distance from particle to the nodes
    Matrix lambda_p; // Lagrange multipliers
    double Beta_p;   // Thermalization parameter

    /*
      Get the distance of the particle to the nodes
    */
    l_Ip = alloc__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        l_Ip.nM[k][l] =
            MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get lambda and beta
    */
    lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[i_GP]);
    Beta_p = MPM_Mesh.Beta.nV[i_GP];

    /*
      Evaluate the shape function
    */
    ShapeFunction_p = p__LME__(l_Ip, lambda_p, Beta_p);

    /*
      Free memory
    */
    free__MatrixLib__(l_Ip);
  } else if (strcmp(ShapeFunctionGP, "aLME") == 0) {
    Matrix l_Ip;     // Distance from particle to the nodes
    Matrix lambda_p; // Lagrange multipliers
    Matrix Beta_p;   // Thermalization parameter

    /*
      Get the distance of the particle to the nodes
    */
    l_Ip = alloc__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        l_Ip.nM[k][l] =
            MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get lambda and beta
    */
    lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[i_GP]);
    Beta_p = memory_to_matrix__MatrixLib__(Ndim, Ndim, MPM_Mesh.Beta.nM[i_GP]);

    /*
      Evaluate the shape function
    */
    ShapeFunction_p = p__aLME__(l_Ip, lambda_p, Beta_p);

    /*
      Free memory
    */
    free__MatrixLib__(l_Ip);
    free(Beta_p.nM);
  }

  else {
    printf("%s : %s %s %s \n", "Error in Get_Operator()", "The shape-function ",
           ShapeFunctionGP, "is not implemented");
    exit(EXIT_FAILURE);
  }

  return ShapeFunction_p;
}

/*********************************************************************/

Matrix compute_dN__MeshTools__(Element GP_Element, Particle MPM_Mesh,
                               Mesh FEM_Mesh) {
  int Ndim = NumberDimensions;
  int i_GP = GP_Element.i_GP; // Current intex of the particle
  int I_p; // Index of a node in the neiborhood of the particle
  int NumNodes_p =
      GP_Element
          .NumberNodes; // Number of nodes in the neiborhood of the particle
  int *Connectivity_p =
      GP_Element
          .Connectivity; // List of nodes in the neiborhood of the particle

  /*
    Matrix with the nodal derivatives
  */
  Matrix Gradient_p;

  if (strcmp(ShapeFunctionGP, "FEM") == 0) {
    Matrix xi_p; // Element coordinates of the Gauss-Point
    Matrix X_I;  // Coordinates of the set of nodes in the neiborhood of the
                 // particle

    /*
      Fill the poligon woth the nodal coordinates of the current element
    */
    X_I = allocZ__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        X_I.nM[k][l] = FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get the element coordinates of the GP
    */
    xi_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.Phi.x_EC.nM[i_GP]);

    /*
      Evaluate the shape function gradient
    */
    Gradient_p = FEM_Mesh.dNdX(xi_p, X_I);

    /* Free memory */
    free__MatrixLib__(X_I);
  }

  else if (strcmp(ShapeFunctionGP, "uGIMP") == 0) {
    Matrix lp;   // Particle voxel
    Matrix l_Ip; // Distance from GP to the nodes
    /*
      Get the distance of the particle to the nodes
    */
    l_Ip = alloc__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        l_Ip.nM[k][l] =
            MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get the GP voxel
    */
    lp.nV = MPM_Mesh.lp.nM[i_GP];

    /*
      Evaluate the shape function gradient
    */
    Gradient_p = dN__GIMP__(l_Ip, lp, FEM_Mesh.DeltaX);

    /*
      Free memory
    */
    free__MatrixLib__(l_Ip);
  }

  else if (strcmp(ShapeFunctionGP, "LME") == 0) {
    Matrix l_Ip;            // Distance from GP to the nodes
    Matrix ShapeFunction_p; // Matrix with the nodal shape functions
    Matrix lambda_p;        // Lagrange multipliers
    double Beta_p;          // Thermalization parameter

    /*
      Get the distance of the particle to the nodes
    */
    l_Ip = alloc__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        l_Ip.nM[k][l] =
            MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get lambda and beta
    */
    lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[i_GP]);
    Beta_p = MPM_Mesh.Beta.nV[i_GP];

    /*
      Evaluate the shape function gradient
    */
    ShapeFunction_p = p__LME__(l_Ip, lambda_p, Beta_p);
    Gradient_p = dp__LME__(l_Ip, ShapeFunction_p);

    /*
      Free memory
    */
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(l_Ip);
  }

  else if (strcmp(ShapeFunctionGP, "aLME") == 0) {
    Matrix l_Ip;            // Distance from GP to the nodes
    Matrix ShapeFunction_p; // Matrix with the nodal shape functions
    Matrix lambda_p;        // Lagrange multipliers
    Matrix Beta_p;          // Thermalization parameter

    /*
      Get the distance of the particle to the nodes
    */
    l_Ip = alloc__MatrixLib__(NumNodes_p, Ndim);
    for (int k = 0; k < NumNodes_p; k++) {
      I_p = Connectivity_p[k];
      for (int l = 0; l < Ndim; l++) {
        l_Ip.nM[k][l] =
            MPM_Mesh.Phi.x_GC.nM[i_GP][l] - FEM_Mesh.Coordinates.nM[I_p][l];
      }
    }

    /*
      Get lambda and beta
    */
    lambda_p = memory_to_matrix__MatrixLib__(Ndim, 1, MPM_Mesh.lambda.nM[i_GP]);
    Beta_p = memory_to_matrix__MatrixLib__(Ndim, Ndim, MPM_Mesh.Beta.nM[i_GP]);

    /*
      Evaluate the shape function gradient
    */
    ShapeFunction_p = p__aLME__(l_Ip, lambda_p, Beta_p);
    Gradient_p = dp__aLME__(l_Ip, ShapeFunction_p);

    /*
      Free memory
    */
    free__MatrixLib__(ShapeFunction_p);
    free__MatrixLib__(l_Ip);
    free(Beta_p.nM);
  }

  else {
    printf("%s : %s %s %s \n", "Error in Get_Operator()", "The shape-function ",
           ShapeFunctionGP, "is not implemented");
    exit(EXIT_FAILURE);
  }

  return Gradient_p;
}

/*********************************************************************/

double *push_forward_dN__MeshTools__(const double *Gradient_n_p,
                                     const double *d_phi, unsigned NumNodes,
                                     int *STATUS) {

  unsigned Ndim = NumberDimensions;
  *STATUS = EXIT_SUCCESS;

  double *Gradient_n1_p = (double *)calloc(NumNodes * Ndim, __SIZEOF_DOUBLE__);
  if (Gradient_n1_p == NULL) {
    fprintf(stderr, "" RED "Error in calloc()" RESET " \n");
    *STATUS = EXIT_FAILURE;
    return Gradient_n1_p;
  }

#if NumberDimensions == 2
  double d_phi_mT[4];
#else
  double d_phi_mT[9];
#endif

  // Compute the adjunt of the incremental deformation gradient
  *STATUS = compute_adjunt__TensorLib__(d_phi_mT, d_phi);
  if (*STATUS == EXIT_FAILURE) {
    fprintf(stderr,
            "" RED "Error in compute_adjunt__TensorLib__()" RESET " \n");
    return Gradient_n1_p;
  }

  // Push-forward the shape function gradient using
  // the incremental deformation gradient
  for (unsigned A = 0; A < NumNodes; A++) {
    for (unsigned i = 0; i < Ndim; i++) {

      Gradient_n1_p[A * Ndim + i] = 0.0;

      for (unsigned j = 0; j < Ndim; j++) {
        Gradient_n1_p[A * Ndim + i] +=
            d_phi_mT[i * Ndim + j] * Gradient_n_p[A * Ndim + j];
      }
    }
  }

  return Gradient_n1_p;
}

/*********************************************************************/