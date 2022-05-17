#include "Nodes/Nodes-Tools.h"

/**************************************************************/

void generate_contour_nodes(Mesh FEM_Mesh) {

  unsigned Nnodes = FEM_Mesh.NumNodesMesh;
  unsigned ContourNodes = 0;
  unsigned SizeNodalLocality_A = 0;
  unsigned NumActiveNodes_A = 0;
  unsigned B;
  ChainPtr NodalLocality_A;

  for (unsigned A = 0; A < Nnodes; A++) {
    if (FEM_Mesh.ActiveNode[A]) {
      SizeNodalLocality_A = FEM_Mesh.SizeNodalLocality_0[A];
      NodalLocality_A = FEM_Mesh.NodalLocality_0[A];
      NumActiveNodes_A = 0;

      while (NodalLocality_A != NULL) {

        B = NodalLocality_A->Idx;

        if (FEM_Mesh.ActiveNode[B]) {
          NumActiveNodes_A++;
        }

        NodalLocality_A = NodalLocality_A->next;
      }

      if (NumActiveNodes_A < 9) {
        FEM_Mesh.BoundaryNode[A] = true;
      } else {
        FEM_Mesh.BoundaryNode[A] = false;
      }

    } else {
      FEM_Mesh.BoundaryNode[A] = false;
    }
  }
}

/**************************************************************/

Mask get_active_nodes__MeshTools__(Mesh FEM_Mesh) {
  int Nnodes = FEM_Mesh.NumNodesMesh;
  int Nactivenodes = 0;
  int *Nodes2Mask = (int *)Allocate_ArrayZ(Nnodes, sizeof(int));
  Mask M;

  for (int A = 0; A < Nnodes; A++) {

    if (FEM_Mesh.ActiveNode[A]) {
      Nodes2Mask[A] = Nactivenodes;
      Nactivenodes++;
    } else {
      Nodes2Mask[A] = -1;
    }
  }

  M.Nactivenodes = Nactivenodes;
  M.Nodes2Mask = Nodes2Mask;

  return M;
}

/**************************************************************/

Mask get_active_dofs__MeshTools__(
  Mask ActiveNodes,
  Mesh FEM_Mesh, 
  int Step,
  int NumTimeStep) {

  /*
    Define auxilar variables
  */
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask * Ndof;
  int Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Id_BCC;        /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  /*
    Generate mask for the static condensation.
  */
  int Nactivenodes = 0;
  int *Nodes2Mask = (int *)Allocate_ArrayZ(Order, sizeof(int));
  Mask Free_and_Restricted_Dofs;

  /*
    Loop over the the boundaries to find the constrained dofs
  */
  for (int i = 0; i < Number_of_BCC; i++) {

    /*
      Get the number of nodes of this boundary
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    /*
      Get the number of dimensions where the BCC it is applied
    */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

    for (int j = 0; j < NumNodesBound; j++) {
      /*
        Get the index of the node
      */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

      /*
        If the boundary condition is under an active node
      */
      if (Id_BCC_mask != -1) {
        /*
          Loop over the dimensions of the boundary condition
        */
        for (int k = 0; k < NumDimBound; k++) {

          /*
            Apply only if the direction is active
          */
          if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + Step] == 1) {
            Id_BCC_mask_k = Id_BCC_mask * Ndof + k;
            Nodes2Mask[Id_BCC_mask_k] = -1;
          }
        }
      }
    }
  }

  /*
    Generate mask using the location of the constrained dofs
  */
  for (int A_i = 0; A_i < Order; A_i++) {
    if (Nodes2Mask[A_i] != -1) {
      Nodes2Mask[A_i] = Nactivenodes;
      Nactivenodes++;
    }
  }

  /*
    Assign information to the output
  */
  Free_and_Restricted_Dofs.Nactivenodes = Nactivenodes;
  Free_and_Restricted_Dofs.Nodes2Mask = Nodes2Mask;

  return Free_and_Restricted_Dofs;
}

/**************************************************************/

Mask generate_Mask_for_static_condensation_upw__MeshTools__(Mask ActiveNodes,
                                                            Mesh FEM_Mesh,
                                                            int Step,
                                                            int NumTimeStep) {

  /*
    Define auxilar variables
  */
  int Ndof = NumberDOF;
  int Nnodes_mask = ActiveNodes.Nactivenodes;
  int Order = Nnodes_mask * Ndof;
  int Number_of_BCC = FEM_Mesh.Bounds.NumBounds;
  int NumNodesBound; /* Number of nodes of the bound */
  int NumDimBound;   /* Number of dimensions */
  int Id_BCC;        /* Index of the node where we apply the BCC */
  int Id_BCC_mask;
  int Id_BCC_mask_k;

  /*
    Generate mask for the static condensation.
  */
  int Nactivenodes = 0;
  int *Nodes2Mask = (int *)Allocate_ArrayZ(Order, sizeof(int));
  Mask Free_and_Restricted_Dofs;

  /*
    Loop over the the boundaries to find the constrained dofs
  */
  for (int i = 0; i < Number_of_BCC; i++) {

    /*
      Get the number of nodes of this boundary
    */
    NumNodesBound = FEM_Mesh.Bounds.BCC_i[i].NumNodes;

    /*
      Get the number of dimensions where the BCC it is applied
    */
    NumDimBound = FEM_Mesh.Bounds.BCC_i[i].Dim;

    for (int j = 0; j < NumNodesBound; j++) {
      /*
        Get the index of the node
      */
      Id_BCC = FEM_Mesh.Bounds.BCC_i[i].Nodes[j];
      Id_BCC_mask = ActiveNodes.Nodes2Mask[Id_BCC];

      /*
        If the boundary condition is under an active node
      */
      if (Id_BCC_mask != -1) {
        /*
          Loop over the dimensions of the boundary condition
        */
        for (int k = 0; k < NumDimBound; k++) {

          /*
            Apply only if the direction is active
          */
          if (FEM_Mesh.Bounds.BCC_i[i].Dir[k * NumTimeStep + Step] == 1) {
            Id_BCC_mask_k = Id_BCC_mask + Nnodes_mask * k;
            Nodes2Mask[Id_BCC_mask_k] = -1;
          }
        }
      }
    }
  }

  /*
    Generate mask using the location of the constrained dofs
  */
  for (int A_i = 0; A_i < Order; A_i++) {
    if (Nodes2Mask[A_i] != -1) {
      Nodes2Mask[A_i] = Nactivenodes;
      Nactivenodes++;
    }
  }

  /*
    Assign information to the output
  */
  Free_and_Restricted_Dofs.Nactivenodes = Nactivenodes;
  Free_and_Restricted_Dofs.Nodes2Mask = Nodes2Mask;

  return Free_and_Restricted_Dofs;
}

/**************************************************************/

void get_set_scalar_field__MeshTools__(
  double * Field_Ap,
  const double * Field, 
  Element Nodes_p,
  Mask ActiveNodes)
/*
  This function performs two operations. First takes the nodal connectivity of
  the particle, and translate it to the mask numeration. Second, generate a
  Matrix with the nodal values. To help in the future computations. Nodal data
  is substracted in the shape (nodesxndofs).
 */
{
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = NumberDimensions;
  int Ap;
  int A_mask;

  for (unsigned A = 0; A < Nnodes; A++) {
  
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];
    
    Field_Ap[A] = Field[A_mask];
  }
}

/**************************************************************/

void get_set_vectorial_field__MeshTools__(
  double * Field_Ap,
  const double * Field, 
  Element Nodes_p,
  Mask ActiveNodes)
/*
  This function performs two operations. First takes the nodal connectivity of
  the particle, and translate it to the mask numeration. Second, generate a
  Matrix with the nodal values. To help in the future computations. Nodal data
  is substracted in the shape (nodesxndofs).
 */
{
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = NumberDimensions;
  int Ap;
  int A_mask;

  for (unsigned A = 0; A < Nnodes; A++) {
  
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];
    
    for (unsigned i = 0; i < Ndim; i++) {
      Field_Ap[A*Ndim + i] = Field[A_mask*Ndim + i];
    }
  }
}

/*********************************************************************/

void gradient_vector_field__MeshTools__(
  double * gradient_vector_field,
  const double * field,
  const double * gradient_Np,
  const ChainPtr ListNodes,
  int * STATUS){

  ChainPtr I = ListNodes;

  unsigned Ndim = NumberDimensions;
  unsigned I_global_Idx = I->Idx;
  unsigned I_local_Idx = 0;

  while (I != NULL) {
   
    for (unsigned i = 0; i < Ndim; i++) {
      for (unsigned j = 0; j < Ndim; j++) {
        gradient_vector_field[i*Ndim + j] += 
        field[I_global_Idx*Ndim + i] * 
        gradient_Np[I_local_Idx*Ndim + j];
      }
    }

    I = I->next;
    I_local_Idx++;
    I_global_Idx = I->Idx;
  }

}

/*********************************************************************/

Matrix get_U_set_field_upw__MeshTools__(Matrix Field_upw, Element Nodes_p,
                                        Mask ActiveNodes) {
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = NumberDimensions;
  Matrix Field_U_Ap = allocZ__MatrixLib__(Nnodes, Ndim);
  int Ap;
  int A_mask;

  for (int A = 0; A < Nnodes; A++) {
    /*
      Get the node in the mass matrix with the mask
    */
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];

    /*
      Get nodal field for particle p
    */
    for (int i = 0; i < Ndim; i++) {
      Field_U_Ap.nM[A][i] = Field_upw.nM[A_mask][i];
    }
  }

  return Field_U_Ap;
}

/*********************************************************************/

Matrix get_Pw_set_field_upw__MeshTools__(Matrix Field_upw, Element Nodes_p,
                                         Mask ActiveNodes) {
  int Nnodes = Nodes_p.NumberNodes;
  int Ndim = NumberDimensions;
  Matrix Field_pw_Ap = allocZ__MatrixLib__(Nnodes, 1);
  int Ap;
  int A_mask;

  for (int A = 0; A < Nnodes; A++) {
    /*
      Get the node in the mass matrix with the mask
    */
    Ap = Nodes_p.Connectivity[A];
    A_mask = ActiveNodes.Nodes2Mask[Ap];

    /*
      Get nodal field for particle p
    */
    Field_pw_Ap.nV[A] = Field_upw.nM[A_mask][Ndim];
  }

  return Field_pw_Ap;
}

/*********************************************************************/

Matrix get_nodes_coordinates__MeshTools__(ChainPtr Element_p,
                                          Matrix Coordinates)
/*
  Get the matrix with the coordinates of an element
*/
{

  int Ndim = NumberDimensions;
  int NumVertex = lenght__SetLib__(Element_p);
  Matrix Element_Coordinates = allocZ__MatrixLib__(NumVertex, Ndim);
  ChainPtr Idx = NULL;
  int I_Idx = 0;

  Idx = Element_p;

  while (Idx != NULL) {

    /* Fill the elelemtn coordinates */
    for (int l = 0; l < Ndim; l++) {
      Element_Coordinates.nM[I_Idx][l] = Coordinates.nM[Idx->Idx][l];
    }

    /* Cycle */
    Idx = Idx->next;
    I_Idx++;
  }

  return Element_Coordinates;
}

/*********************************************************************/

double point_distance__MeshTools__(Matrix End, Matrix Init)
/*
  Distance between two points.
*/
{
  int Ndim = NumberDimensions;

  if ((End.N_rows != Init.N_rows) || (End.N_cols != Init.N_cols)) {
    printf("%s : %s \n", "Error in point_distance__MeshTools__",
           "Inputs arrays are not equal");
    exit(EXIT_FAILURE);
  }

  double DIST = 0;

  for (int i = 0; i < Ndim; i++) {

    DIST += pow(End.nV[i] - Init.nV[i], 2);
  }

  DIST = sqrt(DIST);

  return DIST;
}

/*********************************************************************/

Matrix compute_distance__MeshTools__(ChainPtr Set, Matrix X0,
                                     Matrix Coordinates) {
  int Ndim = NumberDimensions;
  int SizeSet = lenght__SetLib__(Set);
  int I = 0;

  /* Allocate output */
  Matrix Set_Coordinates = allocZ__MatrixLib__(SizeSet, Ndim);

  /* Loop in the set */
  ChainPtr Aux_Set = Set;
  while (Aux_Set != NULL) {
    /* Get coordinates local coodinates of each node in the set */
    for (int i = 0; i < Ndim; i++) {
      Set_Coordinates.nM[I][i] = X0.nV[i] - Coordinates.nM[Aux_Set->Idx][i];
    }
    /* Update index */
    I++;
    Aux_Set = Aux_Set->next;
  }

  return Set_Coordinates;
}

/*********************************************************************/

Matrix get_set_field_old__MeshTools__(Matrix Nodal_Field, Element GP_Element) {

  int *Element_Connectivity = GP_Element.Connectivity;
  int NumNodes = GP_Element.NumberNodes;
  Matrix Element_Field;
  int Ndim = NumberDimensions;
  int Ie;

  /* Allocate a matrix to store the nodal quatities in the element */
  Element_Field = alloc__MatrixLib__(NumNodes, Ndim);

  /* Loop over the nodes of the element */
  for (int I = 0; I < NumNodes; I++) {
    /* Get the node of the element */
    Ie = Element_Connectivity[I];
    /* Fill each dimension of the nodal quantitie */
    for (int i = 0; i < Ndim; i++) {
      Element_Field.nM[I][i] = Nodal_Field.nM[Ie][i];
    }
  }

  return Element_Field;
}

/*********************************************************************/

int get_closest_node__MeshTools__(Matrix X_p, ChainPtr Nodes,
                                  Matrix Coordinates)
/*
  Ordenate recursively and array with distances and get a chain
  with the positions in orden
*/
{
  /* Number of dimensions */
  int Ndim = NumberDimensions;
  /* Index of the current node */
  int I;
  /* Corrdinates of the current node */
  Matrix X_I;
  /* Distance from the particle to the node */
  double Distance_I;
  /* Interator pointer */
  ChainPtr Node_I = NULL;
  /* Get the value of the maximum distance from the particle to the pointer */
  double DistMin;
  int I_DistMin;

  /* Initialize interator with the first node */
  Node_I = Nodes;

  /* Get the index of the first node */
  I = Node_I->Idx;

  /* Get the coordinates of the first node */
  X_I = memory_to_matrix__MatrixLib__(Ndim, 1, Coordinates.nM[I]);

  /* Get the distance from the node to the particle */
  Distance_I = point_distance__MeshTools__(X_p, X_I);

  /* Get the distance from the node to the particle */
  DistMin = Distance_I;
  I_DistMin = Node_I->Idx;

  /* Search in the reamaining nodes */
  Node_I = Node_I->next;

  while (Node_I != NULL) {

    /* Get the index of the node */
    I = Node_I->Idx;

    /* Get the coordinates of the node */
    X_I = memory_to_matrix__MatrixLib__(Ndim, 1, Coordinates.nM[I]);

    /* Get the distance from the node to the particle */
    Distance_I = point_distance__MeshTools__(X_p, X_I);

    /* Get the max distance of the matrix */
    if (Distance_I < DistMin) {
      DistMin = Distance_I;
      I_DistMin = Node_I->Idx;
    }

    /* Continue iterating */
    Node_I = Node_I->next;
  }

  return I_DistMin;
}

/*********************************************************************/

double interpolate_scalar_magnitude__MeshTools__(Matrix A, Matrix N_p) {
  int Nnodes_p = A.N_rows;
  double A_p = 0;

  for (int I = 0; I < Nnodes_p; I++) {
    A_p += A.nV[I] * N_p.nV[I];
  }

  return A_p;
}

/**************************************************************/

Tensor interpolate_vectorial_magnitude__MeshTools__(Matrix Nodal_Variable,
                                                    Matrix ShapeFunction_p)
/*

*/
{

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = Nodal_Variable.N_rows;
  double ShapeFunction_pI;
  Tensor Variable_I;
  Tensor Variable_p = alloc__TensorLib__(1);

  for (int I = 0; I < Nnodes_p; I++) {

    /*
      Assign from matrix
    */
    Variable_I = memory_to_tensor__TensorLib__(Nodal_Variable.nM[I], 1);
    ShapeFunction_pI = ShapeFunction_p.nV[I];

    /*
      Compute nodal contribution
    */
    for (int i = 0; i < Ndim; i++) {
      Variable_p.n[i] += ShapeFunction_pI * Variable_I.n[i];
    }
  }

  return Variable_p;
}

/**************************************************************/

Tensor interpolate_vectorial_magnitude_gradient__MeshTools__(Matrix Variable_p,
                                                             Matrix gradient_p)
/*

*/
{

  /* Variable definition */
  int Ndim = NumberDimensions;
  int Nnodes_p = Variable_p.N_rows;
  Tensor Variable_pI;
  Tensor gradient_I;
  Tensor Variable__o__gradient_Ip;
  Tensor grad_variable_p = alloc__TensorLib__(2);

  for (int I = 0; I < Nnodes_p; I++) {

    /*
      Assign from matrix
    */
    Variable_pI = memory_to_tensor__TensorLib__(Variable_p.nM[I], 1);
    gradient_I = memory_to_tensor__TensorLib__(gradient_p.nM[I], 1);

    Variable__o__gradient_Ip =
        dyadic_Product__TensorLib__(Variable_pI, gradient_I);

    /*
      Compute nodal contribution
    */
    for (int i = 0; i < Ndim; i++) {
      for (int j = 0; j < Ndim; j++) {
        grad_variable_p.N[i][j] += Variable__o__gradient_Ip.N[i][j];
      }
    }

    free__TensorLib__(Variable__o__gradient_Ip);
  }

  return grad_variable_p;
}

/*********************************************************************/
